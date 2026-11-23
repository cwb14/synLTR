#!/usr/bin/env python3
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor
import os
import glob
import itertools
import time
from typing import List, Set, Tuple, Dict

# -------------------------
# Helpers
# -------------------------

BLAST_DB_SUFFIXES = [".cds.njs", ".cds.ntf", ".cds.nto", ".cds.not", ".cds.ndb", ".cds.nhr", ".cds.nin", ".cds.nsq"]
LAST_DB_SUFFIXES  = [".des", ".sds", ".tis", ".ssp", ".bck", ".suf", ".prj"]
PROT_DB_SUFFIXES  = [".pep.dmnd"]

OUTPUT_CHOICES = ["pdf", "lifted", "anchors", "filtered", "last", "jcvi_list", "db", "all", "none"]
PAIR_OUTPUT_KEYS = ("pdf", "lifted", "anchors", "filtered", "last")
MAX_RETRIES = 2  # dumb retry count
RETRY_SLEEP_SEC = 2  # small backoff between retries

def ensure_list_file(list_file_path: str) -> None:
    """Create jcvi_list.txt if missing using *.bed or *.pep prefixes."""
    if os.path.isfile(list_file_path):
        return

    bed_files = glob.glob('*.bed')
    pep_files = glob.glob('*.pep')

    if bed_files:
        raw_files = bed_files
    elif pep_files:
        raw_files = pep_files
    else:
        raise FileNotFoundError(
            "No .bed or .pep files found in the current directory. "
            "Cannot generate jcvi_list.txt automatically."
        )

    accessions = sorted({os.path.splitext(os.path.basename(f))[0] for f in raw_files})
    pairs = itertools.combinations(accessions, 2)
    with open(list_file_path, 'w') as outf:
        count_pairs = 0
        for a, b in pairs:
            outf.write(f"{a}\t{b}\n")
            count_pairs += 1

    print(f"Generated '{list_file_path}' with {len(accessions)} accessions ({count_pairs} pairs).")

def read_pairs(list_path: str) -> List[Tuple[str, str]]:
    with open(list_path, 'r') as infile:
        lines = [ln.strip() for ln in infile if ln.strip()]
    pairs: List[Tuple[str, str]] = []
    for ln in lines:
        parts = ln.split()
        if len(parts) != 2:
            raise ValueError(f"Invalid line in {list_path}: {ln!r}")
        pairs.append((parts[0], parts[1]))
    return pairs

def build_expected_outputs_for_pair(col1: str, col2: str) -> Dict[str, str]:
    """Return mapping for the 5 common outputs tied to a pair (no dbs here)."""
    return {
        "pdf":      f"{col1}.{col2}.pdf",
        "lifted":   f"{col1}.{col2}.lifted.anchors",
        "anchors":  f"{col1}.{col2}.anchors",
        "filtered": f"{col1}.{col2}.last.filtered",
        "last":     f"{col1}.{col2}.last",
    }

def build_expected_db_files_for_accession(acc: str, mode: str) -> List[str]:
    """Return db files for a given accession under a given mode ('prot', 'blast', 'last')."""
    if mode == "prot":
        return [f"{acc}{suf}" for suf in PROT_DB_SUFFIXES]
    elif mode == "blast":
        return [f"{acc}{suf}" for suf in BLAST_DB_SUFFIXES]
    else:  # last (default)
        return [f"{acc}{suf}" for suf in LAST_DB_SUFFIXES]

def safe_unlink(path: str, dry_run: bool = False) -> None:
    if not os.path.exists(path):
        return
    if dry_run:
        print(f"[dry-run] would delete: {path}")
        return
    try:
        os.remove(path)
        print(f"deleted: {path}")
    except Exception as e:
        print(f"warning: failed to delete {path}: {e}")

def file_nonempty(path: str) -> bool:
    return os.path.isfile(path) and os.path.getsize(path) > 0

def choose_check_keys(args_keep: List[str] or None) -> Set[str]:
    """
    Decide which per-pair outputs we should verify to judge success.
    Logic:
      - if --keep is None: check ALL standard per-pair outputs
      - if --keep includes 'all': check ALL standard per-pair outputs
      - otherwise: check those per-pair outputs named in --keep
      - if that set is empty (e.g. user only kept 'db' or 'jcvi_list' or 'none'),
        fall back to checking 'anchors' as a reasonable core artifact.
    """
    if args_keep is None:
        keys = set(PAIR_OUTPUT_KEYS)
    elif "all" in args_keep:
        keys = set(PAIR_OUTPUT_KEYS)
    else:
        keys = {k for k in args_keep if k in PAIR_OUTPUT_KEYS}
        if not keys:
            # fallback so the retry heuristic still does something useful
            keys = {"anchors"}
    return keys

def outputs_ok_for_pair(col1: str, col2: str, check_keys: Set[str]) -> bool:
    expected = build_expected_outputs_for_pair(col1, col2)
    for key in check_keys:
        path = expected.get(key)
        if not path:
            continue
        if not file_nonempty(path):
            return False
    return True

# -------------------------
# Main
# -------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Run jcvi.compara.catalog ortholog in parallel and optionally clean outputs.'
    )
    parser.add_argument('-p', metavar='N', type=int, default=1,
                        help='Number of processes to run in parallel (default: 1)')
    parser.add_argument('--cpus', metavar='N', type=int, default=1,
                        help='Number of CPUs to use (default: 1)')
    parser.add_argument('--blast', action='store_true',
                        help='Use BLAST (adds --align_soft blast)')
    parser.add_argument('--prot', action='store_true',
                        help='Use DIAMOND protein (adds --align_soft diamond_blastp --dbtype prot)')
    parser.add_argument('--cscore', type=float, default=0.99,
                        help='Pass-through value for --cscore (default: 0.99)')
    parser.add_argument('--keep', nargs='+', choices=OUTPUT_CHOICES,
                        help=('Which outputs to KEEP after runs. Choices: '
                              'pdf lifted anchors filtered last jcvi_list db all none. '
                              'If omitted, keeps everything. '
                              'If provided (and not "all"), non-selected known artifacts will be deleted.'))
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be deleted by --keep without deleting.')
    args = parser.parse_args()

    # Validate mutually exclusive modes
    if args.blast and args.prot:
        raise SystemExit("Error: --blast and --prot are mutually exclusive.")

    # Ensure min values
    args.cpus = max(1, args.cpus)
    args.p = max(1, args.p)

    # Determine mode
    mode = "prot" if args.prot else ("blast" if args.blast else "last")

    # Paths to list files
    list_file_path = 'jcvi_list.txt'
    missing_list_file_path = 'jcvi_list_missing.txt'

    # Ensure list file exists or generate it
    ensure_list_file(list_file_path)

    # Which file to use
    file_to_use = missing_list_file_path if os.path.isfile(missing_list_file_path) else list_file_path

    # Read pairs
    pairs = read_pairs(file_to_use)

    # Decide which outputs we will verify for success
    check_keys = choose_check_keys(args.keep)

    # Function to run a command for a given line, with dumb retry on missing/empty outputs
    def run_command(pair: Tuple[str, str]) -> None:
        col1, col2 = pair
        expected_map = build_expected_outputs_for_pair(col1, col2)

        def _invoke():
            command = [
                'python', '-m', 'jcvi.compara.catalog', 'ortholog',
                col1, col2, '--no_strip_names', f'--cpus={args.cpus}',
                '--notex', f'--cscore={args.cscore}'
            ]
            if args.blast:
                command.extend(['--align_soft', 'blast'])
            elif args.prot:
                command.extend(['--align_soft', 'diamond_blastp', '--dbtype', 'prot'])
            # Run
            subprocess.run(command, check=False)

        # Attempt run + verification, with retries
        attempt = 0
        while True:
            attempt += 1
            if attempt > 1:
                # Small courtesy pause before retrying
                time.sleep(RETRY_SLEEP_SEC)

            _invoke()

            if outputs_ok_for_pair(col1, col2, check_keys):
                # Success
                return

            if attempt <= MAX_RETRIES:
                # Best-effort cleanup of any zero-byte files among the check set before retry
                for key in check_keys:
                    path = expected_map.get(key)
                    if path and os.path.isfile(path) and os.path.getsize(path) == 0:
                        try:
                            os.remove(path)
                        except Exception:
                            pass
                print(f"[retry] Detected missing/empty outputs for {col1},{col2} "
                      f"(checked: {', '.join(sorted(check_keys))}); retry {attempt}/{MAX_RETRIES}.")
                continue
            else:
                # Give up after MAX_RETRIES
                missing = []
                for key in check_keys:
                    path = expected_map.get(key)
                    if not path:
                        continue
                    if not os.path.isfile(path):
                        missing.append(f"{key} (missing)")
                    elif os.path.getsize(path) == 0:
                        missing.append(f"{key} (empty)")
                status = "; ".join(missing) if missing else "unknown failure state"
                print(f"[warn] Giving up on {col1},{col2} after {MAX_RETRIES} retries: {status}.")
                return

    # Run in parallel
    with ThreadPoolExecutor(max_workers=args.p) as executor:
        list(executor.map(run_command, pairs))

    # -------------------------
    # Cleanup phase (optional)
    # -------------------------
    if args.keep is None:
        # No cleanup requested: keep everything
        return

    keep_set: Set[str] = set()
    # Expand 'all' / 'none'
    keep_all = "all" in args.keep
    keep_none = ("none" in args.keep) and (len(args.keep) == 1)

    # Always consider only known artifacts.
    known_artifacts: Set[str] = set()

    # Common per-pair artifacts
    for col1, col2 in pairs:
        expected = build_expected_outputs_for_pair(col1, col2)
        known_artifacts.update(expected.values())
        if keep_all or keep_none:
            # We'll add to keep_set later if keep_all; if keep_none, we don't add now.
            pass
        else:
            if "pdf" in args.keep:
                keep_set.add(expected["pdf"])
            if "lifted" in args.keep:
                keep_set.add(expected["lifted"])
            if "anchors" in args.keep:
                keep_set.add(expected["anchors"])
            if "filtered" in args.keep:
                keep_set.add(expected["filtered"])
            if "last" in args.keep:
                keep_set.add(expected["last"])

    # jcvi_list
    known_artifacts.add(list_file_path)
    if keep_all or ("jcvi_list" in args.keep):
        keep_set.add(list_file_path)

    # DB artifacts: by convention jcvi builds DBs for the second accession in each pair
    second_accessions = sorted({b for (_, b) in pairs})
    for acc in second_accessions:
        for db_file in build_expected_db_files_for_accession(acc, mode):
            known_artifacts.add(db_file)
            if keep_all or ("db" in args.keep):
                keep_set.add(db_file)

    # If keep_all: keep everything we recognize
    if keep_all:
        keep_set = set(known_artifacts)

    # If keep_none: keep nothing (but still only delete known artifacts)
    if keep_none:
        keep_set = set()

    # Compute deletions = known - keep
    to_delete = sorted(known_artifacts - keep_set)

    if not to_delete:
        if args.dry_run:
            print("[dry-run] nothing to delete based on the requested --keep set.")
        return

    print(f"{'[dry-run] ' if args.dry_run else ''}cleanup: {len(to_delete)} file(s) to remove...")
    for path in to_delete:
        safe_unlink(path, dry_run=args.dry_run)

if __name__ == "__main__":
    main()
