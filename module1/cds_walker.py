#!/usr/bin/env python3
import sys
import argparse
import logging
from typing import List, Tuple, Optional
from multiprocessing import Pool
from Bio import SeqIO
from Bio.Seq import Seq

# A companion script to 'liftover.py'.
# Takes miniprot CDS and proteins and removes 1-2bp nt insertions that cause frameshift.
# Adds a cautious fallback for 3–4 bp insertions that also keeps the *next* AA in-frame.
# Outputs a modified version of the CDS that directly translates to the protein.
# Useful for KaKs.

def _aa_of(codon: str, table: int) -> Optional[str]:
    """Translate a 3bp codon->AA (single letter). Returns None if codon is incomplete."""
    if len(codon) != 3:
        return None
    return str(Seq(codon).translate(table=table))

def _next_aa_matches(cds_seq: str, i_after: int, aa_next: Optional[str], table: int) -> bool:
    """
    Check whether the next amino acid (aa_next) matches the next codon starting at i_after.
    If aa_next is None (we're at the last AA), consider it a match (nothing to check).
    """
    if aa_next is None:
        return True
    if i_after + 3 > len(cds_seq):
        return False
    aa = _aa_of(str(cds_seq[i_after:i_after+3]), table)
    return aa == aa_next or (aa_next == 'X')  # allow X in target as wildcard for the *next* AA too

def correct_cds(cds_seq, pep_seq, table=1):
    """
    Walk through cds_seq vs. pep_seq; remove any 1–2 bp insertions
    in cds_seq that rescue the translation frame. If that fails,
    cautiously try a 3–4 bp removal but only accept it if it also keeps the
    *following* amino acid in-frame (mitigates false positives).

    'X' in the protein is treated as a wildcard AA *only if* the following AA
    also remains in-frame (so we don't drift).

    Returns (fixed_seq:str, modifications:list of (pos, removed_seq)).
    """
    i = 0  # position in cds_seq (nt)
    j = 0  # position in pep_seq (aa)
    fixed_codons = []
    mods = []
    Lc, Lp = len(cds_seq), len(pep_seq)

    while j < Lp:
        aa_target = pep_seq[j]
        aa_next = pep_seq[j+1] if (j + 1) < Lp else None

        # 0-skip (current frame) acceptance rule:
        # - exact AA match; or
        # - if target is 'X', only accept if the *next* AA will match from i+3
        if i + 3 <= Lc:
            codon = str(cds_seq[i:i+3])
            aa0 = _aa_of(codon, table)
            if aa0 == aa_target or (
                aa_target == 'X' and _next_aa_matches(cds_seq, i + 3, aa_next, table)
            ):
                fixed_codons.append(codon)
                i += 3
                j += 1
                continue

        # Try skipping 1 or 2 bp (original behavior), with support for 'X' as wildcard
        rescued = False
        for skip in (1, 2):
            if i + skip + 3 <= Lc:
                codon2 = str(cds_seq[i + skip:i + skip + 3])
                aa2 = _aa_of(codon2, table)
                if (aa2 == aa_target) or (
                    aa_target == 'X' and _next_aa_matches(cds_seq, i + skip + 3, aa_next, table)
                ):
                    removed = str(cds_seq[i:i + skip])
                    mods.append((i, removed))
                    i += skip
                    fixed_codons.append(codon2)
                    i += 3
                    j += 1
                    rescued = True
                    break
        if rescued:
            continue

        # Fallback: try skipping 3 or 4 bp.
        # Preference: smaller skip; also *require* that the next AA matches after the skip,
        # to avoid false positives that push us out of frame.
        candidates: List[Tuple[int, str]] = []  # (skip, codon_after_skip)
        for skip in (3, 4):
            if i + skip + 3 <= Lc:
                codon_f = str(cds_seq[i + skip:i + skip + 3])
                aa_f = _aa_of(codon_f, table)
                # Accept if exact match, or 'X' with next AA remaining in-frame.
                if (aa_f == aa_target) or (
                    aa_target == 'X' and _next_aa_matches(cds_seq, i + skip + 3, aa_next, table)
                ):
                    # Require next AA to match (if exists). This is the mitigation step.
                    if _next_aa_matches(cds_seq, i + skip + 3, aa_next, table):
                        candidates.append((skip, codon_f))

        if candidates:
            # Choose the smallest skip that also keeps next AA in-frame (already enforced)
            candidates.sort(key=lambda x: x[0])
            skip, codon_f = candidates[0]
            removed = str(cds_seq[i:i + skip])
            mods.append((i, removed))
            i += skip
            fixed_codons.append(codon_f)
            i += 3
            j += 1
            continue

        # Nothing worked; raise a clear error as before
        raise ValueError(
            f"Cannot align aa #{j+1} ({aa_target}) "
            f"at cds position {i+1} (next codon != target)"
        )

    return "".join(fixed_codons), mods

def setup_logging(verbose: bool):
    level = logging.DEBUG if verbose else logging.WARNING
    logging.basicConfig(
        stream=sys.stderr,
        level=level,
        format="%(levelname)s: %(message)s"
    )

def _wrap_fasta(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def _worker(task) -> Tuple[str, Optional[str], List[Tuple[int, str]], Optional[str]]:
    """
    Worker to process a single (rid, cds_str, pep_str, table) tuple.

    Returns:
      (rid, fixed_seq_or_None, mods, error_msg_or_None)
    """
    rid, cds_str, pep_str, table = task
    try:
        fixed_seq, mods = correct_cds(cds_str, pep_str, table=table)
        return rid, fixed_seq, mods, None
    except Exception as e:
        return rid, None, [], str(e)

def main():
    p = argparse.ArgumentParser(
        description=(
            "Remove 1–2 bp insertions in CDS so that it translates exactly to the given protein; "
            "if needed, cautiously try 3–4 bp removals that also keep the following AA in-frame. "
            "Supports parallel processing."
        )
    )
    p.add_argument("-c", "--cds", required=True,
                   help="FASTA of coding sequences")
    p.add_argument("-p", "--pep", required=True,
                   help="FASTA of proteins (same IDs as CDS)")
    p.add_argument("-o", "--out", default="-",
                   help="Output corrected CDS FASTA (default: stdout)")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Show informational messages about frame corrections")
    p.add_argument("-t", "--threads", type=int, default=1,
                   help="Number of CDS/PEP pairs to process in parallel (default: 1)")
    p.add_argument("--table", type=int, default=1,
                   help="NCBI translation table (default: 1)")
    args = p.parse_args()

    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    # load into dicts
    try:
        cds_dict = SeqIO.to_dict(SeqIO.parse(args.cds, "fasta"))
        pep_dict = SeqIO.to_dict(SeqIO.parse(args.pep, "fasta"))
    except Exception as e:
        logger.error(f"Failed to read input FASTAs: {e}")
        sys.exit(1)

    # Build ordered task list using peptide order; warn on missing CDS
    tasks = []
    for rid, pep_rec in pep_dict.items():
        if rid not in cds_dict:
            logger.warning(f"No CDS record for {rid}")
            continue
        cds_rec = cds_dict[rid]
        tasks.append((rid, str(cds_rec.seq), str(pep_rec.seq), args.table))

    if not tasks:
        logger.error("No matching CDS/PEP pairs to process.")
        sys.exit(1)

    # Prepare output handle
    out_handle = sys.stdout if args.out == "-" else open(args.out, "w")

    try:
        # Process tasks: in parallel if threads > 1, otherwise serial
        if args.threads and args.threads > 1:
            # chunksize heuristic: larger chunks reduce overhead
            # Keep output order stable with imap
            with Pool(processes=args.threads) as pool:
                # Aim for ~8 chunks per worker
                chunksize = max(1, len(tasks) // (args.threads * 8) or 1)
                for rid, fixed_seq, mods, err in pool.imap(_worker, tasks, chunksize=chunksize):
                    if err is not None:
                        logger.error(f"{rid}: {err}")
                        continue
                    # report
                    if not mods:
                        logger.info(f"{rid}: no changes needed")
                    else:
                        for pos, removed in mods:
                            logger.info(
                                f"{rid}: removed {len(removed)} bp at CDS pos {pos+1}: '{removed}'"
                            )
                    # write fasta
                    out_handle.write(f">{rid}\n{_wrap_fasta(fixed_seq)}\n")
        else:
            # Serial path (no multiprocessing)
            for rid, cds_str, pep_str, table in tasks:
                try:
                    fixed_seq, mods = correct_cds(cds_str, pep_str, table=table)
                except Exception as e:
                    logger.error(f"{rid}: {e}")
                    continue
                if not mods:
                    logger.info(f"{rid}: no changes needed")
                else:
                    for pos, removed in mods:
                        logger.info(
                            f"{rid}: removed {len(removed)} bp at CDS pos {pos+1}: '{removed}'"
                        )
                out_handle.write(f">{rid}\n{_wrap_fasta(fixed_seq)}\n")
    finally:
        if args.out != "-":
            out_handle.close()

if __name__ == "__main__":
    main()
