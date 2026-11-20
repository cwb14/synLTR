"""Samples sheet: the tool's main input, one row per genome. Fails fast on any problem."""
import os
import re
from dataclasses import dataclass
from pathlib import Path

REQUIRED = ("label",)
OPTIONAL = ("fasta", "annotation", "group", "role", "library", "reads")
ROLES = ("outgroup", "reference")
LABEL_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
PAIR_SEP = "__"  # joins two labels in cache file names (<A>__<B>.<route>.npz)
RESERVED = ("all",)  # labels taken by the outputs: "all" is qc.tsv's global scope


class SampleSheetError(ValueError):
    """The samples sheet is unreadable or inconsistent."""


@dataclass(frozen=True)
class Sample:
    label: str
    fasta: Path | None  # genome FASTA, or its .fai where a command accepts one
    annotation: str | None  # LTRquest prefix or file; checked by the caller, not here
    group: str | None
    role: str | None
    library: str | None = None  # LTR library as text: LTRquest prefix, FASTA, or "FASTA,TABLE"; parsed by io.library
    reads: Path | None = None  # read alignments (BAM) of this sample


def _resolve(value: str, base: Path) -> Path:
    """`value` as an absolute path, relative ones against `base`; symlinks are kept (a FASTA's index may sit
    beside its link)."""
    p = Path(value)
    return Path(os.path.abspath(p if p.is_absolute() else base / p))


def read_samples(path: Path, *, min_rows: int = 2) -> list[Sample]:
    """Parse a tab-delimited samples sheet (UTF-8, a byte-order mark allowed); relative paths resolve against
    its directory. A sheet needs at least `min_rows` rows."""
    path = Path(path)
    if not path.is_file():
        raise SampleSheetError(f"samples sheet not found: {path}")
    base = path.resolve().parent
    header, rows = None, []
    with open(path, encoding="utf-8-sig") as fh:
        for n, line in enumerate(fh, 1):
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith("#"):
                continue
            fields = [f.strip() for f in line.split("\t")]
            if header is None:
                header = fields
                unknown = [c for c in header if c not in REQUIRED + OPTIONAL]
                if unknown:
                    raise SampleSheetError(
                        f"{path}: unknown column(s) {unknown}; the first line must be the header "
                        f"(label, fasta, annotation, optionally group, role, library and reads, tab-separated); "
                        f"allowed: "
                        f"{list(REQUIRED + OPTIONAL)}")
                if "fasta" not in header and "library" not in header and "reads" not in header:
                    raise SampleSheetError(f"{path}: header lacks required column(s) ['fasta']; a sheet needs "
                                           "fasta, library or reads")
                missing = [c for c in REQUIRED if c not in header]
                if missing:
                    raise SampleSheetError(f"{path}: header lacks required column(s) {missing}")
                if len(set(header)) != len(header):
                    raise SampleSheetError(f"{path}: duplicate column names in header")
            elif len(fields) != len(header):
                raise SampleSheetError(
                    f"{path}: line {n} has {len(fields)} fields, header has {len(header)}")
            else:
                rows.append((n, dict(zip(header, fields))))
    if header is None:
        raise SampleSheetError(f"{path}: no header line")

    def opt(r, k):
        v = r.get(k, ".")
        return None if v in (".", "") else v

    samples, seen, bad_labels, bad_seps, bad_roles = [], set(), [], [], []
    for n, r in rows:
        label = r["label"]
        if not LABEL_RE.match(label):
            bad_labels.append(f"{label!r} (line {n})")
        elif PAIR_SEP in label or label[0] == "_" or label[-1] == "_":
            bad_seps.append(f"{label!r} (line {n})")
        if label in RESERVED:
            raise SampleSheetError(f"{path}: label {label!r} (line {n}) is reserved: it names qc.tsv's global "
                                   "scope")
        if label in seen:
            raise SampleSheetError(f"{path}: duplicate label {label!r} (line {n})")
        if opt(r, "fasta") is None and opt(r, "library") is None and opt(r, "reads") is None:
            raise SampleSheetError(f"{path}: line {n} ({label}) has no FASTA; every row needs a fasta, "
                                   "library or reads")
        seen.add(label)
        role = opt(r, "role")
        if role is not None and role not in ROLES:
            bad_roles.append(f"{role!r} (line {n})")
        ann = opt(r, "annotation")
        fa, lib, rd = opt(r, "fasta"), opt(r, "library"), opt(r, "reads")
        samples.append(Sample(label, None if fa is None else _resolve(fa, base),
                              None if ann is None else str(_resolve(ann, base)), opt(r, "group"), role,
                              None if lib is None else ",".join(str(_resolve(p, base)) for p in lib.split(",")),
                              None if rd is None else _resolve(rd, base)))
    if bad_labels:
        raise SampleSheetError(f"{path}: invalid label(s) {', '.join(bad_labels)}; use only A-Z a-z 0-9 _ . -")
    if bad_seps:  # with no '__' inside and no '_' at either end, the join splits back into one pair only
        raise SampleSheetError(
            f"{path}: label(s) {', '.join(bad_seps)} contain '{PAIR_SEP}' or start or end with '_'; "
            f"'{PAIR_SEP}' is reserved: it joins two labels in cache file names")
    if bad_roles:
        raise SampleSheetError(f"{path}: unknown role(s) {', '.join(bad_roles)}; allowed: {list(ROLES)}")
    missing = [str(s.fasta) for s in samples if s.fasta is not None and not s.fasta.is_file()]
    if missing:
        raise SampleSheetError(f"{path}: FASTA file(s) not found: {', '.join(missing)}")
    missing = [str(s.reads) for s in samples if s.reads is not None and not s.reads.is_file()]
    if missing:
        raise SampleSheetError(f"{path}: reads file(s) not found: {', '.join(missing)}")
    if any(s.reads for s in samples):
        refs = [s for s in samples if s.role == "reference"]
        if len(refs) != 1 or refs[0].fasta is None or refs[0].annotation is None:
            raise SampleSheetError(f"{path}: a sheet with reads needs exactly one row with role reference, a "
                                   "fasta and an annotation")
    if len(samples) < min_rows:
        raise SampleSheetError(f"{path}: need at least {min_rows} samples, found {len(samples)}")
    if not any(s.annotation for s in samples) and not any(s.library or s.reads for s in samples):
        raise SampleSheetError(f"{path}: no sample has an annotation; at least one is required")
    return samples
