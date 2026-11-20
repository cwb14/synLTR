"""Target site duplications: the annotated TSD, and the one inferred from the flanks.

Genome coordinates: an Element is 1-based closed, `Fasta.fetch` takes 0-based half-open positions, so the left
flank of an element starting at `start` ends at 0-based position `start - 1` and the right flank begins at `end`.
"""
from __future__ import annotations

import re

from ltrtrace.base.constants import TSD_MAX_LEN, TSD_MIN_LEN
from ltrtrace.io.annotation.records import TSD_RULES, UNKNOWN, check_rule
from ltrtrace.io.elements import Element
from ltrtrace.io.seqio import Fasta

_TSD = re.compile(rf"[ACGT]{{{TSD_MIN_LEN},{TSD_MAX_LEN}}}")


def clean_tsd(x: str | None) -> str | None:
    """An annotated TSD in upper case if it is 3-8 bp of ACGT (genome orientation), else None."""
    x = (x or "").upper()
    return x if _TSD.fullmatch(x) else None


_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_ACGT = re.compile(r"[ACGT]+")


def canon_tsd(tsd: str | None) -> str:
    """A TSD in its canonical orientation: the lexicographically smaller of the sequence and its reverse
    complement, upper case, so one site read from either strand gives one string. "" for anything that is not
    a nonempty ACGT string (missing, ".", "NA", other letters)."""
    tsd = (tsd or "").upper()
    if not _ACGT.fullmatch(tsd):
        return ""
    return min(tsd, tsd[::-1].translate(_COMPLEMENT))


def annotated_tsd(x: str | None, rule: str = "keep", null: set[str] = UNKNOWN) -> str | None:
    """An annotated TSD under a tool's rule, in upper case.

    "keep" (LTRpav): any value that is not missing or in `null` is kept. "clean" (ltrkin): ACGT, 3-8 bp
    (`clean_tsd`). "acgtn" (phasing): surrounding blanks dropped, any length, only the letters ACGTN, and "."
    and NA in any case are none.
    """
    check_rule("tsd_rule", rule, TSD_RULES)
    if rule == "clean":
        return clean_tsd(x)
    if rule == "acgtn":
        x = (x or "").strip().upper()
        return x if x and x not in (".", "NA") and set(x) <= set("ACGTN") else None
    return None if x is None or x in null else x.upper()


def infer_tsd(fasta: Fasta, e: Element, order: tuple[int, ...]) -> str | None:
    """The direct repeat flanking the element, trying the lengths of `order` in turn (the first that repeats
    wins), or None. A flank with an N or running off a sequence end does not count."""
    return infer_tsd_at(fasta, e.seqid, e.start, e.end, order)


def infer_tsd_at(fasta: Fasta, seqid: str, start: int, end: int, order: tuple[int, ...]) -> str | None:
    """`infer_tsd` for the element `seqid:start-end` (1-based closed)."""
    zstart = start - 1
    for n in order:
        if zstart - n < 0:
            continue
        left = fasta.fetch(seqid, zstart - n, zstart)
        right = fasta.fetch(seqid, end, end + n)
        if len(left) == n and len(right) == n and left == right and "N" not in left:
            return left
    return None
