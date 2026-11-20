"""Helpers every annotation loader shares: skipped records (counts by reason, and a warning naming the file and
line of each of the first malformed records) and the reading of classification strings."""
from __future__ import annotations

import logging
import math
import re
from collections import Counter
from pathlib import Path

from ltrtrace.base.constants import RECORD_WARNINGS

log = logging.getLogger("ltrtrace.io.annotation")

UNKNOWN = {"", ".", "unknown", "NA"}  # values a GFF3 uses for "not known"
CLASS_RULES = ("split", "last")
TSD_RULES = ("keep", "clean", "acgtn")  # see tsd.annotated_tsd


def check_rule(name: str, value: str, allowed: tuple[str, ...]) -> None:
    """Raise ValueError unless `value` is one of `allowed`: an unknown rule never acts as the default."""
    if value not in allowed:
        raise ValueError(f"{name} must be one of {allowed}, not {value!r}")


def known(v: str | None) -> str | None:
    """`v`, or None when it is missing or one of the UNKNOWN spellings."""
    return None if v is None or v in UNKNOWN else v


def number(v: str | None) -> float | None:
    """A finite number from text; None when absent, unreadable or not finite."""
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def own_length(v: str | None) -> int | None:
    """An LTR's own length from text: a positive whole number, else None (not given)."""
    x = number(v)
    return None if x is None or x < 1 or x != int(x) else int(x)


def class_fields(cls: str | None, rule: str, superfamily: str | None = None, clade: str | None = None,
                 unknown: set[str] = UNKNOWN, dante: bool = False) -> tuple[str | None, str | None]:
    """(superfamily, clade) from a classification string and the explicit `superfamily` and `clade` attributes
    (raw text, as the annotation gives them).

    rule "split" (LTRpav): the explicit attributes count unless missing or in `unknown`; a classification
    `LTR/Copia/Ale` fills only what they leave unknown, superfamily from its second field and clade from its
    third. Rule "last" (ltrkin): the superfamily is the attribute if known; the clade is the attribute unless it
    is missing, empty or ".", and then the last field of the classification (split on `/` or `|`), kept as
    written ("unknown" too). `dante` reads a DANTE_LTR classification `Class_I|LTR|Ty3/gypsy|chromovirus|Tekay`:
    the superfamily is the field after LTR and the clade is the last field.
    """
    check_rule("class_rule", rule, CLASS_RULES)
    sf = None if superfamily is None or superfamily in unknown else superfamily
    if rule == "last":
        if clade in (None, "", ".") and cls:
            clade = re.split(r"[/|]", cls)[-1]
        sf = known(superfamily)
        return sf, None if clade in (None, "", ".") else clade
    clade = None if clade is None or clade in unknown else clade
    if cls and (sf is None or clade is None):
        if dante:
            parts = [p for p in cls.split("|") if not p.startswith("Class_")]
            fields = (parts[1] if len(parts) > 1 else None, parts[-1] if len(parts) > 2 else None)
        else:
            parts = cls.split("/")
            fields = (parts[1] if len(parts) > 1 else None, parts[2] if len(parts) > 2 else None)
        sf = sf if sf is not None else (None if fields[0] in unknown else fields[0])
        clade = clade if clade is not None else (None if fields[1] in unknown else fields[1])
    return sf, clade


class RecordSkips:
    """A loader's skip counts by reason (`counts`). `malformed` skips a record that cannot be read, with a
    warning naming its file and line for the first RECORD_WARNINGS of them (later ones are only counted);
    `filtered` counts a record the loader leaves out by design (no LTR pair, LTRquest status not pass)."""

    def __init__(self) -> None:
        self.counts: Counter = Counter()
        self._warned = 0

    def malformed(self, path: Path | str, line: int, reason: str, why: str) -> None:
        self.counts[reason] += 1
        if self._warned < RECORD_WARNINGS:
            self._warned += 1
            more = "; further malformed records are only counted" if self._warned == RECORD_WARNINGS else ""
            log.warning("%s:%d: record skipped, %s: %s%s", path, line, reason, why, more)

    def filtered(self, reason: str) -> None:
        self.counts[reason] += 1
