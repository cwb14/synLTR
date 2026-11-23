"""Read and write LTR-RT GFF3 (LTRquest, EDTA, LTR_retriever, LTRharvest and DANTE_LTR conventions), and the one
GFF3 escaping every writer of the package uses (escape_seqid, escape_value; unescape, the reader's inverse).

GFF3 is 1-based closed, like the elements: no conversion happens here."""
from __future__ import annotations

import gzip
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote

from synltr.base.constants import GFF3_SOURCE
from synltr.io.annotation.records import (
    CLASS_RULES,
    TSD_RULES,
    RecordSkips,
    check_rule,
    class_fields,
    known,
    number,
    own_length,
)
from synltr.io.annotation.tsd import annotated_tsd
from synltr.io.elements import NO_LTR, Element, make_eid
from synltr.io.files import open_input
from synltr.seq.divergence import jc69_from_identity

_VALUE_ESC = re.compile(r"[;=&,%\x00-\x1f\x7f]")  # GFF3 column 9: reserved characters, tabs, line breaks, controls
_SEQID_ESC = re.compile(r"[^a-zA-Z0-9.:^*$@!+_?|-]")  # GFF3 column 1: characters outside this set are %XX
_COORD_LIKE = re.compile(r"^[^:]+:\d+[-.]+\d+")  # a Name that is only the element's coordinates
K_KEYS = ("K2P_d", "k2p", "K2P", "ltr_identity")  # attributes that can give k; ltr_identity gives a JC69 distance
_LTR = "long_terminal_repeat"
_HOST = "LTR_retrotransposon"  # a type ending in it is an LTR-RT (EDTA: Copia_LTR_retrotransposon, ...)
_DANTE = "transposable_element"  # DANTE_LTR's LTR-RT: read when it carries a Final_Classification
_DANTE_CLASS = "Final_Classification"
LTR_RULES = ("two", "outermost", "optional")
NAME_RULES = ("name", "id")
HOST_RULES = ("type", "k2p")
_K2P_HOSTS = ("LTR_retrotransposon", "repeat_region")  # host_rule "k2p": these types exactly, carrying a K2P key
_K2P_KEYS = ("K2P_d", "k2p", "K2P")


def _attrs(col9: str) -> tuple[dict[str, str], list[str]]:
    """Column 9 as (attribute -> value, the Parent ids), percent-decoded; Parent is split on its raw commas,
    then each id decoded."""
    out, parents = {}, []
    for part in col9.strip().split(";"):
        if "=" in part:
            key, val = part.split("=", 1)
            key = key.strip()
            out[key] = unescape(val)
            if key == "Parent":
                parents = [unescape(x) for x in val.split(",") if x]
    return out, parents


def _num(a: dict[str, str], *keys: str) -> float | None:
    """The finite number under the first of `keys` present in `a`; None if absent, unreadable or not finite."""
    return next((number(a[k]) for k in keys if k in a), None)


def _divergence(a: dict[str, str], k_keys: tuple[str, ...]) -> tuple[float | None, str]:
    """(k, kind) from the first of `k_keys` present: K2P from K2P_d, k2p or K2P; the JC69 distance from EDTA's
    ltr_identity. k is None if the value is unreadable or not finite (an identity of nan gives nan, as in
    LTRpav)."""
    key = next((k for k in k_keys if k in a), None)
    if key is None:
        return None, "K2P"
    if key != "ltr_identity":
        return number(a[key]), "K2P"
    try:
        identity = float(a[key])
    except ValueError:
        return None, "JC69"
    return jc69_from_identity(identity), "JC69"


def _ltrs_of(host: tuple, ltrs_by_parent: dict, ltr_rule: str, key: str | None) -> list[tuple[int, int]]:
    """The LTRs linked to a host: its children (under ID `key`), else (EDTA, LTR_retriever) the LTRs under its
    repeat_region that lie inside it."""
    seqid, start, end, _, _, parents, _ = host
    ltrs = list(ltrs_by_parent.get((seqid, key), [])) if key else []
    short = len(ltrs) < 2 if ltr_rule == "outermost" else not ltrs
    if short:
        for p in (parents if ltr_rule == "outermost" else parents[:1]):
            ltrs += [x for x in ltrs_by_parent.get((seqid, p), []) if x[0] >= start and x[1] <= end]
    return ltrs


def load_gff3(path: Path, genome: str, ltr_rule: str = "two", class_rule: str = "split", tsd_rule: str = "keep",
              name_family: bool = False, k_keys: tuple[str, ...] = ("K2P_d", "ltr_identity"),
              name_rule: str = "name", host_rule: str = "type") -> tuple[list[Element], Counter]:
    """Return (elements, skip counts).

    Only LTR-RT lines and long_terminal_repeat lines are parsed: a line of any other feature type, however
    malformed, is never read. Under `host_rule` "type" (LTRpav, ltrkin) a type ending in LTR_retrotransposon is an
    LTR-RT, and so is DANTE_LTR's transposable_element when it carries a Final_Classification; under "k2p" (phasing)
    an LTR_retrotransposon or repeat_region line (those types exactly) that carries K2P_d, k2p or K2P is one. One of
    these with fewer than 9 columns (malformed_columns) or a start and end that are not integers with 1 <= start <=
    end (malformed_coordinates) is skipped with a warning naming the file and line (RecordSkips). Seqids and
    attribute values are percent-decoded.

    `ltr_rule` "two" (LTRpav): an LTR-RT without exactly two LTRs is no_ltr_pair. "outermost" is ltrkin's
    reader: types are matched in any case, every transposable_element is an LTR-RT, the outermost two of two or
    more LTRs bound the element (its span extended to their ends), and a repeated element ID, or coordinate ID
    for an element without one, is an error. "optional" (phasing): exactly two LTRs give the extents, and an
    LTR-RT with any other number is kept with both extents NO_LTR. `class_rule` is read by `class_fields` (under
    "last" a Final_Classification is read on any element), `tsd_rule` by `annotated_tsd` (rules other than "keep"
    also read a `TSD` attribute). With `name_family` (ltrkin) the family is the family or Family attribute as written,
    else the Name unless that is only coordinates. `k_keys` lists the attributes that give k, in order of
    precedence: LTRpav's are K2P_d then ltr_identity (a JC69 distance); phasing's are K2P_d, k2p, K2P. A record
    with a Classification and no `family` attribute is EDTA's: its Name is the family and its ID the element
    name. `name_rule` "name" (LTRpav) names an element by its Name, else its ID; "id" (ltrkin) by its ID, else its
    line's coordinates seqid:start-end (as written, before the outermost rule extends the span).
    """
    check_rule("ltr_rule", ltr_rule, LTR_RULES)
    check_rule("name_rule", name_rule, NAME_RULES)
    check_rule("class_rule", class_rule, CLASS_RULES)
    check_rule("tsd_rule", tsd_rule, TSD_RULES)
    check_rule("host_rule", host_rule, HOST_RULES)
    if not set(k_keys) <= set(K_KEYS):
        raise ValueError(f"k_keys must be among {K_KEYS}, not {k_keys}")
    path = Path(path)
    outer = ltr_rule == "outermost"
    skips = RecordSkips()
    hosts: list[tuple[str, int, int, str, dict, list[str], bool]] = []
    ltrs_by_parent: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    with open_input(path) as fh:
        for n, line in enumerate(fh, 1):
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\r\n").split("\t")
            ftype = f[2] if len(f) > 2 else ""
            kind = ftype.lower() if outer else ftype
            if host_rule == "k2p":
                is_dante, is_host = False, kind in _K2P_HOSTS
            else:
                is_dante = kind == _DANTE and (outer or _DANTE_CLASS in line)
                is_host = kind.endswith(_HOST.lower() if outer else _HOST) or is_dante
            if kind != _LTR and not is_host:
                continue  # a feature type the loader does not use
            if len(f) < 9:
                skips.malformed(path, n, "malformed_columns", f"{ftype} line with {len(f)} tab-separated columns, "
                                "not 9")
                continue
            try:
                start, end = int(f[3]), int(f[4])
            except ValueError:
                start = end = 0
            if not 1 <= start <= end:
                skips.malformed(path, n, "malformed_coordinates", f"{ftype} start {f[3]!r} and end {f[4]!r} are "
                                "not integers with 1 <= start <= end")
                continue
            seqid = unescape(f[0])
            a, parents = _attrs(f[8])
            if kind == _LTR:
                for p in parents:
                    ltrs_by_parent[(seqid, p)].append((start, end))
            elif host_rule == "type" or any(k in a for k in _K2P_KEYS):
                hosts.append((seqid, start, end, f[6], a, parents, is_dante))
    elements: list[Element] = []
    seen: set[str] = set()
    for host in hosts:
        seqid, start, end, strand, a, parents, dante = host
        eid_attr = a.get("ID")
        key = eid_attr or (f"{seqid}:{start}-{end}" if outer else None)
        if outer:
            if key in seen:  # IDs link the LTRs: a repeated one merges LTRs of different elements
                raise ValueError(f"{path}: element ID {key} occurs more than once (concatenated GFF3 files?); "
                                 "element IDs must be unique")
            seen.add(key)
        ltrs = _ltrs_of(host, ltrs_by_parent, ltr_rule, key)
        if outer:
            ltrs = sorted(set(ltrs))
            if len(ltrs) < 2:
                skips.filtered("no_ltr_pair")
                continue
            left, right = ltrs[0], ltrs[-1]
            start, end = min(start, left[0]), max(end, right[1])
        elif len(ltrs) == 2:
            left, right = sorted(ltrs)
        elif ltr_rule == "optional":
            left = right = NO_LTR
        else:
            skips.filtered("no_ltr_pair")
            continue
        tsd = annotated_tsd(a.get("tsd") or (a.get("TSD") if tsd_rule == "clean" else None), tsd_rule)
        cls = a.get("classification") or a.get("Classification")
        if not cls and (dante or class_rule == "last"):
            cls = a.get(_DANTE_CLASS)
        sf, clade = class_fields(cls, class_rule, a.get("superfamily"), a.get("clade"), dante=dante)
        name, family = a.get("Name"), known(a.get("family"))
        if "family" not in a and "Classification" in a:
            # EDTA: Name is the family (TE_00000123) and ID the unique element
            name, family = eid_attr, known(a.get("Name"))
        if name_rule == "id":  # ltrkin: the element's identity is its ID (or line coordinates), which must be unique
            name = eid_attr or f"{seqid}:{host[1]}-{host[2]}"
        if name_family:  # ltrkin: the family as written, else the Name unless that is only coordinates
            family = a.get("family") or a.get("Family")
            if not family:
                nm = a.get("Name", "")
                family = nm if nm and not _COORD_LIKE.match(nm) else None
            family = None if family == "." else family
        k, k_kind = _divergence(a, k_keys)
        elements.append(Element(
            genome=genome, eid=make_eid(genome, seqid, start, end),
            name=name or eid_attr or make_eid(genome, seqid, start, end),
            seqid=seqid, start=start, end=end,
            strand=strand if strand in ("+", "-") else ".",
            ltr_left=left, ltr_right=right,
            tsd=tsd, tsd_source="annot" if tsd else None,
            k=k, family=family, superfamily=sf, clade=clade,
            confidence=known(a.get("confidence")),
            ltr_left_len=own_length(a.get("ltr5_len")), ltr_right_len=own_length(a.get("ltr3_len")),
            k_kind=k_kind, k_se=_num(a, "k2p_se"), n_sites=_num(a, "n_sites", "aln_len"),
            n_ts=_num(a, "ti"), n_tv=_num(a, "tv"),
        ))
    return elements, skips.counts


def escape_seqid(seqid: str) -> str:
    """A seqid as GFF3 writes it in column 1 and in ##sequence-region: every character outside
    [a-zA-Z0-9.:^*$@!+_?|-] as %XX of its UTF-8 bytes (PanSN HG002#1#chr1 -> HG002%231%23chr1)."""
    return _SEQID_ESC.sub(lambda m: "".join(f"%{b:02X}" for b in m.group().encode()), seqid)


def escape_value(value: str) -> str:
    """An attribute value as GFF3 writes it in column 9: ; = & , %, tabs, line breaks and the other control
    characters as %XX."""
    return _VALUE_ESC.sub(lambda m: f"%{ord(m.group()):02X}", value)


def unescape(text: str) -> str:
    """A GFF3 seqid or attribute value with every %XX decoded: the inverse of escape_seqid and escape_value."""
    return unquote(text)


def _k_attribute(e: Element) -> tuple[str, str | None]:
    """The attribute that carries k: K2P_d, or for a JC69 distance the ltr_identity it came from (so that reading
    the file back keeps the label)."""
    if e.k is None:
        return "K2P_d", None
    if e.k_kind == "JC69":
        return "ltr_identity", repr(1.0 - 0.75 * (1.0 - math.exp(-4.0 * e.k / 3.0)))
    return "K2P_d", repr(e.k)


def write_gff3(elements: list[Element], path: Path) -> None:
    """Write LTRquest-convention GFF3 (gzip if the path ends in .gz), its source column GFF3_SOURCE (synltr and its
    version); k as K2P_d, or as ltr_identity when it is a JC69 distance."""
    lines = ["##gff-version 3\n"]
    for e in elements:
        eid = escape_value(e.eid)
        opt = [("Name", e.name), ("tsd", e.tsd), _k_attribute(e),
               ("family", e.family), ("superfamily", e.superfamily), ("clade", e.clade),
               ("confidence", e.confidence)]
        attrs = f"ID={eid}" + "".join(f";{k}={escape_value(v)}" for k, v in opt if v is not None)
        head = f"{escape_seqid(e.seqid)}\t{GFF3_SOURCE}\t"
        tail = f"\t.\t{e.strand}\t.\t"
        lines.append(f"{head}LTR_retrotransposon\t{e.start}\t{e.end}{tail}{attrs}\n")
        for tag, (s, t) in (("lLTR", e.ltr_left), ("rLTR", e.ltr_right)):
            lines.append(f"{head}long_terminal_repeat\t{s}\t{t}{tail}ID={eid}.{tag};Parent={eid}\n")
    p = Path(path)
    if str(p).endswith(".gz"):
        with open(p, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as gz:
            gz.write("".join(lines).encode())
    else:
        p.write_text("".join(lines))
