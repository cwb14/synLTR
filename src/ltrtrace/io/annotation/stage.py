"""The elements stage: one genome's annotation loaded and normalized once, then read back from its table.

<work>/elements/<label>.tsv (io.elements: the 17 element columns, 1-based closed) and <label>.skipped.json (the
loader's skip counts), keyed by the FASTA and annotation checksums (`input_checksums`), the loader's library
versions and the code of this module's import closure (the stage, the loaders and the element model). The
elements are always read back from the table, so a fresh run and a rerun hand the later stages the same elements;
their k_kind (K2P or JC69), which elements.tsv does not hold, is kept in the stage's done marker and
restored on read-back.
"""
from __future__ import annotations

import json
import logging
from collections import Counter
from dataclasses import replace
from pathlib import Path

from ltrtrace.base.checkpoint import (
    atomic_path,
    atomic_write,
    clear_done,
    code_fingerprint,
    done_meta,
    file_checksum,
    is_done,
    mark_done,
    marker,
    stage_key,
)
from ltrtrace.base.log import LtrtraceError
from ltrtrace.io.annotation import RecordSkips, annotation_files, load_annotation, loader_versions
from ltrtrace.io.elements import Element, read_elements, write_elements
from ltrtrace.io.samples import Sample
from ltrtrace.io.seqio import Fasta

log = logging.getLogger(__name__)


def input_checksums(s: Sample, sums: Path) -> dict[str, str]:
    """sha256 of the genome's inputs: "fasta"; "annotation" for a GFF3/BED file, or "annotation:<name>" for each
    table of an LTRquest prefix (none if no file matches: load_annotation then says why)."""
    out = {"fasta": file_checksum(s.fasta, sums)}
    if s.annotation is not None:
        files = annotation_files(s.annotation)
        one = files == [Path(s.annotation)]
        for f in files:
            out["annotation" if one else f"annotation:{f.name}"] = file_checksum(f, sums)
    return out


def elements_stage(label: str, spec: str, fasta: Fasta, work: Path,
                   inputs: dict[str, str]) -> tuple[list[Element], Counter]:
    """The cached elements stage of genome `label`: its normalized elements and its loader's skip counts,
    through <work>/elements/<label>.tsv and <label>.skipped.json, keyed by `inputs` (the FASTA and annotation
    checksums), the loader's library versions and this module's code. The elements are read back from the file,
    each with its k_kind from the done marker's meta: one label for the genome, or one per element in file order
    when the annotation mixes them."""
    d = Path(work) / "elements"
    tsv, counts = d / f"{label}.tsv", d / f"{label}.skipped.json"
    done = marker(tsv)
    key = stage_key("elements", inputs, {"genome": label, **loader_versions(spec)},
                    code_fingerprint("io.annotation.stage"))
    meta = done_meta(done) if is_done(done, key) else None
    if meta is not None and "k_kind" in meta:
        n = sum(json.loads(counts.read_text()).values())
        log.info("%s: elements done, skip (%s; %d annotation records skipped)", label, tsv, n)
    else:
        clear_done(done)
        skipped = RecordSkips()
        found = list(load_annotation(spec, label, fasta, skipped=skipped))
        d.mkdir(parents=True, exist_ok=True)
        with atomic_path(tsv) as tmp:
            write_elements(found, tmp)
        atomic_write(counts, json.dumps(dict(sorted(skipped.counts.items()))) + "\n")
        kinds = [e.k_kind for e in found]
        meta = {"k_kind": kinds[0] if len(set(kinds)) == 1 else kinds}
        mark_done(done, key, [tsv, counts], meta)
        log.info("%s: %d elements (%d top-level), %d annotation records skipped -> %s", label, len(found),
                 sum(e.level == 0 for e in found), sum(skipped.counts.values()), tsv)
    kind = meta["k_kind"]
    elements = list(read_elements(tsv))
    kinds = kind if isinstance(kind, list) else [kind] * len(elements)
    if len(kinds) != len(elements):
        raise LtrtraceError(f"{done}: {len(kinds)} k_kind labels for {len(elements)} elements; delete {done} and "
                            "rerun")
    return ([replace(e, k_kind=k) for e, k in zip(elements, kinds)],
            Counter(json.loads(counts.read_text())))
