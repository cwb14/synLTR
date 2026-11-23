"""manifest.json: the record of a finished run in its output folder (what was run, on what, with which tools), for
reproducibility; `synltr plot` reads which run it was, its samples sheet and its route from it, and `synltr
source` the reference of a read genotyping run.

A command that finishes a run returns its record (`Finished`); the command line writes it with `record_run`. A
command that will write one removes the folder's earlier manifest.json once its own checks have passed and before it
writes anything (`remove_manifest`), so a run that fails leaves no record of other results. The manifest holds the
command line as given; the parsed options (paths absolute) with the run's own values over them, and beside --paf the
checksum of each PAF file the anchors keyed on (`paf_checksums`); the evidence route (assembly, library or reads;
none for subgenomes and source); each sample of the sheet with the checksums of its inputs that the run's stages keyed
on (`sample_checksums`, from base.checkpoint.keyed_checksums: nothing is hashed for the manifest, and nothing is read
back from the checksum cache) and its annotation loader's skip counts; the code fingerprint of the command's module
and its imports; Python, its packages and every external tool the engine can run, with their versions; the start and
end times (UTC) and the host."""
from __future__ import annotations

import json
import platform
import socket
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path

from synltr import __version__
from synltr.base.checkpoint import atomic_write, code_fingerprint, tool_version
from synltr.base.constants import MANIFEST, TOOL_VERSION_ARGS
from synltr.base.log import UsageError
from synltr.io.samples import Sample

_PYTHON_DEPS = ("numpy", "numba", "pysam", "mappy", "scipy", "matplotlib")
_OPTIONAL_DEPS = ("msprime",)  # recorded only when installed (the simulator's extra)


@dataclass(frozen=True)
class Finished:
    """A finished run in `outdir`, as its command returns it for manifest.json: the samples sheet's rows it read,
    its evidence route (None for subgenomes and source), the checksums of their inputs its stages keyed on (label ->
    {input: sha256}; `sample_checksums`), each annotated genome's loader skip counts (label -> {reason: n}), and the
    options it resolved to other values than the parsed ones with what it records beside them (`params`, e.g. a
    default --ploidy, or the checksums of the --paf files)."""

    outdir: Path
    samples: tuple[Sample, ...]
    route: str | None = None
    checksums: Mapping[str, Mapping[str, str]] = field(default_factory=dict)
    skipped: Mapping[str, Mapping[str, int]] = field(default_factory=dict)
    params: Mapping[str, object] = field(default_factory=dict)


def _py_version(name: str) -> str | None:
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return None


def _now() -> datetime:
    return datetime.now(UTC).replace(microsecond=0)


def write_manifest(path: Path, argv: list[str], params: dict, samples: list[Sample],
                   checksums: Mapping[str, Mapping[str, str]], started: datetime, finished: datetime,
                   skipped: Mapping[str, Mapping[str, int]] | None = None, units: Iterable[str] = (),
                   route: str | None = None) -> None:
    """Write the manifest at `path` atomically. checksums: genome label -> {input name: checksum}; skipped: genome
    label -> its annotation loader's skip counts by reason (malformed records included; {} for a genome without
    one); route: the evidence route. Dependencies: Python, its packages (msprime only when installed) and every
    external tool the engine can run (checkpoint.tool_version; None when absent); code: the code fingerprint of
    each unit in `units` (checkpoint.code_fingerprint), which the stage keys hold."""
    optional = {n: v for n in _OPTIONAL_DEPS if (v := _py_version(n)) is not None}
    deps = {"python": platform.python_version(), **{n: _py_version(n) for n in _PYTHON_DEPS}, **optional,
            **{exe: tool_version(exe) for exe in TOOL_VERSION_ARGS}}
    skipped = skipped or {}
    doc = {
        "tool": "synltr", "version": __version__, "code": {u: code_fingerprint(u) for u in units},
        "command": list(argv),
        "params": params,
        "route": route,
        "samples": [{"label": s.label, "fasta": None if s.fasta is None else str(s.fasta), "annotation": s.annotation,
                     "library": s.library, "reads": None if s.reads is None else str(s.reads), "group": s.group,
                     "role": s.role, "checksums": dict(checksums.get(s.label, {})),
                     "skipped": {k: int(v) for k, v in sorted(skipped.get(s.label, {}).items()) if v}}
                    for s in samples],
        "dependencies": deps, "started": started.isoformat(), "finished": finished.isoformat(),
        "host": socket.gethostname(),
    }
    atomic_write(Path(path), json.dumps(doc, indent=2, sort_keys=True, default=str) + "\n")


def record_run(run: Finished, argv: list[str], options: Mapping[str, object], started: datetime,
               unit: str | None) -> Path:
    """Write manifest.json of the finished run `run` in run.outdir (`write_manifest`) and return its path: `argv`
    as given, `options` (the parsed options as JSON values) with run.params over them, run's samples, checksums,
    skip counts and route, the code fingerprint of `unit` (the command's module; None: none), `started` and now."""
    path = Path(run.outdir) / MANIFEST
    path.parent.mkdir(parents=True, exist_ok=True)
    write_manifest(path, argv, {**options, **run.params}, list(run.samples), run.checksums, started, _now(),
                   skipped=run.skipped, units=() if unit is None else (unit,), route=run.route)
    return path


def read_manifest(folder: Path, remedy: str, *, absent: str | None = None) -> dict:
    """The manifest of the finished run in `folder`, its "params" a dict ({} when it records none). Without one, a
    UsageError of one line: `absent` when given, else one naming the folder, then `remedy` (what to run); with one that
    cannot be read or is not a manifest, one naming the file, then `remedy`."""
    path = Path(folder) / MANIFEST
    try:
        doc = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        raise UsageError(absent or f"{folder} holds no finished run: no {MANIFEST}, which a run writes when it "
                                   f"finishes; {remedy}") from None
    except (OSError, ValueError) as exc:
        raise UsageError(f"{path} is unreadable ({exc}); {remedy}") from None
    if not isinstance(doc, dict) or not isinstance(doc.setdefault("params", {}), dict):
        raise UsageError(f"{path} is not an synltr manifest; {remedy}")
    return doc


def remove_manifest(outdir: Path) -> None:
    """Remove the manifest.json of an earlier run in `outdir`. A command that will write one calls this once its own
    checks have passed and before it writes anything: a usage error leaves a finished run's record as it was, and a
    run that fails after this leaves no record that describes other results."""
    (Path(outdir) / MANIFEST).unlink(missing_ok=True)


def sample_checksums(samples: Iterable[Sample], keyed: Mapping[str, str],
                     kinds: Iterable[str]) -> dict[str, dict[str, str]]:
    """The checksums of the samples' inputs of `kinds` ("fasta", "annotation", "library": those the run's route reads)
    as its stages keyed on them: `keyed` (base.checkpoint.keyed_checksums, absolute path -> sha256). Nothing is read or
    hashed here, and an input no stage keyed on is left out. Names: "fasta"; "annotation" for a GFF3 or BED file, or
    "annotation:<file name>" for each table of an LTRquest prefix; "library:<file name>" for each FASTA and table of
    a library. Returns label -> {name: sha256}."""
    from synltr.io.annotation import annotation_files
    from synltr.io.library import library_files

    kinds = set(kinds)
    out = {}
    for s in samples:
        files: list[tuple[str, Path]] = []
        if "fasta" in kinds and s.fasta is not None:
            files.append(("fasta", Path(s.fasta)))
        if "annotation" in kinds and s.annotation is not None:
            found = annotation_files(s.annotation)
            one = found == [Path(s.annotation)]
            files += [("annotation" if one else f"annotation:{f.name}", f) for f in found]
        if "library" in kinds and s.library is not None:
            fastas, tables = library_files(s.library)
            files += [(f"library:{f.name}", Path(f)) for f in (*fastas, *tables)]
        out[s.label] = {name: keyed[k] for name, f in files if (k := str(Path(f).resolve())) in keyed}
    return out


def paf_checksums(folder: Path, keyed: Mapping[str, str]) -> dict[str, str]:
    """The checksums of the files of a --paf folder that the run's stages keyed on (`keyed`: absolute path ->
    sha256), by file name: the whole-genome alignments the anchors were built from."""
    return {f.name: keyed[k] for f in sorted(Path(folder).iterdir()) if (k := str(f.resolve())) in keyed}


def manifest_samples(doc: Mapping) -> list[Sample]:
    """The samples sheet's rows a manifest (`read_manifest`) records."""
    def path(value):
        return None if value is None else Path(value)

    return [Sample(s["label"], path(s.get("fasta")), s.get("annotation"), s.get("group"), s.get("role"),
                   s.get("library"), path(s.get("reads"))) for s in doc.get("samples", [])]
