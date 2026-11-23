"""Stage checkpoints: cached input checksums, stage keys and done markers.

A stage runs as: key = stage_key(...); if is_done(marker, key): load the outputs; else clear_done(marker),
write the outputs, mark_done(marker, key, outputs). The key hashes the stage name, its inputs' checksums, its
parameters, the synltr version and the fingerprint of the code the stage runs (`code_fingerprint`), so a
change in any of them reruns the stage: an outdir finished by older code is recomputed, never reused. The
versions of the external tools a stage runs (`tool_version`) go into its parameters. The marker (JSON) is
written last and atomically; it records the key, any metadata the stage adds, and each output's size,
relative to the marker's directory.

Code fingerprints: a stage is keyed on the import closure of the module (or package) holding its code
(`code_files`): that code and every synltr module it imports, directly or not, so an edit to any code the stage
runs reruns it and an edit elsewhere never does. Reads are keyed on `bam_key`, which avoids checksumming a
multi-Gb BAM.

A run's record of its inputs (calls.manifest) holds the checksums its stages keyed on, collected as file_checksum
returns them (`keyed_checksums`), never read back from the checksum cache, which may have changed since (the file
edited, another process sharing the cache).
"""
from __future__ import annotations

import ast
import contextlib
import functools
import hashlib
import json
import logging
import os
import shutil
import subprocess
import threading
from collections.abc import Iterator
from datetime import UTC, datetime
from pathlib import Path

from synltr import __version__
from synltr.base import constants

log = logging.getLogger("synltr.checkpoint")

_ROOT = Path(__file__).resolve().parent.parent  # the synltr package: code_fingerprint resolves units under it


def code_fingerprint(unit: str) -> str:
    """sha256 (hex) of the code a stage of `unit` runs: the files of `code_files(unit)`, each as its path relative
    to the package, its length and its bytes, in path order. Cached per process."""
    return _fingerprint(_ROOT, unit)


def code_files(unit: str) -> list[str]:
    """The import closure of `unit` inside synltr, as sorted paths relative to the package: the *.py files of the
    unit (a dotted name under synltr, such as "routes.assembly.pair": a package folder, all of it, or a single
    module) and of every synltr module they import, transitively. Imports are read from the source (a static
    scan of imports anywhere in a module, functions included; the body of `if TYPE_CHECKING:` is skipped);
    `from pkg import name` counts the module pkg.name when there is one, else pkg's __init__, which otherwise
    counts only where something is imported from it. An unknown unit is a ValueError."""
    return [p.relative_to(_ROOT).as_posix() for p in _closure(_ROOT, unit)]


@functools.cache
def _fingerprint(root: Path, unit: str) -> str:
    h = hashlib.sha256()
    for path in _closure(root, unit):
        data = path.read_bytes()
        h.update(path.relative_to(root).as_posix().encode() + b"\0" + len(data).to_bytes(8, "little"))
        h.update(data)
    return h.hexdigest()


def _closure(root: Path, unit: str) -> list[Path]:
    """The files of `unit` in the package directory `root` and of the synltr modules they import, transitively,
    in path order."""
    seen: set[Path] = set()
    todo = _source_files(root, unit)
    while todo:
        path = todo.pop()
        if path not in seen:
            seen.add(path)
            todo.extend(_imported(root, path))
    return sorted(seen, key=lambda p: p.relative_to(root).as_posix())


def _source_files(root: Path, unit: str) -> list[Path]:
    """The *.py files of `unit` (dotted, under synltr) in the package directory `root`."""
    path = root.joinpath(*unit.split("."))
    if path.is_dir():
        return sorted(path.rglob("*.py"))
    if path.with_suffix(".py").is_file():
        return [path.with_suffix(".py")]
    raise ValueError(f"no synltr unit {unit!r} under {root}")


def _imported(root: Path, path: Path) -> list[Path]:
    """The files of the synltr modules that the module at `path` (under `root`) imports."""
    package = list(path.relative_to(root).parts[:-1])  # where a relative import starts: the module's package
    found = []
    for node in _imports(ast.parse(path.read_bytes())):
        if isinstance(node, ast.Import):
            names = [a.name for a in node.names]
        else:
            base = node.module or ""
            if node.level:
                base = ".".join(["synltr", *package[:len(package) - node.level + 1], *([base] if base else [])])
            names = [f"{base}.{a.name}" if _module_file(root, f"{base}.{a.name}") else base for a in node.names]
        found += [f for f in (_module_file(root, name) for name in names) if f is not None]
    return found


def _module_file(root: Path, name: str) -> Path | None:
    """The file of the module `name` (dotted, absolute) when it is an synltr module under `root`, else None."""
    head, *parts = name.split(".")
    if head != "synltr":
        return None
    path = root.joinpath(*parts)
    if parts and path.with_suffix(".py").is_file():
        return path.with_suffix(".py")
    return path / "__init__.py" if (path / "__init__.py").is_file() else None


def _imports(node: ast.AST):
    """Every import statement under `node`, skipping the body of `if TYPE_CHECKING:` (its else branch counts)."""
    for child in ast.iter_child_nodes(node):
        if isinstance(child, ast.If) and _is_type_checking(child.test):
            for stmt in child.orelse:
                yield from [stmt] if isinstance(stmt, (ast.Import, ast.ImportFrom)) else _imports(stmt)
        elif isinstance(child, (ast.Import, ast.ImportFrom)):
            yield child
        else:
            yield from _imports(child)


def _is_type_checking(test: ast.expr) -> bool:
    return (isinstance(test, ast.Name) and test.id == "TYPE_CHECKING") or (
        isinstance(test, ast.Attribute) and test.attr == "TYPE_CHECKING")


def bam_key(path: Path) -> str:
    """A key for the BAM at `path` without reading its records: "<size>-<mtime_ns>-<sha256 of the header>".
    The header (reference names and lengths, read groups, program lines) is read through pysam; the alignments
    are not, so the key is cheap for a multi-Gb file but does not see an edit that keeps the size and the
    modification time."""
    import pysam

    path = Path(path)
    st = path.stat()
    with pysam.AlignmentFile(str(path), "rb", check_sq=False) as bam:
        header = str(bam.header)
    return f"{st.st_size}-{st.st_mtime_ns}-{hashlib.sha256(header.encode()).hexdigest()}"


def tool_version(exe: str) -> str | None:
    """`exe`'s version: the first line it prints for its version arguments (constants.TOOL_VERSION_ARGS, else
    --version); None when it is not on PATH (or does not answer in 10 s). Cached per binary file."""
    path = shutil.which(exe)
    if path is None:
        return None
    real = os.path.realpath(path)
    try:
        st = os.stat(real)
    except OSError:
        return None
    args = constants.TOOL_VERSION_ARGS.get(exe, ["--version"])
    return _tool_version(path, tuple(args), st.st_size, st.st_mtime_ns)


@functools.cache
def _tool_version(path: str, args: tuple[str, ...], size: int, mtime_ns: int) -> str | None:
    """tool_version of the executable `path` (its binary's `size` and `mtime_ns` key the cache)."""
    try:
        r = subprocess.run([path, *args], capture_output=True, text=True, timeout=constants.TOOL_VERSION_TIMEOUT,
                           check=False,
                           stdin=subprocess.DEVNULL)
    except (OSError, subprocess.SubprocessError):
        return None
    lines = (r.stdout + "\n" + r.stderr).strip().splitlines()
    return lines[0].strip() if lines else None


def sha256_file(path: Path) -> str:
    """sha256 (hex) of the file at `path`, streamed in CHECKSUM_BLOCK pieces, never loaded whole."""
    h = hashlib.sha256()
    buf = memoryview(bytearray(constants.CHECKSUM_BLOCK))
    with open(path, "rb") as fh:
        while n := fh.readinto(buf):
            h.update(buf[:n])
    return h.hexdigest()


def file_checksum(path: Path, cache: Path) -> str:
    """sha256 (hex) of the file at `path`, cached in the JSON file `cache` by (absolute path, size, mtime_ns), and
    kept by every open `keyed_checksums` collection.

    The file is streamed, never loaded whole. A damaged cache is rebuilt, with a warning. Each cache write is
    atomic; concurrent writers may drop each other's new entries, which only costs a recomputation.
    """
    path, cache = Path(path).resolve(), Path(cache)
    st = path.stat()
    entries = _read_json(cache)
    hit = entries.get(str(path))
    if isinstance(hit, dict) and hit.get("size") == st.st_size and hit.get("mtime_ns") == st.st_mtime_ns:
        sha = hit["sha256"]
    else:
        sha = sha256_file(path)
        entries[str(path)] = {"size": st.st_size, "mtime_ns": st.st_mtime_ns, "sha256": sha}
        _write_json(cache, entries)
    _keep(str(path), sha)
    return sha


_KEYED: list[dict[str, str]] = []  # the open keyed_checksums collections of this process


@contextlib.contextmanager
def keyed_checksums() -> Iterator[dict[str, str]]:
    """Collect the checksums the stages key on while the block runs: every sha256 `file_checksum` returns in this
    process, cache hits included, by absolute path (symlinks resolved). This is a run's record of its inputs as its
    stages read them, whatever the shared checksum cache holds later. A file whose checksum changes inside the block
    (it changed between two stages) keeps the latest, with one warning naming it."""
    keyed: dict[str, str] = {}
    _KEYED.append(keyed)
    try:
        yield keyed
    finally:
        _KEYED[:] = [k for k in _KEYED if k is not keyed]


def _keep(path: str, sha: str) -> None:
    """Record `sha` as the checksum of the file `path` in every open keyed_checksums collection; one warning when it
    replaces another."""
    changed = False
    for keyed in tuple(_KEYED):
        changed |= keyed.get(path, sha) != sha
        keyed[path] = sha
    if changed:
        log.warning("%s changed during the run: its checksum differs between two stages, which may have read different "
                    "versions of it; rerun once it stays as it is", path)


def stage_key(stage: str, inputs: dict[str, str], params: dict, fingerprint: str) -> str:
    """sha256 (hex) of the canonical JSON (sorted keys) of the stage, its inputs, its params, the version and
    the fingerprint of its code (`code_fingerprint(unit)`).

    `inputs` maps input names to checksums; `params` holds JSON values (Path and numpy scalars converted),
    among them the versions of the external tools the stage runs.
    """
    payload = {"stage": stage, "inputs": inputs, "params": params, "version": __version__, "code": fingerprint}
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=_jsonable)
    return hashlib.sha256(blob.encode()).hexdigest()


def is_done(marker: Path, key: str) -> bool:
    """True if `marker` records `key` and every output it lists still exists with its recorded size."""
    rec = _intact(marker)
    return rec is not None and rec.get("key") == key


def done_outputs(marker: Path) -> list[Path] | None:
    """The outputs `marker` lists (resolved against its directory), whatever its key, while every one still exists
    with its recorded size; None without a readable marker or when an output is gone or changed size."""
    rec = _intact(marker)
    return None if rec is None else [Path(marker).parent / out["path"] for out in rec["outputs"]]


def _intact(marker: Path) -> dict | None:
    """The marker's record when it is readable and every output it lists still exists with its recorded size."""
    marker = Path(marker)
    try:
        with open(marker) as fh:
            rec = json.load(fh)
        if all((marker.parent / out["path"]).stat().st_size == out["size"] for out in rec["outputs"]):
            return rec
    except (OSError, ValueError, KeyError, TypeError, AttributeError):  # absent, unreadable or malformed: not done
        pass
    return None


def mark_done(marker: Path, key: str, outputs: list[Path], meta: dict | None = None) -> None:
    """Record that the stage keyed `key` wrote `outputs` (call it last, after every output is complete), with
    the stage's `meta` (JSON values, such as its code fingerprint; read back by `done_meta`).

    Written atomically (temp file + rename). Output paths are stored relative to the marker's directory, so a
    work directory stays valid when moved. A missing output is an error: nothing is written.
    """
    marker = Path(marker)
    base = os.path.abspath(marker.parent)
    recs = []
    for out in outputs:
        try:
            size = Path(out).stat().st_size
        except FileNotFoundError:
            raise FileNotFoundError(f"stage output {out} is missing; {marker} not written") from None
        recs.append({"path": os.path.relpath(os.path.abspath(out), base), "size": size})
    _write_json(marker, {"key": key, "outputs": recs, "version": __version__, "meta": meta or {},
                         "finished": datetime.now(UTC).isoformat(timespec="seconds")})


def done_meta(marker: Path) -> dict | None:
    """The meta mark_done recorded beside the key ({} if none); None without a readable marker."""
    try:
        with open(marker) as fh:
            meta = json.load(fh).get("meta") or {}
        return meta if isinstance(meta, dict) else None
    except (OSError, ValueError, AttributeError):
        return None


def clear_done(marker: Path) -> None:
    """Remove `marker` before a stage rewrites its outputs.

    Otherwise a run interrupted mid-write would leave the old marker beside new, partial outputs, trusted
    again once the old key returns.
    """
    Path(marker).unlink(missing_ok=True)


def _jsonable(obj):
    if isinstance(obj, Path):
        return str(obj)
    if hasattr(obj, "item") and hasattr(obj, "dtype"):  # numpy scalar
        return obj.item()
    raise TypeError(f"stage parameter of type {type(obj).__name__} is not JSON-serializable")


def _read_json(path: Path) -> dict:
    try:
        with open(path) as fh:
            data = json.load(fh)
    except FileNotFoundError:
        return {}
    except (OSError, ValueError) as exc:
        log.warning("checksum cache %s is unreadable (%s); rebuilding it", path, exc)
        return {}
    if not isinstance(data, dict):
        log.warning("checksum cache %s is malformed; rebuilding it", path)
        return {}
    return data


def marker(path: Path) -> Path:
    """The done marker of the file `path`: beside it, named <name>.done.json."""
    return path.with_name(path.name + ".done.json")


@contextlib.contextmanager
def atomic_path(path: Path) -> Iterator[Path]:
    """A temporary path beside `path`, one per process and thread, renamed onto `path` when the block ends
    without an error and removed otherwise: a reader never sees a partial file. The caller writes the file."""
    tmp = path.with_name(f".{path.name}.{os.getpid()}.{threading.get_ident()}.tmp")
    tmp.unlink(missing_ok=True)
    try:
        yield tmp
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def atomic_write(path: Path, data: bytes | str) -> None:
    """Write `data` (text as UTF-8) to `path` through a temporary file renamed into place."""
    with atomic_path(path) as tmp:
        if isinstance(data, str):
            tmp.write_text(data, encoding="utf-8")
        else:
            tmp.write_bytes(data)


def _write_json(path: Path, data: dict) -> None:
    """Write `data` as JSON to `path` atomically (atomic_write)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    atomic_write(path, json.dumps(data, indent=1, sort_keys=True) + "\n")
