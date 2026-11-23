"""Import layers: a module imports lower ranks, its own unit, or another rank-0 unit (no cycles)."""

import ast
from pathlib import Path

# (unit prefix, rank); the longest matching prefix wins
UNITS = [
    ("base", 0), ("io", 0), ("seq", 0), ("stats", 0), ("plot.style", 0), ("calls.states", 0),
    ("calls.records", 0),
    ("routes", 1), ("routes.library", 1), ("routes.assembly", 1), ("routes.reads", 1),
    ("calls", 2),
    ("analyses", 3), ("analyses.geneflow", 3), ("analyses.subgenomes", 3), ("analyses.source", 3),
    ("plot", 4), ("sim", 4),
    ("cli", 5), ("__main__", 6),
]
ROOT_PACKAGE = "synltr"


def unit_of(name: str) -> tuple[str, int]:
    """The (unit, rank) of a module name relative to the package (e.g. 'calls.vcf' -> ('calls', 2))."""
    best = max((u for u in UNITS if name == u[0] or name.startswith(u[0] + ".")), key=lambda u: len(u[0]),
               default=None)
    if best is None:
        raise AssertionError(f"module {name!r} belongs to no layer; add it to UNITS")
    return best


def module_files(root: Path) -> dict[str, Path]:
    """Map each module name relative to the package ('calls.vcf'; packages by their folder) to its file."""
    out = {}
    for path in sorted(root.rglob("*.py")):
        rel = path.relative_to(root).with_suffix("")
        parts = list(rel.parts)
        if parts[-1] == "__init__":
            parts = parts[:-1]
        if parts:
            out[".".join(parts)] = path
    return out


def _is_type_checking(test: ast.expr) -> bool:
    return (isinstance(test, ast.Name) and test.id == "TYPE_CHECKING") or (
        isinstance(test, ast.Attribute) and test.attr == "TYPE_CHECKING")


def _imports(node: ast.AST):
    """Yield every import node under `node`, skipping the body of `if TYPE_CHECKING:` (its else branch counts)."""
    for child in ast.iter_child_nodes(node):
        if isinstance(child, ast.If) and _is_type_checking(child.test):
            for stmt in child.orelse:
                yield from [stmt] if isinstance(stmt, (ast.Import, ast.ImportFrom)) else _imports(stmt)
        elif isinstance(child, (ast.Import, ast.ImportFrom)):
            yield child
        else:
            yield from _imports(child)


def imported_modules(name: str, path: Path, known: set[str]) -> set[str]:
    """The package-internal modules (relative names) that module `name` imports."""
    package = name.split(".") if path.name == "__init__.py" else name.split(".")[:-1]
    found = set()
    for node in _imports(ast.parse(path.read_text())):
        targets = []
        if isinstance(node, ast.Import):
            targets = [a.name for a in node.names]
        else:
            if node.level:
                base = package[: len(package) - (node.level - 1)]
                full = ".".join([ROOT_PACKAGE, *base] + ([node.module] if node.module else []))
            else:
                full = node.module or ""
            targets = []
            for a in node.names:  # `from pkg import name` imports module pkg.name when there is one
                sub = f"{full}.{a.name}"
                targets.append(sub if sub[len(ROOT_PACKAGE) + 1:] in known else full)
        for t in targets:
            if t == ROOT_PACKAGE or not t.startswith(ROOT_PACKAGE + "."):
                continue
            rel = t[len(ROOT_PACKAGE) + 1:]
            if rel in known:
                found.add(rel)
    return found


def layer_violations(root: Path) -> list[str]:
    """Every import that breaks the layer rule, plus any cycle among rank-0 units."""
    files = module_files(root)
    known = set(files)
    bad, rank0_edges = [], {}
    for name, path in files.items():
        unit, rank = unit_of(name)
        for dep in imported_modules(name, path, known):
            dunit, drank = unit_of(dep)
            if dunit == unit or drank < rank:
                continue
            if drank == 0 and rank == 0:
                rank0_edges.setdefault(unit, set()).add(dunit)
                continue
            bad.append(f"{name} (rank {rank}) imports {dep} (rank {drank})")
    # rank-0 units must form a DAG
    state = {}

    def visit(u, trail):
        state[u] = 1
        for v in rank0_edges.get(u, ()):
            if state.get(v) == 1:
                bad.append("rank-0 import cycle: " + " -> ".join([*trail, u, v]))
            elif v not in state:
                visit(v, [*trail, u])
        state[u] = 2

    for u in list(rank0_edges):
        if u not in state:
            visit(u, [])
    return bad


def test_package_respects_layers(src_root):
    assert layer_violations(src_root) == []


def _plant(tmp_path, files: dict[str, str]) -> Path:
    root = tmp_path / ROOT_PACKAGE
    for rel, text in files.items():
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
    return root


def test_layers_catch_upward_import(tmp_path):
    root = _plant(tmp_path, {"__init__.py": "", "base/__init__.py": "", "base/x.py": "from synltr import cli\n",
                             "cli.py": ""})
    assert layer_violations(root) == ["base.x (rank 0) imports cli (rank 5)"]


def test_layers_catch_route_importing_route(tmp_path):
    root = _plant(tmp_path, {"routes/library/__init__.py": "", "routes/reads/__init__.py": "",
                             "routes/reads/a.py": "from synltr.routes.library import b\n",
                             "routes/library/b.py": ""})
    assert layer_violations(root) == ["routes.reads.a (rank 1) imports routes.library.b (rank 1)"]


def test_layers_catch_rank0_cycle_but_allow_dag(tmp_path):
    cyc = _plant(tmp_path / "c", {"io/a.py": "import synltr.seq.b\n", "seq/b.py": "from synltr import io\n",
                                  "io/__init__.py": "", "seq/__init__.py": ""})
    assert any("cycle" in v for v in layer_violations(cyc))
    dag = _plant(tmp_path / "d", {"io/a.py": "import synltr.seq.b\n", "seq/b.py": "",
                                  "io/__init__.py": "", "seq/__init__.py": ""})
    assert layer_violations(dag) == []


def test_layers_ignore_type_checking_imports(tmp_path):
    src = "from typing import TYPE_CHECKING\nif TYPE_CHECKING:\n    from synltr import cli\n"
    assert layer_violations(_plant(tmp_path, {"base/x.py": src, "base/__init__.py": "", "cli.py": ""})) == []


def test_layers_resolve_relative_imports(tmp_path):
    root = _plant(tmp_path, {"base/__init__.py": "", "base/x.py": "from .. import cli\n", "cli.py": ""})
    assert layer_violations(root) == ["base.x (rank 0) imports cli (rank 5)"]


def test_layers_scan_type_checking_else_branch(tmp_path):
    src = ("from typing import TYPE_CHECKING\nif TYPE_CHECKING:\n    pass\nelse:\n    from synltr import cli\n"
           "if TYPE_CHECKING:\n    pass\nelse:\n    if True:\n        import synltr.cli\n")
    root = _plant(tmp_path, {"base/x.py": src, "base/__init__.py": "", "cli.py": ""})
    assert layer_violations(root) == ["base.x (rank 0) imports cli (rank 5)"]
