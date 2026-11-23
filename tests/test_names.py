"""Names users do not know stay out of what users read: the documents, and every string the package prints or writes
(its docstrings and comments may name the program a rule came from; legacy/ holds those programs)."""
import ast
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
# the four programs the package was built from and its simulator, the package's earlier name, and LTR_peel's
OLD_NAMES = re.compile(r"LTRpav|ltrpav|LTRmix|ltrmix|LTRkin|ltrkin|ltrharvest_phasing|simallo|LTRtrace|ltrtrace|LTRquest")
DOCUMENTS = ("README.md", "CHANGELOG.md", "CITATION.cff", "pyproject.toml", "environment.yml", "docs", "examples")


def _docstrings(tree: ast.AST) -> set[int]:
    """The ids of the docstring nodes of a module, its classes and its functions."""
    found = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Module | ast.ClassDef | ast.FunctionDef | ast.AsyncFunctionDef) and node.body:
            first = node.body[0]
            if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) and isinstance(first.value.value, str):
                found.add(id(first.value))
    return found


def test_documents_name_no_old_program():
    files = []
    for name in DOCUMENTS:
        path = REPO / name
        files += sorted(p for p in path.rglob("*") if p.is_file()) if path.is_dir() else [path]
    assert files
    hits = [f"{p.relative_to(REPO)}:{i}: {m.group()}" for p in files
            for i, line in enumerate(p.read_text().splitlines(), 1) for m in OLD_NAMES.finditer(line)]
    assert not hits, hits


def test_package_strings_name_no_old_program(src_root):
    hits = []
    for path in sorted(src_root.rglob("*.py")):
        tree = ast.parse(path.read_text())
        docstrings = _docstrings(tree)
        hits += [f"{path.relative_to(src_root)}:{node.lineno}: {OLD_NAMES.search(node.value).group()}"
                 for node in ast.walk(tree)
                 if isinstance(node, ast.Constant) and isinstance(node.value, str) and id(node) not in docstrings
                 and OLD_NAMES.search(node.value)]
    assert not hits, hits
