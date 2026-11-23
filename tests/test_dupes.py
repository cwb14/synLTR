"""No two functions of 5+ statements may have the same AST once identifiers, constants and docstrings are ignored."""

import ast
from pathlib import Path

MIN_STATEMENTS = 5


class _Normalize(ast.NodeTransformer):
    """Replace every identifier and constant by a placeholder."""

    def visit_Name(self, node):
        return ast.copy_location(ast.Name(id="_", ctx=node.ctx), node)

    def visit_Attribute(self, node):
        self.generic_visit(node)
        node.attr = "_"
        return node

    def visit_arg(self, node):
        node.arg = "_"
        node.annotation = None
        return node

    def visit_keyword(self, node):
        self.generic_visit(node)
        node.arg = "_" if node.arg else None
        return node

    def visit_Constant(self, node):
        return ast.copy_location(ast.Constant(value=0), node)

    def visit_FunctionDef(self, node):
        self.generic_visit(node)
        node.name = "_"
        node.returns = None
        node.decorator_list = []
        return node

    visit_AsyncFunctionDef = visit_FunctionDef


def _strip_docstring(fn):
    body = fn.body
    if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) \
            and isinstance(body[0].value.value, str):
        fn.body = body[1:] or [ast.Pass()]


def duplicate_functions(root: Path) -> list[tuple[str, str]]:
    """Pairs 'file:function' whose normalized ASTs are equal, for functions of MIN_STATEMENTS or more."""
    seen: dict[str, str] = {}
    pairs = []
    for path in sorted(root.rglob("*.py")):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            fn = ast.parse(ast.unparse(node)).body[0]
            _strip_docstring(fn)
            n = sum(isinstance(x, ast.stmt) for x in ast.walk(fn)) - 1
            if n < MIN_STATEMENTS:
                continue
            key = ast.dump(_Normalize().visit(fn))
            where = f"{path.relative_to(root)}:{node.name}"
            if key in seen:
                pairs.append((seen[key], where))
            else:
                seen[key] = where
    return pairs


def test_package_has_no_duplicate_functions(src_root):
    assert duplicate_functions(src_root) == []


BODY = """
def {name}({arg}):
    '''{doc}'''
    total = {c}
    for item in {arg}:
        if item > {c}:
            total += item
        else:
            total -= 1
    return total
"""


def test_dupes_catch_renamed_copy(tmp_path):
    (tmp_path / "a.py").write_text(BODY.format(name="f", arg="xs", doc="one", c=0))
    (tmp_path / "b.py").write_text(BODY.format(name="g", arg="values", doc="two", c=7))
    assert duplicate_functions(tmp_path) == [("a.py:f", "b.py:g")]


def test_dupes_ignore_short_and_different_functions(tmp_path):
    (tmp_path / "a.py").write_text("def f(x):\n    y = x\n    return y\n")
    (tmp_path / "b.py").write_text("def g(x):\n    y = x\n    return y\n")
    (tmp_path / "c.py").write_text(BODY.format(name="h", arg="xs", doc="", c=0).replace("-= 1", "*= 2"))
    (tmp_path / "d.py").write_text(BODY.format(name="f", arg="xs", doc="", c=0))
    assert duplicate_functions(tmp_path) == []
