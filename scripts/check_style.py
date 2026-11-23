"""Fail on what PEP 8 tooling misses: single-character names and lambdas.

    uv run python scripts/check_style.py

Checks every Python file under src/, scripts/, tests/ and notebooks/, and
the code cells that notebooks/build_notebook.py writes into the notebook.
"""
import ast
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def names(tree):
    """Every name the code BINDS, with its line.

    That is: targets, arguments, loop and comprehension variables, import
    aliases, functions and classes.
    """
    for node in ast.walk(tree):
        if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store):
            yield node.id, node.lineno
        elif isinstance(node, ast.arg):
            yield node.arg, node.lineno
        elif isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            yield node.name, node.lineno
        elif isinstance(node, ast.alias) and node.asname:
            yield node.asname, node.lineno
        elif isinstance(node, ast.ExceptHandler) and node.name:
            yield node.name, node.lineno


def problems(source, label):
    tree = ast.parse(source)
    out = []
    for name, line in names(tree):
        if len(name) == 1 and name != "_":
            out.append(f"{label}:{line}: single-character name '{name}'")
    for node in ast.walk(tree):
        if isinstance(node, ast.Lambda):
            out.append(f"{label}:{node.lineno}: lambda")
    return out


def notebook_cells():
    """The code cells build_notebook.py defines, as (index, source)."""
    tree = ast.parse((ROOT / "notebooks" / "build_notebook.py").read_text())
    found = []
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id == "add_code" and node.args
                and isinstance(node.args[0], ast.Constant)):
            found.append((node.lineno, node.args[0].value))
    return found


def main():
    report = []
    for folder in ("src", "scripts", "tests", "notebooks"):
        for path in sorted((ROOT / folder).rglob("*.py")):
            report += problems(path.read_text(), str(path.relative_to(ROOT)))
    for line, cell in notebook_cells():
        report += problems(cell, f"notebook cell at build_notebook.py:{line}")
    print("\n".join(report) or "no single-character names, no lambdas")
    return 1 if report else 0


if __name__ == "__main__":
    sys.exit(main())
