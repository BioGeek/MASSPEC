"""The code stays PEP 8, without single-character names or lambdas."""
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def test_no_single_character_names_or_lambdas():
    spec = importlib.util.spec_from_file_location(
        "check_style", ROOT / "scripts" / "check_style.py")
    check_style = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(check_style)
    assert check_style.main() == 0
