"""Guard: comparator and IC generator must not import apsis (oracle independence).

rebound_side.py and smoke_pair.py are excluded -- they are the REBOUND-side
driver and integration harness; importing rebound there is intentional and
importing apsis there would be a separate concern caught by integration tests.
"""
import ast
from pathlib import Path

HERE = Path(__file__).resolve().parent

_CHECKED = {"compare.py", "generate_ics.py"}


def test_no_apsis_import():
    offenders = []
    for py in sorted(HERE.glob("*.py")):
        if py.name not in _CHECKED:
            continue
        tree = ast.parse(py.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            mods = []
            if isinstance(node, ast.Import):
                mods = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                mods = [node.module or ""]
            if any(m.split(".")[0] == "apsis" for m in mods):
                offenders.append(py.name)
    assert not offenders, f"oracle files must not import apsis: {offenders}"


if __name__ == "__main__":
    test_no_apsis_import()
    print("ok")
