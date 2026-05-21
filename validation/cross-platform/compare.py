"""Cross-platform parity comparator: Windows reference vs Linux EC2 run.

Reads paired outputs from ``validation/cross-platform/{windows,linux}/`` and
emits a per-scenario table reporting how many ULPs each CSV column drifted
between platforms. The single-number Mercury 1PN gate is parsed from the
captured stdout and compared as a fraction of the expected GR rate.

Run from the repo root:

    python validation/cross-platform/compare.py

Outputs a Markdown table to stdout and (if ``--json out.json``) also a
structured JSON dump for downstream notebook ingestion.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path

ROOT = Path(__file__).parent

PARITY_SCENARIOS = ["kepler", "figure8", "pythagorean", "retrograde"]


def ulp_distance(a: float, b: float) -> float:
    """Distance between two f64 values in ULPs of max(|a|, |b|).

    Returns 0.0 if both are bit-equal (including same-signed zeros and
    same NaN representations). Returns ``math.inf`` if exactly one is
    NaN / inf / sign-mismatch.
    """
    if a == b:
        return 0.0
    if math.isnan(a) and math.isnan(b):
        return 0.0
    if math.isnan(a) != math.isnan(b):
        return math.inf
    if math.isinf(a) or math.isinf(b):
        return math.inf
    scale = max(abs(a), abs(b))
    if scale == 0.0:
        return math.inf
    eps = math.ldexp(1.0, -52)  # 2^-52 = f64 ULP at unity
    return abs(a - b) / (scale * eps)


def read_csv(path: Path) -> tuple[list[str], list[list[float]]]:
    """Header line + numeric rows. Skips ``#``-prefixed comment lines and
    blank lines anywhere in the file."""
    header: list[str] | None = None
    rows: list[list[float]] = []
    with path.open() as f:
        for line_no, raw in enumerate(f, start=1):
            raw = raw.rstrip("\n")
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            if header is None:
                header = raw.split(",")
                continue
            try:
                rows.append([float(x) for x in raw.split(",")])
            except ValueError as e:
                raise ValueError(f"{path}:{line_no} parse error: {e}") from e
    if header is None:
        raise ValueError(f"{path}: no header found (file empty or all comments)")
    return header, rows


@dataclass
class ColumnDiff:
    column: str
    n_rows: int
    n_bit_equal: int
    max_ulp: float
    mean_ulp: float
    max_abs_diff: float
    sample_at_max: tuple[float, float] = (0.0, 0.0)


def _empty_columns() -> list[ColumnDiff]:
    return []


@dataclass
class ScenarioDiff:
    name: str
    n_rows_windows: int
    n_rows_linux: int
    columns: list[ColumnDiff] = field(default_factory=_empty_columns)
    parity_ok: bool = True
    error: str | None = None


@dataclass
class MercuryDiff:
    name: str = "mercury_perihelion"
    windows_rate: float | None = None
    linux_rate: float | None = None
    abs_diff: float | None = None
    ppm_of_gr: float | None = None
    error: str | None = None


def diff_csv(name: str, win_path: Path, lin_path: Path) -> ScenarioDiff:
    if not win_path.exists():
        return ScenarioDiff(name, 0, 0, error=f"missing windows file: {win_path.name}")
    if not lin_path.exists():
        return ScenarioDiff(name, 0, 0, error=f"missing linux file: {lin_path.name}")

    win_header, win_rows = read_csv(win_path)
    lin_header, lin_rows = read_csv(lin_path)

    if win_header != lin_header:
        return ScenarioDiff(
            name,
            len(win_rows),
            len(lin_rows),
            error=f"header mismatch:\n  win: {win_header}\n  lin: {lin_header}",
        )

    n = min(len(win_rows), len(lin_rows))
    if len(win_rows) != len(lin_rows):
        sys.stderr.write(
            f"warn: {name}: row count differs "
            f"(win={len(win_rows)} lin={len(lin_rows)}, comparing first {n})\n"
        )

    result = ScenarioDiff(name, len(win_rows), len(lin_rows))
    for col_idx, col_name in enumerate(win_header):
        bit_equal = 0
        max_ulp = 0.0
        ulp_sum = 0.0
        max_abs = 0.0
        sample = (0.0, 0.0)
        for row_idx in range(n):
            a = win_rows[row_idx][col_idx]
            b = lin_rows[row_idx][col_idx]
            d = ulp_distance(a, b)
            if d == 0.0:
                bit_equal += 1
            if d > max_ulp:
                max_ulp = d
                sample = (a, b)
            ulp_sum += 0.0 if math.isinf(d) else d
            abs_diff = abs(a - b)
            if abs_diff > max_abs:
                max_abs = abs_diff
        result.columns.append(
            ColumnDiff(
                column=col_name,
                n_rows=n,
                n_bit_equal=bit_equal,
                max_ulp=max_ulp,
                mean_ulp=ulp_sum / n if n else 0.0,
                max_abs_diff=max_abs,
                sample_at_max=sample,
            )
        )
    return result


_RATE_RE = re.compile(r"rate\s*=\s*([+-]?\d+\.?\d*)\s*arcsec/century")


def diff_mercury(win_path: Path, lin_path: Path) -> MercuryDiff:
    if not win_path.exists() or not lin_path.exists():
        return MercuryDiff(error="missing one or both mercury_perihelion.txt files")
    win_match = _RATE_RE.search(win_path.read_text())
    lin_match = _RATE_RE.search(lin_path.read_text())
    if not win_match or not lin_match:
        return MercuryDiff(error="could not parse `rate = X arcsec/century` line")
    win_rate = float(win_match.group(1))
    lin_rate = float(lin_match.group(1))
    abs_diff = abs(win_rate - lin_rate)
    return MercuryDiff(
        windows_rate=win_rate,
        linux_rate=lin_rate,
        abs_diff=abs_diff,
        ppm_of_gr=abs_diff / 43.0 * 1e6,
    )


def render_scenario(d: ScenarioDiff) -> str:
    lines = [f"### {d.name}\n"]
    if d.error:
        lines.append(f"**error:** {d.error}\n")
        return "\n".join(lines)
    if d.n_rows_windows != d.n_rows_linux:
        lines.append(
            f"_row count: win={d.n_rows_windows} lin={d.n_rows_linux}_\n"
        )
    lines.append("| column | bit-equal / total | max ULP | mean ULP | max \\|Δ\\| |")
    lines.append("| --- | --- | --- | --- | --- |")
    for col in d.columns:
        lines.append(
            f"| `{col.column}` | "
            f"{col.n_bit_equal}/{col.n_rows} | "
            f"{col.max_ulp:.3g} | "
            f"{col.mean_ulp:.3g} | "
            f"{col.max_abs_diff:.3e} |"
        )
    return "\n".join(lines)


def render_mercury(m: MercuryDiff) -> str:
    lines = ["### mercury_perihelion (federation operator)\n"]
    if m.error is not None:
        lines.append(f"**error:** {m.error}")
        return "\n".join(lines)
    assert m.windows_rate is not None
    assert m.linux_rate is not None
    assert m.abs_diff is not None
    assert m.ppm_of_gr is not None
    lines.append("| windows rate | linux rate | \\|Δ\\| arcsec/century | ppm of GR (43) |")
    lines.append("| --- | --- | --- | --- |")
    lines.append(
        f"| {m.windows_rate:.6f} | "
        f"{m.linux_rate:.6f} | "
        f"{m.abs_diff:.3e} | "
        f"{m.ppm_of_gr:.3f} |"
    )
    return "\n".join(lines)


def main() -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", type=Path, default=ROOT / "windows",
                    help="first run dir (default: ./windows)")
    ap.add_argument("--b", type=Path, default=ROOT / "linux",
                    help="second run dir (default: ./linux)")
    ap.add_argument("--json", type=Path, help="also write structured JSON to PATH")
    args = ap.parse_args()

    a_dir: Path = args.a
    b_dir: Path = args.b

    print("# Parity\n")
    print(f"A: `{a_dir}`")
    print(f"B: `{b_dir}`\n")

    scenarios: list[ScenarioDiff] = []
    for name in PARITY_SCENARIOS:
        d = diff_csv(name, a_dir / f"{name}.csv", b_dir / f"{name}.csv")
        scenarios.append(d)
        print(render_scenario(d))
        print()

    mercury = diff_mercury(
        a_dir / "mercury_perihelion.txt",
        b_dir / "mercury_perihelion.txt",
    )
    print(render_mercury(mercury))
    print()

    print("---")
    print("ULP definition: distance in 2^-52 units of `max(|a|, |b|)`. "
          "`0` means bit-equal; `inf` means sign mismatch or NaN vs finite.")

    if args.json:
        def _clean(s: ScenarioDiff) -> dict[str, object]:
            d = asdict(s)
            for c in d["columns"]:
                if not math.isfinite(c["max_ulp"]):
                    c["max_ulp"] = None
                if not math.isfinite(c["mean_ulp"]):
                    c["mean_ulp"] = None
            return d

        payload: dict[str, object] = {
            "parity": [_clean(s) for s in scenarios],
            "mercury": asdict(mercury),
        }
        args.json.write_text(json.dumps(payload, indent=2))
        print(f"\n[wrote JSON to {args.json}]", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
