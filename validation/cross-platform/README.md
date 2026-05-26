# Cross-platform reproducibility

Bit-identical trajectory reproduction of the v0.1 `apsis` federation portfolio across Windows and Linux on x86_64. The claim and the per-`pow`-implementation analysis are recorded in [`paper/notebooks/2026-05-20-cross-platform-determinism.md`](../../paper/notebooks/2026-05-20-cross-platform-determinism.md) and [`paper/notebooks/2026-05-22-controller-pow-implementations.md`](../../paper/notebooks/2026-05-22-controller-pow-implementations.md).

## Layout

| Path                  | Purpose                                                                                                                  |
| --------------------- | ------------------------------------------------------------------------------------------------------------------------ |
| `run_linux_side.sh`   | Executes the seven parity scenarios + Mercury 1PN long-horizon under the same toolchain pin as the Windows reference run |
| `compare.py`          | Loads paired CSV outputs from `windows/` and `linux/`, reports per-column ULP distance + the Mercury rate diff           |
| `windows/`            | Captured Windows reference outputs (committed): seven `*.csv` parity scenarios, `mercury_perihelion.txt`, `meta.txt`     |
| `linux/`              | Captured Linux outputs from the same source commit and lockfile (gitignored; produced by `run_linux_side.sh`)            |

## Scenarios

The seven scenarios cover the full v0.1 federation:

| CSV / file                                  | Source `cargo` example                                          | Crate exercised                  |
| ------------------------------------------- | --------------------------------------------------------------- | -------------------------------- |
| `kepler.csv`                                | `rebound_parity_kepler`                                         | `apsis` (IAS15)                  |
| `figure8.csv`                               | `rebound_parity_figure8`                                        | `apsis` (IAS15)                  |
| `pythagorean.csv`                           | `rebound_parity_pythagorean`                                    | `apsis` (IAS15)                  |
| `retrograde.csv`                            | `rebound_parity_retrograde`                                     | `apsis` (IAS15)                  |
| `mercurius_outer_solar.csv`                 | `rebound_parity_mercurius`                                      | `apsis` (MERCURIUS hybrid)       |
| `whfast_outer_solar.csv`                    | `whfast_outer_solar`                                            | `apsis` (Wisdom–Holman)          |
| `central_observable_inversion_long.csv`     | `central_observable_inversion_long`                             | `apsis-central` + `apsis` (IAS15) |
| `mercury_perihelion.txt`                    | `mercury_perihelion`                                            | `apsis-1pn` + `apsis` (IAS15)    |

## Workflow

The reference run is captured on a Windows host; the Linux side runs against the same source commit and `Cargo.lock`; ULP-level comparison happens on the developer machine where both captures are available.

### 1. Capture the Windows reference (one-time, then commit)

Each cargo example writes its CSV into `validation/cross-platform/windows/` by default. `mercury_perihelion` prints to stdout and is redirected. After the runs, `meta.txt` records the toolchain version and the `Cargo.lock` hash.

```bash
cargo run --release --example rebound_parity_kepler -p apsis --quiet
cargo run --release --example rebound_parity_figure8 -p apsis --quiet
cargo run --release --example rebound_parity_pythagorean -p apsis --quiet
cargo run --release --example rebound_parity_retrograde -p apsis --quiet
cargo run --release --example rebound_parity_mercurius -p apsis --quiet
cargo run --release --example whfast_outer_solar -p apsis --quiet
cargo run --release --example central_observable_inversion_long -p apsis-central --quiet
cargo run --release --example mercury_perihelion -p apsis-1pn --quiet \
  > validation/cross-platform/windows/mercury_perihelion.txt
```

### 2. Run the Linux side at the same commit

Any Linux x86_64 host with `rustup` works (e.g. an EC2 spot instance). The side script produces a tarball at `/tmp/apsis-xplat-linux.tar.gz`; `scp` it back to the developer machine and extract into `validation/cross-platform/linux/`.

```bash
git clone <repo> apsis && cd apsis
git checkout <commit-recorded-in-windows-meta>
cargo build --release --workspace
bash validation/cross-platform/run_linux_side.sh
```

### 3. Compare

```bash
python validation/cross-platform/compare.py
```

`--json out.json` writes the per-scenario verdict for downstream notebook ingestion.

## Output

`compare.py` emits a Markdown table to stdout with per-scenario columns:

- bit-equal rows / total rows
- maximum ULP distance across any single column
- file-size delta in bytes

A scenario where every column agrees at 0 ULP is the success criterion. The single-scalar Mercury 1PN gate (parsed from `mercury_perihelion.txt`) reports the difference in `Δω` rate against the same GR prediction value, in fractional units of the predicted rate.

## Why this directory commits captured outputs

The seven parity CSVs and `mercury_perihelion.txt` under `windows/` are committed because the cross-platform bit-equal claim is verifiable from the captured files plus the source commit + `Cargo.lock`. A reviewer can rerun `run_linux_side.sh` on any Linux x86_64 host with the pinned toolchain, place the output in `linux/`, and reproduce the comparison without coordinating with the original author. `linux/` is gitignored — it is reconstructed locally per reviewer host.
