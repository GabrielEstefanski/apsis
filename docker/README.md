# Validation container

A reproducible Linux environment for apsis's validation harnesses and figure
pipeline: the Rust toolchain plus a Python environment with the figure and
REBOUNDx-parity dependencies. Its reason to exist is that **REBOUNDx does not
build on Windows/MSVC** (C99 VLAs in `gr_full.c`), so the REBOUNDx-parity
checks and any reboundx-dependent work need a Linux toolchain; the container
makes that reproducible on any host with Docker.

## Scope (important)

This is the reproducible **validation / figures** environment. It is **not**
the platform for the §3.5 cross-platform bit-exactness claim — that claim is
about native Windows vs native Linux (`libm` giving bit-identical trajectories
across real platforms), and an all-Linux container would mask exactly what it
demonstrates. Run the §3.5 cross-platform protocol natively; use this
container for everything else (figures, REBOUNDx parity, notebook scripts).

## Build

```bash
docker build -f docker/Dockerfile -t apsis-validation .
```

Pins: Ubuntu 24.04, Rust 1.89.0 (the workspace `rust-version`), and the Python
deps in `docker/requirements.txt` (rebound 4.6.0 + reboundx 4.6.2, compiled
fresh; numpy/matplotlib pinned to the figure-producing versions).

## Use

The repo is bind-mounted at run time, so the source is always your current
checkout (nothing is copied into the image):

```bash
docker run --rm -it -v "$PWD:/apsis" apsis-validation
```

Then, inside the container:

```bash
make figures                                          # re-render paper figures
cd validation/reboundx-parity/gr-mercury && python run.py   # REBOUNDx gr parity
cargo test --release -p apsis-1pn --tests -- --ignored      # release gates
```

Or run a single command without an interactive shell:

```bash
docker run --rm -v "$PWD:/apsis" apsis-validation make figures
```

## File ownership on Linux hosts

The container runs as root and the repo is bind-mounted, so on a **Linux
host** anything written back into the tree — `make figures` (PDFs under
`paper/figures/`), `run.py` (`out/*.csv`) — lands root-owned, and you need
`sudo` to edit or `git checkout` it afterwards. (On macOS and Windows/Docker
Desktop the bind-mount remaps ownership, so this does not bite.)

To keep outputs owned by you, run as your host UID:

```bash
docker run --rm -it --user "$(id -u):$(id -g)" -v "$PWD:/apsis" apsis-validation
```

The Python figure/parity flows need only `--user`. A non-root UID cannot write
cargo's default `CARGO_HOME` (`/root/.cargo`), so for `cargo` flows also point
it at a writable path:

```bash
docker run --rm -it --user "$(id -u):$(id -g)" -e CARGO_HOME=/tmp/.cargo \
  -v "$PWD:/apsis" apsis-validation \
  cargo test --release -p apsis-1pn --tests -- --ignored
```

## Notes

- **git provenance:** the image sets `safe.directory` for `/apsis` (`git
  config --system`) so the build scripts' `git rev-parse` works on the
  bind-mounted repo despite the root/host-UID ownership split. Without it git
  refuses the repo ("dubious ownership") and records ship an empty
  `apsis.git_sha`.
- **reboundx install:** the Dockerfile uses `pip install --no-cache-dir` so
  reboundx is compiled against its co-installed rebound. Reusing a cached
  wheel bakes a stale `librebound` RPATH that fails to load — do not drop
  `--no-cache-dir`.
- **rebound version:** reboundx 4.6.2 requires rebound 4.6.0 (it uninstalls
  rebound 5 if present). When reboundx ships a rebound-5–compatible release,
  bump both in `docker/requirements.txt`.
- **Pinning:** `scipy`/`pandas` are floored, not pinned exactly; run
  `pip freeze` in the built image and pin them for full reproducibility.
