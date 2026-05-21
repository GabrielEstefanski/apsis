#!/usr/bin/env bash
#
# Cross-platform parity — Linux side of the apsis lockfile-pins-physics gate.
#
# Run inside the apsis repo checked out to the same commit as the Windows
# reference (see ../windows/meta.txt). Assumes the Rust toolchain is
# installed (rustup) and the workspace already builds (`cargo build --release
# --workspace`).
#
# Produces:
#   /tmp/apsis-xplat/linux/
#     ├── meta.txt
#     ├── kepler.csv
#     ├── figure8.csv
#     ├── pythagorean.csv
#     ├── retrograde.csv
#     └── mercury_perihelion.txt
#
# And a tarball: /tmp/apsis-xplat-linux.tar.gz
#
# scp the tarball back to your dev machine and extract into
# validation/cross-platform/linux/.
#
# Usage:
#   cd ~/apsis
#   bash validation/cross-platform/run_linux_side.sh

set -eu

OUT_DIR="${OUT_DIR:-/tmp/apsis-xplat/linux}"
mkdir -p "$OUT_DIR"

echo "── meta ───────────────────────────────────────"

{
    echo "=== rustc ==="
    rustc --version
    echo
    echo "=== cargo ==="
    cargo --version
    echo
    echo "=== OS ==="
    if [ -r /etc/os-release ]; then
        . /etc/os-release
        echo "$PRETTY_NAME"
    else
        uname -srm
    fi
    echo
    echo "=== Kernel ==="
    uname -a
    echo
    echo "=== CPU ==="
    if command -v lscpu >/dev/null 2>&1; then
        lscpu | grep -E "Architecture|Model name|CPU\(s\):|CPU MHz|CPU max MHz|Flags" | head -20
    else
        grep -m1 "model name" /proc/cpuinfo
        grep -c ^processor /proc/cpuinfo | awk '{print "CPU count:", $1}'
    fi
    echo
    echo "=== glibc ==="
    ldd --version | head -1
    echo
    echo "=== Cargo.lock SHA256 ==="
    sha256sum Cargo.lock | awk '{print toupper($1)}'
    echo
    echo "=== git ==="
    git rev-parse HEAD
} > "$OUT_DIR/meta.txt"

cat "$OUT_DIR/meta.txt"

echo
echo "── building (release) ─────────────────────────"
cargo build --release --workspace --quiet

echo
echo "── parity scenarios ───────────────────────────"

run_scenario() {
    local name="$1"
    local example="$2"
    local crate="$3"
    local out="$OUT_DIR/${name}.csv"
    echo "  $name → $out"
    cargo run --release --example "$example" -p "$crate" --quiet -- --output "$out"
}

run_scenario "kepler"       "rebound_parity_kepler"       "apsis"
run_scenario "figure8"      "rebound_parity_figure8"      "apsis"
run_scenario "pythagorean"  "rebound_parity_pythagorean"  "apsis"
run_scenario "retrograde"   "rebound_parity_retrograde"   "apsis"

echo
echo "── mercury 1PN ────────────────────────────────"
cargo run --release --example mercury_perihelion -p apsis-1pn --quiet \
    > "$OUT_DIR/mercury_perihelion.txt"
tail -5 "$OUT_DIR/mercury_perihelion.txt"

echo
echo "── tarball ────────────────────────────────────"
TARBALL="/tmp/apsis-xplat-linux.tar.gz"
tar -czf "$TARBALL" -C "$(dirname "$OUT_DIR")" "$(basename "$OUT_DIR")"
ls -lh "$TARBALL"

echo
echo "Done. scp back with:"
echo "  scp -i <key.pem> ubuntu@<public-ip>:$TARBALL ."
