"""Open an Apsis Record and print its provenance + event log.

Demonstrates the read-only Python API. The header is the canonical
TOML provenance block; parse with `tomllib` for structured access.

Run::

    python examples/inspect_record.py mercury_1pn.apsis
"""

from __future__ import annotations

import sys
import tomllib

import apsis


def main() -> None:
    if len(sys.argv) != 2:
        print("usage: inspect_record.py <path-to-.apsis>", file=sys.stderr)
        sys.exit(2)

    rec = apsis.Record(sys.argv[1])
    header = tomllib.loads(rec.header)

    print(f"apsis      = {header['apsis']['version']}")
    print(f"git_sha    = {header['apsis']['git_sha']}")
    print(f"created    = {header['apsis']['created_utc']}")
    print(f"seed       = {header['reproducibility']['seed']}")
    print(f"lock_hash  = {header['reproducibility']['cargo_lock_blake3'][:16]}…")
    print(f"integrator = {header['integrator']['kind']}")
    print(f"kernel     = {header['kernel']['variant']}")
    print(f"n_bodies   = {header['bodies']['count']}")

    operators = header.get("operators", [])
    if operators:
        print(f"operators  = {len(operators)}")
        for op in operators:
            req = op.get("requirements", {})
            print(
                f"  {op['name']:24s} v{op['version']:8s} "
                f"exactness={req.get('kernel_exactness', '-'):8s} "
                f"continuity={req.get('kernel_continuity', '-')}"
            )

    events = rec.events()
    print(f"events     = {len(events)}")
    for ev in events[:10]:
        print(" ", ev)

    print(f"snapshots  = {rec.snapshot_count()}")


if __name__ == "__main__":
    main()
