"""Open an Apsis Record and print its provenance + event log.

Demonstrates the read-only Python API. The header is the canonical TOML
provenance block; this example prints it verbatim so the script runs on
Python 3.10+ without a TOML parser (``tomllib`` is 3.11+, ``tomli`` is
a third-party dep). Scripts that need structured access should add
``tomllib``/``tomli`` themselves.

Run::

    python examples/inspect_record.py mercury_1pn.apsis
"""

from __future__ import annotations

import sys

import apsis


def main() -> None:
    if len(sys.argv) != 2:
        print("usage: inspect_record.py <path-to-.apsis>", file=sys.stderr)
        sys.exit(2)

    rec = apsis.Record(sys.argv[1])

    print("=== header (TOML) ===")
    print(rec.header)

    print("=== events ===")
    events = rec.events()
    print(f"total: {len(events)}")
    for ev in events[:10]:
        print(" ", ev)

    print("=== snapshots ===")
    print(f"count: {rec.snapshot_count()}")


if __name__ == "__main__":
    main()
