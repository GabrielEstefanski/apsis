"""Type stubs for the public ``apsis`` package surface.

The stubs here are written to match the runtime API exposed by
``python/apsis/__init__.py`` and re-exported from the Rust extension
module ``apsis._native``. They drive IDE autocompletion and
``mypy --strict`` type checking; CI gates a release on the stubs and
the runtime agreeing.

Each subsequent PR that adds a class or free function to the Rust
side adds a matching stub here in the same commit; type checking is
not a follow-up task. The naming convention mirrors the runtime
exactly — no aliases, no rebadging — so a researcher reading the IDE
hover and running the same name at the REPL sees identical
signatures.
"""

__version__: str
