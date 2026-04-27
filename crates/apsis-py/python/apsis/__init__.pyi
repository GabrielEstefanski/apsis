"""Type stubs for the public ``apsis`` package surface.

The class-level signatures are owned by ``apsis/_native.pyi``, which
mirrors the runtime API exposed by the Rust extension module
``apsis._native``. This file simply re-exports the public symbols to
match what ``apsis/__init__.py`` does at runtime; it is the single
type-checker entry point a researcher's IDE consults when they write
``import apsis``.

Each subsequent PR that adds a class or free function to the Rust
side adds a matching declaration in ``_native.pyi`` and a re-export
line here in the same commit; type checking is not a follow-up task.
"""

from apsis._native import (
    Body as Body,
    IntegratorKind as IntegratorKind,
    System as System,
    __version__ as __version__,
)

__all__ = [
    "Body",
    "IntegratorKind",
    "System",
    "__version__",
]
