"""Smoke tests for the binding's import surface.

These tests are deliberately minimal — they exercise *only* the
boundary between Python and the Rust extension module, not the
behaviour of the underlying integrators or force models. Behavioural
correctness is the responsibility of the parent crate's test suite
(``crates/apsis/tests``, ``crates/apsis-1pn/tests``) and of the
cross-implementation parity portfolio under ``validation/``.

Each new ``#[pyclass]`` or ``#[pyfunction]`` added on the Rust side
acquires a matching one-liner here that asserts the symbol is
importable through ``apsis`` and has the expected runtime type.
That is the full scope of this file: prove the façade is wired up.
"""

from __future__ import annotations

import re

import apsis


def test_module_imports() -> None:
    """``import apsis`` succeeds once the extension module is built."""
    assert apsis is not None


def test_version_string_matches_semver() -> None:
    """``apsis.__version__`` is sourced from the workspace ``Cargo.toml``.

    The test verifies the format rather than a specific value so the
    binding does not require an update on every workspace version
    bump; the upstream ``[workspace.package].version`` field remains
    the single source of truth.
    """
    version = apsis.__version__
    assert isinstance(version, str)
    assert re.fullmatch(r"\d+\.\d+\.\d+(?:[-+].+)?", version), (
        f"unexpected version string: {version!r}"
    )
