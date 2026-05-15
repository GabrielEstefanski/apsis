# apsis-python

PyO3 cdylib backing the `apsis` Python distribution. Built by maturin
via the root `pyproject.toml`; not built directly via `cargo build`.

```bash
# from repository root
maturin develop --release
```

## Features

| feature | submodule  | crate     |
|---------|------------|-----------|
| `gr`    | `apsis.gr` | apsis-1pn |
