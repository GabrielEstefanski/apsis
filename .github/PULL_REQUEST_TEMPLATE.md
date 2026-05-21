## Summary

What this PR changes and why. Link the issue this closes if applicable
(`Closes #N`).

## Test plan

- [ ] `cargo fmt --all` clean
- [ ] `cargo clippy --workspace --all-targets -- -D warnings` clean
- [ ] `cargo nextest run --workspace` passes
- [ ] (Python changes only) `maturin develop --release` succeeds and `pytest python-tests/` passes
- [ ] (Numerical claims) measurement reported with units; gate value justified

## References

- Linked issue(s):
- ADR (if architectural):
- Lab notebook (if experimental):
