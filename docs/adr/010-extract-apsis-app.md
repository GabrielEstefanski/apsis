# ADR-010 — Extract apsis-app to Separate Repository

**Status:** Accepted
**Date:** 2026-05-15
**Supersedes (in part):** ADR-005 §"Crate inventory"

---

## Context

`apsis-app` shipped as a workspace member from the project's first
commit: an egui/wgpu interactive shell that consumes the apsis
core to render and edit simulations live. By v0.1.0-alpha.1 the
shell carried:

- ~12 KLOC of UI / rendering code (panels, camera, GPU pipeline,
  trail rendering, design tokens)
- Heavy graphics dependency tree (`eframe`, `egui`, `wgpu`, `winit`,
  `pollster`, `egui-phosphor`, `egui_kittest`)
- Its own validation surface (kittest UI snapshots) that is unrelated
  to the library's physics validation portfolio

The library's published contribution is structural (federated
operator architecture). The shell is downstream — consumes the public
API, demonstrates it visually, but does not extend it. Keeping the
shell in the main workspace had three concrete costs:

1. **CI surface.** `cargo build --workspace` compiled the full graphics
   stack on every PR. mold + nextest + profile.test mitigated wall
   time but the dependency download / compile remained.
2. **Paper scope.** Reviewer reading the repo for the federation thesis
   had to disambiguate which crates are part of the validated surface
   and which are demonstration. README and docs/overview repeatedly
   noted "not part of the library's validated surface" to manage this.
3. **Release coupling.** Bumping the workspace version forced a coupled
   release of UI work-in-progress with library work, even when the
   library was paper-ready and the UI was mid-refactor.

A v0.2 release that ships the apsis library cleanly to JOSS / arxiv
is incompatible with carrying an in-flight UI shell on the same
release cadence.

## Decision

Move `apsis-app` to a standalone repository at
[`GabrielEstefanski/apsis-app`](https://github.com/GabrielEstefanski/apsis-app).
The new repository depends on the apsis core via Cargo's git
dependency (pinned to a tag) until apsis is published to crates.io,
at which point the dependency switches to a version pin.

The apsis core repository removes:

- `crates/apsis-app/` (entire directory)
- Workspace member + default-member entries for apsis-app in root
  `Cargo.toml`
- UI dependencies from `[workspace.dependencies]`: `eframe`, `egui`,
  `egui-phosphor`, `wgpu`, `winit`, `pollster`, `egui_kittest`
- README architecture-table row for apsis-app; replaced with a
  paragraph pointing at the new repository
- `docs/overview.md` architecture-table row and prose mention; same
  replacement pattern

The pre-extraction state of `apsis-app` is preserved on the apsis
core repository as the lightweight tag `apsis-app-pre-extraction`
for recoverability if the new repository is lost or corrupted.

## Consequences

**Architectural wins:**

- Apsis core repository is now exclusively the federated library
  surface (core + operators + Python distribution + capsule transport).
  "What is apsis" is unambiguous.
- CI shrinks materially. Graphics stack no longer compiles on apsis
  PRs.
- Federation thesis strengthens: even the visual frontend is now a
  downstream consumer of the public API, not a privileged in-tree
  shell. Demonstrates the contract works for downstream projects.
- Apsis library can release on its own cadence (paper-ready);
  apsis-app releases independently when the UI is ready.

**Migration cost:**

- Anyone with a script or tool that ran `cargo run -p apsis-app`
  needs to clone the new repository. No published wheel / binary
  consumers exist yet, so the migration cost is internal-only.
- ADR-001 and other historical docs reference paths under
  `crates/apsis-app/` that no longer exist in this repository. ADRs
  are immutable historical record (see `project_adr_discipline`); the
  references are not edited retroactively.
- The federation memory entry tracking app-side backlog
  (`project_app_federation_demo`, `project_post_fo_backlog`,
  `project_inspector_backlog`, `project_app_libs_backlog`) carries
  to the new repository's issue tracker on its own schedule.

**Apsis-app downstream contract:**

The new repository depends on apsis through the public extension
API only: `Body`, `System`, `IntegratorKind`, the perturbation
registration surface, the apsis-1pn / apsis-radiation / apsis-central
operator crates. It does not reach into `pub(crate)` types or
internal modules. If the apsis-app codebase relied on any
non-public surface at the time of extraction, that reliance is a
bug in apsis-app that surfaces at the next compatible apsis
release.

**Future re-merge anti-pattern:**

A future PR proposing "let's bring the GUI back in-tree to make
demos easier" should be rejected. The federation thesis depends on
the GUI being a downstream consumer; in-tree GUI re-couples the
release cadence and dilutes the library's published scope. If the
apsis core needs visualization for paper figures, the right path
is a paper-specific Python script under `paper/` that uses
matplotlib, not a re-imported egui shell.

## References

- New repository: <https://github.com/GabrielEstefanski/apsis-app>
- Preservation tag: `apsis-app-pre-extraction` on this repository
- ADR-005 (federated perturbation operators) — the federation
  thesis the extraction reinforces
- ADR-009 (consolidated Python distribution) — same pattern at the
  Python binding boundary
