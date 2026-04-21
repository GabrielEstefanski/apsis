# ADR-001 — Substituir `steps_per_frame` por modelo de wall-time budget

**Data:** 2026-04-21
**Status:** Aceito — implementação pendente em `feat/wall-budget`

---

## Contexto

`steps_per_frame: u32` controla quantos passos de física são executados por frame. Funciona para integradores de custo fixo, mas quebra com IAS15 (adaptativo):

| Integrador | Force evals / step | ×100k steps/frame |
|---|---|---|
| VelocityVerlet | 1 | 100k evals |
| Yoshida-4 | 4 | 400k evals |
| IAS15 | ~14–23 | ~2.3M evals |

Sintomas: UI trava ao pausar/deletar com IAS15 + alto steps_per_frame; slider exige que o usuário conheça o custo interno do integrador.

**Quickwin aplicado em `develop`:** `MAX_BATCH_WALL_MS = 33` em `physics_thread.rs` — o batch quebra após 33ms de wall-clock mesmo que steps_per_frame não tenha sido atingido. Resolve o freeze imediato, não resolve a abstração.

---

## Decisão

Substituir `steps_per_frame` por **wall-time budget por batch** (`batch_budget_ms: u32`).

O physics thread roda steps até o orçamento de wall-clock ser consumido. O integrador determina quantos steps cabem — não o usuário.

```
// Antes
while i < steps_per_frame { system.step(); i += 1; }

// Depois
let deadline = Instant::now() + Duration::from_millis(batch_budget_ms);
while Instant::now() < deadline { system.step(); }
```

O display `yr/s` já existente passa a ser o feedback primário de velocidade.

---

## Alternativas descartadas

| Alternativa | Motivo |
|---|---|
| Manter quickwin como permanente | Não resolve abstração errada |
| Target sim rate (yr/s) | Mais complexo, requer estimativa de custo por integrador |

---

## Plano de implementação (`feat/wall-budget`)

**Passo 1 — `physics_thread.rs`**
- `PhysicsCmd::SetStepsPerFrame(u32)` → `SetBatchBudgetMs(u32)`
- `steps_per_frame: u32` → `batch_budget_ms: u32` no loop interno
- Inner loop: `while i < steps_per_frame` → `while Instant::now() < deadline`
- Manter `MAX_BATCH_WALL_MS` como hard cap de segurança acima do budget máximo

**Passo 2 — `PhysicsHandle`**
- `set_steps_per_frame` → `set_batch_budget_ms`
- Atualizar todos os call-sites em `ui.rs`

**Passo 3 — `SimulationApp`**
- `steps_per_frame: u32` → `batch_budget_ms: u32` (default sugerido: 8ms)

**Passo 4 — `playbar.rs`**
- Slider `×N steps` → `Xms` (range: 1–100ms)
- Remover hint "↑dt for speed" do IAS15 (torna-se desnecessário)

**Passo 5 — Snapshot / config**
- Verificar se `steps_per_frame` é persistido em `.grav` — se sim, campo ignorado na leitura (não quebra saves antigos)

**Passo 6 — Testes + smoke**
- `cargo test` — nenhum teste deve referenciar steps_per_frame diretamente
- VV: comportamento de throughput igual ou melhor
- IAS15: UI responsiva com qualquer budget

**Passo 7 — Merge e limpeza**
- Remover `MAX_BATCH_WALL_MS` do develop (redundante após merge) ou manter como double safety

---

## Consequências

**Positivo:** UI nunca congela; slider agnóstico ao integrador; `yr/s` é o feedback natural.

**Atenção:** `steps_per_frame` some da API — checklist de call-sites necessário antes de começar. Comportamento de throughput muda para VV em cenas leves (roda mais steps que antes num mesmo budget).

---

## Referências

- REBOUND: `reb_integrate(sim, tmax)` — integra até tempo simulado, não N steps
- `src/core/physics_thread.rs:620` — loop atual + `MAX_BATCH_WALL_MS`
- `src/app/panel/playbar.rs` — slider de velocidade atual
