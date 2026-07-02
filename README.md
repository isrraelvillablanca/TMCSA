# TMCSA
Teoría de la Membrana Constitutiva en Sistemas Abiertos — Scripts R canónicos de verificación empírica en 60 especies (5 reinos). Isrrael Villablanca Fuentes, Junio 2026.
# TMCSA — Theory of the Constitutive Membrane in Open Systems

**Author:** Isrrael Villablanca Fuentes  
**Affiliation:** Independent researcher, Región de Ñuble, Chile  
**Year:** 2026  
**License:** Licencia dual — código bajo MIT (LICENSE-CODE.txt), documentos bajo CC-BY-4.0 (LICENSE-DOCS.txt)

---

## What is TMCSA?

The Theory of the Constitutive Membrane in Open Systems (TMCSA, from Spanish *Teoría de la Membrana Constitutiva en Sistemas Abiertos*) is a theoretical biology framework that models aging, senescence, and the cessation of life using a single differential equation [EM-2] whose **form does not change across species** — only the species-specific parameter set R(species) varies.

The master equation [EM-2] has three regimes:

- **Regime A (saturating repair):** standard senescence — damage accumulates until the vitality threshold Dc=0.60 is crossed
- **Regime B (proportional repair):** negligible senescence — damage stabilizes at X_eq = η_basal/ρ_rep, far below Dc
- **Regime C (conditional reversion):** rejuvenation — conditional on Regime A, triggered by specific stressors (documented in *Turritopsis dohrnii*)

The framework has been verified empirically across **60 species from 5 kingdoms** (Mammalia, Aves, Reptilia, Echinodermata, Fungi, Plantae, Bacteria), spanning 5 orders of magnitude in longevity — from 6-day cereal crops to 10,000-year aspen clones — using the same equation form throughout.

---

## Repository Contents

### Canonical scripts (10A–10E)
Full verification across all 60 species, organized by biological group:

| Script | Version | Species (n) | Regime | Status |
|--------|---------|-------------|--------|--------|
| `TMCSA_canonico_10A_EM2_v7.R` | v7 | 11 | A | ✓ Verified in R — deterministic branch corrected |
| `TMCSA_canonico_10B_EM2_v8.R` | v8 | 14 (12A + 2B) | A + B | ✓ Verified in R |
| `TMCSA_canonico_10C_EM2_v8.R` | v8 | 13 (plants) | A | ✓ Verified in R |
| `TMCSA_canonico_10D_EM2_v2.R` | v2 | 10 (diverse animals) | A | ✓ Verified in R |
| `TMCSA_canonico_10E_EM2_v5.R` | v5 | 10 (extreme plants) | A | ✓ Verified in R |

### Individual Regime B scripts
Calibrated individually before integration into the canonical group:

| Script | Species | Error | Notes |
|--------|---------|-------|-------|
| `TMCSA_RegB_Spurpuratus_v3.R` | *Strongylocentrotus purpuratus* | 1.6% | MS=180 |
| `TMCSA_RegB_Sfranciscanus_v3.R` | *S./M. franciscanus* | 0.5% | MS=360 |
| `TMCSA_EM2_Hydra_canonico_v2.R` | *Hydra vulgaris* | 4.7% | ρ_rep LEVEL 1/2 (David & Campbell 1972, T_cycle=3.5d → 72.29/yr); η_basal pending direct oxidative damage measurement |

### Individual species scripts
| Script | Species | Error | Notes |
|--------|---------|-------|-------|
| `TMCSA_Drimys_canonico_v2.R` | *Drimys winteri* (canelo, Chilean native tree) | 3.4% | W_met=0.75 (seasonal latency); includes real growth data panels (Navarro 1993) |

### Graphics script
| Script | Version | Generates |
|--------|---------|-----------|
| `TMCSA_AnexoA_graficas_v10.R` | v10 | 5 PNG figures: S(t), t_cese distribution, Gompertz h(t) log scale, V(t)/D(t) Regime A, V(t)/D(t) Regime B, Safety Margin P11, global parity |

### Analysis documents
| Document | Version | Contents |
|----------|---------|----------|
| `TMCSA_identificabilidad_paso1a4_v4.docx` | v4 | Structural and practical identifiability, subset selection, local sensitivity analysis for [EM-2]. Central finding: t_cese is blind to {a,b,μ,n} across all three main species families — T_activo must be included as standard observable. Recommended by Consensus (July 2026). |

---

## How to Run

All scripts require **R** (version ≥ 4.0). No additional packages beyond base R are needed for the canonical scripts. The graphics script requires base R graphics only.

```r
# Parse before running (mandatory — rule B.12)
Rscript -e 'invisible(parse(file="TMCSA_canonico_10A_EM2_v7.R"))'

# Run
Rscript TMCSA_canonico_10A_EM2_v7.R
```

Each canonical script produces:
- Console output with median t_cese and error for each species
- PNG figures (saved to working directory)

**Reproducibility note:** All scripts use fixed seeds (`set.seed(42)` by default). Results are reported as averages of 5 seeds for species with external cessation mode (long-tail distributions), following the methodological rule documented in Manual Operativo [EM-2] §6bis.4. This rule was formalized after a real case where a single seed produced 7.4% error vs. 1.6% with 5-seed average — documented in script headers.

---

## Key Results

| Indicator | Value |
|-----------|-------|
| Master equation | [EM-2] — 3 regimes, 15 parameters in R(species) |
| Species verified | 60 (5 kingdoms, 5 orders of magnitude in longevity) |
| Median error across species | <5% (most); <10% all except *Dracaena* ~12% (documented long-tail variance) |
| Regime B (quantitative) | *Myotis brandtii* 1.0%, *S. purpuratus* 1.6%, *S. franciscanus* 0.5% |
| Universal candidate | Dc=0.60 — damage threshold, consistent across 60 species |
| Independent convergences | Karin & Raz et al. 2019 (Nat. Comm.); Kogan et al. 2015 (Sci. Rep.) |

---

## Algebraic Properties

The framework derives 15 algebraic properties (P1–P15) from [EM-2]:

- **P1–P10:** Regime A (saturating repair) — fixed point stability, Gompertz emergence, maturity law
- **P11–P14:** Regime B (proportional repair) — bifurcation threshold ρ_rep_critical = η_basal/(Xc·Dc), Safety Margin scaling with longevity
- **P15:** Regime C (conditional reversion) — ρ depends on stressor type, not only on species identity

**Key identifiability finding (July 2026):** t_cese is structurally blind to parameters {a, b, μ, n} across all main species families. This is not a hidden limitation — it follows algebraically from property P2 (stable fixed point V* = (a/b)^(1/(μ-1))). T_activo (integral of V(t)) must be included as a standard observable for calibration of these four parameters. See `TMCSA_identificabilidad_paso1a4_v4.docx`.

---

## Independent Convergences

Two research groups arrived independently at the same bifurcation structure as [EM-2]:

1. **Karin, Alon & Raz et al. (2019, Nature Communications)** — modeled senescent cell accumulation and found that sufficiently high repair rates stabilize the system, preventing senescence. Same bifurcation as P11.

2. **Kogan, Molodtsov, Menshikov, Shmookler Reis & Fedichev (2015, Scientific Reports)** — analyzed gene network stability motivated explicitly by naked mole-rat and long-lived sea urchins (the same Regime B reference species in TMCSA). Found the same threshold structure.

Neither group knew of TMCSA. The convergence is structural, not coincidental.

---

## Citation

If you use these scripts in your research, please cite:

> Villablanca Fuentes, I. (2026). *Teoría de la Membrana Constitutiva en Sistemas Abiertos (TMCSA)*. Independent researcher, Región de Ñuble, Chile. GitHub: https://github.com/isrraelvillablanca/TMCSA

---

## Contact

Isrrael Villablanca Fuentes  
Independent researcher — Región de Ñuble, Chile
