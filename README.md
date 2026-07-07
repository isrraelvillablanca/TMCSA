# TMCSA — Theory of the Constitutive Membrane in Open Systems

**Author:** Isrrael Villablanca Fuentes  
**Affiliation:** Independent researcher, Región de Ñuble, Chile  
**Year:** 2026  
**License:** MIT (scripts) — see LICENSE file

---

## What is TMCSA?

The Theory of the Constitutive Membrane in Open Systems (TMCSA, from Spanish *Teoría de la Membrana Constitutiva en Sistemas Abiertos*) is a theoretical biology framework that models aging, senescence, and the cessation of life using a single differential equation [EM-2] whose **form does not change across species** — only the species-specific parameter set R(species) varies.

The master equation [EM-2] has three regimes:

- **Regime A (saturating repair):** standard senescence — damage accumulates until the vitality threshold Dc=0.60 is crossed
- **Regime B (proportional repair):** negligible senescence — damage stabilizes at X_eq = η_basal/ρ_rep, far below Dc
- **Regime C (conditional reversion):** rejuvenation — conditional on Regime A, triggered by specific stressors (documented in *Turritopsis dohrnii*)

The framework has been verified empirically across **60 species from 5 kingdoms** (Mammalia, Aves, Reptilia, Echinodermata, Fungi, Plantae, Bacteria), spanning 5 orders of magnitude in longevity — from 6-day cereal crops to 10,000-year aspen clones — using the same equation form throughout.

*Note: beyond the 60 verified species, the repository also includes exploratory Regime B contributions (Level 3, "partial support") for species with no direct cellular data yet — see the Individual Regime B scripts table below. These do not count toward the 60.*

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
| `TMCSA_EM2_Tuatara_exploratorio_v1.R` | *Sphenodon punctatus* (tuatara) | 9.5% | **EXPLORATORY — Level 3, "partial support" (Contribution 30).** η_basal, ρ_rep derived by analogy with *Myotis brandtii*, corrected by W_met=0.55 (measured, Wilson & Lee 1970). Not a verified calibration. |
| `TMCSA_EM2_Balaena_mysticetus_RegB_v3.R` | *Balaena mysticetus* (bowhead whale) | 1.9% | **EXPLORATORY — Level 3, "partial support" (Contribution 31).** Target=211a, Level 2 (George et al. 1999, aspartic acid racemization, ±35a uncertainty). η_basal, ρ_rep by analogy with *Myotis brandtii*, W_met=0.35. |
| `TMCSA_EM2_Aldabrachelys_RegB_v1.R` | *Aldabrachelys gigantea* (Aldabra giant tortoise) | 4.3% | **EXPLORATORY — Level 3, "partial support" (Contribution 32).** Target=192a, Level 2/3 ("Jonathan", documentary/historical verification, weaker than the whale's biochemical method). W_met≈0.85 — unlike the tuatara, no evidence of hypometabolism relative to other ectotherms; longevity mechanism not directly identified. |

### Individual species scripts
| Script | Species | Error | Notes |
|--------|---------|-------|-------|
| `TMCSA_EM2_Drimys_canonico_v2.R` | *Drimys winteri* (canelo, Chilean native tree) | 3.4% | W_met=0.75 (seasonal latency); includes real growth data panels (Navarro 1993) |

### Sex differentiation script
| Script | Species | Result | Notes |
|--------|---------|--------|-------|
| `TMCSA_EM2_Hsapiens_Sexo_v5.R` | *Homo sapiens* (male vs female) | F−M diff = 5.0a, error 0.5% vs CDC 2023 target ~5a | Correct channel: δη_H (males produce SC faster, not δβ). N=2000 per sex. Consistent with [EF-R6] Family A (γ_H/γ_M=1.039). Universal structure: η_sex = η_species × (1 + δη_sex) |

### Global sensitivity analysis (Sobol)
Global sensitivity analysis (Sobol indices, `soboljansen`) across four species families, testing whether parametric identifiability depends on the dominant cessation mode rather than the repair regime:

| Script | Family | Species | Key finding |
|--------|--------|---------|-------------|
| `TMCSA_SensibilidadGlobal_Sobol_Hsapiens_v4.R` | 1 — internal cessation (Regime A) | *Homo sapiens* | η, β dominant directly (S1=0.26–0.33); a,b,μ,n act only through interaction (ST≈0.18–0.19) |
| `TMCSA_SensibilidadGlobal_Sobol_Hglaber_v2.R` | 2 — external cessation (Regime A, Γ_ext dominant) | *Heterocephalus glaber* | Structural non-identifiability — no parameter resolvable, direct or interactive, with any observable tested. Confirmed not a sampling-size artifact (N_ind=20→50 did not resolve it) |
| `TMCSA_SensibilidadGlobal_Sobol_Spurpuratus_v3.R` | 3 — mixed (Regime B, Γ_ext perturbed) | *Strongylocentrotus purpuratus* | Γ_ext dominant (S1=0.62); internal parameters retain real structured total effect (ST≈0.28–0.34), unlike Family 2. Conditional-identifiability post-hoc (X_eq by Γ_ext tercile) tested via five independent methods — all negative, closed as a definitive null result. |
| `TMCSA_SensibilidadGlobal_Sobol_Turritopsis_v1.R` | C — discrete reset mechanism | *Turritopsis dohrnii* | umbral_stress shows the single most concentrated direct effect across all four families (S1=0.833) — converges independently with prior profile-likelihood work on the same parameter pair |

**Unifying finding:** parametric identifiability in [EM-2] depends on which cessation mode dominates (internal vs. external), not on the repair regime (A vs. B vs. C). See `TMCSA_SensibilidadGlobal_CuatroFamilias_Consolidado_v3.docx`.

### Graphics script
| Script | Version | Generates |
|--------|---------|-----------|
| `TMCSA_AnexoA_graficas_v10.R` | v10 | 5 PNG figures: S(t), t_cese distribution, Gompertz h(t) log scale, V(t)/D(t) Regime A, V(t)/D(t) Regime B, Safety Margin P11, global parity |

### Analysis documents
| Document | Version | Contents |
|----------|---------|----------|
| `TMCSA_identificabilidad_paso1a4_v4.docx` | v4 | Structural and practical identifiability, subset selection, local sensitivity analysis for [EM-2]. Central finding: t_cese is blind to {a,b,μ,n} across all three main species families — T_activo must be included as standard observable. Recommended by Consensus (July 2026). |
| `TMCSA_ValidacionPredictiva_Calment_Kimura_Documento_v3.docx` | v3 | Held-out predictive validation: η_F calibrated against Jeanne Calment's record, η_H derived exclusively from the independent EF-R6 biochemical channel, predicting Jiroemon Kimura's record with 0.65% error. Includes controversy declaration on the Calment record and pending candidates (gorilla, chimpanzee, bowhead whale). |
| `TMCSA_SensibilidadGlobal_CuatroFamilias_Consolidado_v3.docx` | v3 | Consolidated global sensitivity (Sobol) results across all four families, with the unifying identifiability hypothesis and its status per family. |
| `TMCSA_SensibilidadGlobal_Sobol_Familia3_Documento_v3.docx` | v3 | Full documentation of Family 3 (regime B, perturbed Γ_ext), including the five independent attempts to test conditional identifiability, all negative — a definitive closed result. |
| `TMCSA_RegimenB_NuevasEspecies_Documento_v6.docx` | v6 | Regime B corpus extension: candidate search and rejection process, and full documentation of Contributions 30–32 (tuatara, bowhead whale, Aldabra giant tortoise). |

---

## How to Run

All scripts require **R** (version ≥ 4.0). Scripts in the global sensitivity analysis section additionally require the `sensitivity` package (`install.packages("sensitivity")`). No other packages beyond base R are needed.

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

**Runtime note:** global sensitivity scripts (Sobol) are computationally heavier than calibration scripts — between 1.5 and 6 hours depending on family, at N_base=1000 (12,000 total model evaluations). Each prints elapsed time during the run and saves partial progress every 500 evaluations.

---

## Key Results

| Indicator | Value |
|-----------|-------|
| Master equation | [EM-2] — 3 regimes, 15 parameters in R(species) |
| Species verified | 60 (5 kingdoms, 5 orders of magnitude in longevity) |
| Median error across species | <5% (most); <10% all except *Dracaena* ~12% (documented long-tail variance) |
| Regime B (quantitative) | *Myotis brandtii* 1.0%, *S. purpuratus* 1.6%, *S. franciscanus* 0.5% |
| Regime B (exploratory, Level 3) | tuatara 9.5%, bowhead whale 1.9%, Aldabra giant tortoise 4.3% — not verified calibrations, see notes above |
| Universal candidate | Dc=0.60 — damage threshold, consistent across 60 species |
| Held-out predictive validation | 1 pair verified (Calment/Kimura), 0.65% error, no demographic data used to derive the predicted value |
| Global sensitivity (Sobol) | 4 families tested; identifiability shown to depend on cessation mode, not repair regime |
| Independent convergences | Karin & Raz et al. 2019 (Nat. Comm.); Kogan et al. 2015 (Sci. Rep.) |

---

## Algebraic Properties

The framework derives 15 algebraic properties (P1–P15) from [EM-2]:

- **P1–P10:** Regime A (saturating repair) — fixed point stability, Gompertz emergence, maturity law
- **P11–P14:** Regime B (proportional repair) — bifurcation threshold ρ_rep_critical = η_basal/(Xc·Dc), Safety Margin scaling with longevity
- **P15:** Regime C (conditional reversion) — ρ depends on stressor type, not only on species identity

**Key identifiability finding (July 2026):** t_cese is structurally blind to parameters {a, b, μ, n} across all main species families. This is not a hidden limitation — it follows algebraically from property P2 (stable fixed point V* = (a/b)^(1/(μ-1))). T_activo (integral of V(t)) must be included as a standard observable for calibration of these four parameters. See `TMCSA_identificabilidad_paso1a4_v4.docx`.

**Global sensitivity confirms and extends this (July 2026):** across four species families with distinct cessation mechanisms (internal collapse, external mortality, mixed, discrete reset), the identifiability of R(species) parameters depends on which cessation mode dominates, not on which of the three repair regimes governs internal dynamics. See `TMCSA_SensibilidadGlobal_CuatroFamilias_Consolidado_v3.docx`.

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
