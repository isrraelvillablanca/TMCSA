# ==============================================================================
# TMCSA -- SENSIBILIDAD GLOBAL (SOBOL) -- Familia 5 / regimen C (Turritopsis)
# Especie representativa: Turritopsis dohrnii (R_Turritopsis, script canonico
# 10A v7, mecanismo de reinicio discreto por transdegeneracion)
# Archivo: TMCSA_SensibilidadGlobal_Sobol_Turritopsis_v1.R
#
# DISEÑO: 8 parametros perturbados -- a, b, mu, n (parametros de la ecuacion
# de vitalidad, candidatos a universales, +-30%, igual que Familia 1), eta,
# p_stress, umbral_stress (propios de regimen C, +-30%), y Gamma_ext (+-50%,
# el parametro que determina el balance interno/externo, mismo criterio que
# Familia 3).
#
# DOS OBSERVABLES, CALCULADOS CON UN TRATAMIENTO DELIBERADAMENTE DISTINTO DE
# Gamma_ext (precision de Isrrael, jul 2026):
#
#   1. resets/individuo -- calculado en REGIMEN LABORATORIO, Gamma_ext FIJO
#      EN 0, sin importar el valor de Gamma_ext que traiga esa fila del
#      diseño. Es el mismo observable que ya uso exitosamente el profile
#      likelihood (Paso 2, TMCSA_Identificabilidad_Paso1a3_v3.docx, S7):
#      encontro un minimo claro para (p_stress, umbral_stress) alrededor de
#      los valores reales (umbral=0.10, p_stress=20, error=5.41%).
#
#   2. fraccion sin colapso interno a t_obs=4a -- calculado usando el valor
#      de Gamma_ext QUE SI trae esa fila del diseño (perturbado +-50%). Aqui
#      el balance interno/externo si importa: a Gamma_ext alto, mas
#      individuos mueren por depredacion antes de que el mecanismo de
#      reinicio o el colapso interno puedan decidir nada.
#
# RAZON DE ESTA SEPARACION: si el mismo Gamma_ext perturbado se usara para
# ambos observables, el Sobol no podria separar el efecto de Gamma_ext sobre
# resets (que DEBERIA salir cero, porque el observable de laboratorio nunca
# tiene depredacion) del efecto sobre la fraccion de supervivencia (que SI
# deberia ser real). Mezclarlos contaminaria la interpretacion de ambos.
#
# HIPOTESIS A TESTEAR: dado que el profile likelihood ya demostro que
# (p_stress, umbral_stress) son identificables via resets/individuo, se
# espera que el Sobol muestre ST real (no ruido) para ambos parametros en
# el observable de resets -- la primera vez que profile likelihood y Sobol
# apuntarian en la misma direccion para el mismo par.
#
# ADVERTENCIA DECLARADA -- eta esta calibrado circularmente en el script
# original (comentario del script: "calibrado para que rho sea necesario y
# suficiente", no una medicion independiente). Cualquier indice de
# sensibilidad de eta debe leerse con esa limitacion: no tiene el mismo
# respaldo biologico independiente que p_stress o umbral_stress.
#
# REGLA B.12 -- parsear antes de correr:
#   Rscript -e 'invisible(parse(file="TMCSA_SensibilidadGlobal_Sobol_Turritopsis_v1.R"))'
# ==============================================================================

if (!requireNamespace("sensitivity", quietly=TRUE)) install.packages("sensitivity")
library(sensitivity)

N_base <- 1000
N_ind  <- 20

set.seed(123)
Dc <- 0.60

# --- Punto calibrado real de Turritopsis dohrnii (R_Turritopsis, 10A v7) ----
BASE <- list(a=0.0988, b=0.0960, mu=1.0, n=1.5,
             eta=1.0e-02, p_stress=20.0, umbral_stress=0.10,
             Gamma_ext=2.0)
nombres <- names(BASE)
p <- length(nombres)

t_lat   <- 0.1      # NIVEL 3, sin fuente citada -- fijo, no perturbado (no forma parte de los 8)
Xc      <- 24.9932  # EQ-3, t_activo=1a -- fijo
beta    <- 54.75    # heredado -- fijo
eps     <- 0.26463  # heredado -- fijo
CV_eta  <- 0.1      # PASO5: clonal laboratorio + rejuvenece como Hydra -- fijo
t_obs   <- 4.0       # criterio de verificacion, laboratorio

# --- Rangos: +-30% para todos, EXCEPTO Gamma_ext (+-50%) --------------------
rango_bajo  <- sapply(BASE, function(v) v * 0.7)
rango_alto  <- sapply(BASE, function(v) v * 1.3)
rango_bajo["Gamma_ext"] <- BASE$Gamma_ext * 0.5
rango_alto["Gamma_ext"] <- BASE$Gamma_ext * 1.5

cat("Rangos usados por parametro:\n")
for (nm in nombres) cat(sprintf("  %-14s [%.6g, %.6g]\n", nm, rango_bajo[nm], rango_alto[nm]))

generar_muestra <- function(N) {
  m <- matrix(NA, nrow=N, ncol=p)
  colnames(m) <- nombres
  for (j in 1:p) m[,j] <- runif(N, rango_bajo[j], rango_alto[j])
  as.data.frame(m)
}

X1 <- generar_muestra(N_base)
X2 <- generar_muestra(N_base)

# ------------------------------------------------------------------------------
# MOTOR -- identico en estructura a simular_turritopsis (script canonico 10A
# v7), separado en dos funciones: una fuerza Gamma_ext=0 (resets), la otra
# usa el Gamma_ext perturbado de la fila (fraccion).
# ------------------------------------------------------------------------------
correr_turritopsis <- function(params, con_gamma_ext, N_ind) {
  a <- params[["a"]]; b <- params[["b"]]; mu <- params[["mu"]]; n <- params[["n"]]
  eta <- params[["eta"]]; p_stress <- params[["p_stress"]]
  umbral_stress <- params[["umbral_stress"]]
  gamma_ext <- if (con_gamma_ext) params[["Gamma_ext"]] else 0.0

  kappa_X  <- 0.5
  dt_sem   <- 7/365
  beta_sem <- beta * dt_sem
  eps_sem  <- eps * sqrt(7)
  sig      <- sqrt(log(1 + CV_eta^2))
  muln     <- log(eta) - sig^2/2
  T_sto    <- t_obs * 1.2
  N_s      <- round(T_sto / dt_sem)

  t_cese_sto <- rep(NA, N_ind); n_resets <- rep(0, N_ind)

  for (sim in 1:N_ind) {
    eta_i <- rlnorm(1, muln, sig)
    X <- 0.01; V_i <- 1.0; t_fin <- NA; nr <- 0
    for (sem in 1:N_s) {
      ta <- (sem - 1) * dt_sem

      if (con_gamma_ext && runif(1) < gamma_ext * dt_sem) { t_fin <- ta; break }

      D_i <- X / Xc
      if (D_i > umbral_stress && runif(1) < p_stress * dt_sem) {
        X <- 0.001; V_i <- 1.0; nr <- nr + 1; next
      }
      xi <- rnorm(1, 0, 1)
      dX <- eta_i * max(0, ta - t_lat) * 7 - beta_sem * X/(kappa_X + X) + eps_sem * xi
      X  <- max(0, X + dX)
      D_i <- X / Xc
      base_f <- 1 - min(D_i/Dc, 1.0)
      fric <- if (base_f > 0) base_f^n else 0
      dV <- a * V_i * fric - b * V_i^mu - max(0, D_i - Dc)
      V_i <- min(1.10, max(0, V_i + dV * dt_sem))
      if (X >= Xc || V_i <= 0.01) { t_fin <- ta; break }
    }
    t_cese_sto[sim] <- t_fin; n_resets[sim] <- nr
  }
  list(resets_prom = mean(n_resets),
       frac_sin_colapso = mean(t_cese_sto >= t_obs | is.na(t_cese_sto)))
}

simular_par <- function(params, N_ind=20, seed=1) {
  set.seed(seed)
  res_lab <- correr_turritopsis(params, con_gamma_ext=FALSE, N_ind=N_ind)  # Gamma_ext=0 FIJO
  set.seed(seed + 1000000)
  res_nat <- correr_turritopsis(params, con_gamma_ext=TRUE,  N_ind=N_ind)  # Gamma_ext perturbado
  list(resets = res_lab$resets_prom, frac_sin_colapso = res_nat$frac_sin_colapso)
}

# ------------------------------------------------------------------------------
# DISENO SOBOL -- dos objetos (resets, fraccion), mismo diseno verificado.
# ------------------------------------------------------------------------------
cat(sprintf("\nFamilia 5 (Turritopsis, regimen C) -- p=%d, N_base=%d, N_ind=%d -> evaluaciones = %d\n",
            p, N_base, N_ind, (p+2)*N_base))

sob_resets   <- soboljansen(model=NULL, X1=X1, X2=X2, nboot=200)
sob_fraccion <- soboljansen(model=NULL, X1=X1, X2=X2, nboot=200)
stopifnot(isTRUE(all.equal(as.matrix(sob_resets$X), as.matrix(sob_fraccion$X))))
cat("Verificado: los dos disenos son identicos fila por fila.\n")

X_diseno <- sob_resets$X
n_eval <- nrow(X_diseno)
Y_resets   <- numeric(n_eval)
Y_fraccion <- numeric(n_eval)

t_inicio <- Sys.time()
for (i in 1:n_eval) {
  res <- simular_par(X_diseno[i, ], N_ind=N_ind, seed=i)
  Y_resets[i]   <- res$resets
  Y_fraccion[i] <- res$frac_sin_colapso
  if (i %% 100 == 0 || i == n_eval) {
    transcurrido <- as.numeric(difftime(Sys.time(), t_inicio, units="secs"))
    cat(sprintf("  %d/%d  (%.1fs transcurridos, %.3fs/eval)\n",
                i, n_eval, transcurrido, transcurrido/i))
  }
  if (i %% 500 == 0) {
    saveRDS(list(X_diseno=X_diseno[1:i, , drop=FALSE],
                 Y_resets=Y_resets[1:i], Y_fraccion=Y_fraccion[1:i], ultimo_i=i),
            "TMCSA_Sobol_Turritopsis_progreso_completo.rds")
  }
}
cat(sprintf("\nTiempo total: %.1f min\n", as.numeric(difftime(Sys.time(), t_inicio, units="mins"))))

saveRDS(list(X_diseno=X_diseno, Y_resets=Y_resets, Y_fraccion=Y_fraccion),
        "TMCSA_Sobol_Turritopsis_diseno_completo_FINAL.rds")
cat("Diseno completo y resultados guardados en TMCSA_Sobol_Turritopsis_diseno_completo_FINAL.rds\n")

sob_resets   <- tell(sob_resets,   Y_resets)
sob_fraccion <- tell(sob_fraccion, Y_fraccion)

cat("\n", strrep("=", 78), "\nRESULTADOS -- SOBOL resets/individuo (Gamma_ext=0 FIJO)\n", strrep("=", 78), "\n", sep="")
print(sob_resets)
cat("\n", strrep("=", 78), "\nRESULTADOS -- SOBOL fraccion sin colapso interno (Gamma_ext perturbado)\n", strrep("=", 78), "\n", sep="")
print(sob_fraccion)

hacer_tabla <- function(sob, sufijo) {
  df <- data.frame(
    parametro = nombres,
    S1 = sob$S$original, S1_min = sob$S$`min. c.i.`, S1_max = sob$S$`max. c.i.`,
    ST = sob$T$original, ST_min = sob$T$`min. c.i.`, ST_max = sob$T$`max. c.i.`
  )
  setNames(df, c("parametro", paste0(c("S1","S1_min","S1_max","ST","ST_min","ST_max"), "_", sufijo)))
}
tabla_final <- merge(hacer_tabla(sob_resets,"resets"), hacer_tabla(sob_fraccion,"fraccion"), by="parametro")
write.csv(tabla_final, "TMCSA_SensibilidadGlobal_Turritopsis_resultados.csv", row.names=FALSE)
cat("\nGuardado en TMCSA_SensibilidadGlobal_Turritopsis_resultados.csv\n")

cat("\n[ADVERTENCIA] eta esta calibrado circularmente en el script original\n")
cat("(\"calibrado para que rho sea necesario y suficiente\", no medicion\n")
cat("independiente). Su indice de sensibilidad no tiene el mismo respaldo\n")
cat("biologico que p_stress o umbral_stress.\n")
