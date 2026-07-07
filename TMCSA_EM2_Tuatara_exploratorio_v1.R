# ==============================================================================
# TMCSA [EM-2] -- Sphenodon punctatus (tuatara) -- Regimen B -- CALIBRACION
# EXPLORATORIA v1 -- aporte CUALITATIVO, NO VERIFICADO (igual que Hydra antes
# del dato de David & Campbell)
# ==============================================================================
# Molde base: Myotis brandtii (parametros reales de TMCSA_AnexoA_graficas_v10.R,
# funcion simular_B): eta_basal=0.06, rho_rep=10.0, Xc=3.0, a/b/mu/n universales,
# eps=0.180*(Xc/17).
#
# ORIGEN DE CADA PARAMETRO:
#
#   Gamma_ext = ln(2)/91 = 0.00762/año -- NIVEL 1.
#     Derivado directamente del target empirico (ver abajo). No requiere dato
#     externo adicional -- ya es consecuencia del target real.
#
#   t_cese_target = 91 años -- NIVEL 1.
#     Fuente: N.J. Nelson (Victoria University of Wellington), comunicacion
#     personal citada en Moore et al., Biological Conservation -- longevidad
#     maxima REGISTRADA EN TERRENO (no la media modelada de 137a de Reinke
#     et al. 2022, que es una salida estadistica, no una edad observada).
#
#   W_met = 0.55 -- NIVEL 1.
#     Fuente: Wilson & Lee (1970), Comparative Biochemistry and Physiology --
#     consumo de O2 en reposo = 55% del valor predicho para un lagarto de
#     igual masa corporal (5-36 grados C).
#
#   eta_basal = 0.0044/año -- NIVEL 3, CANDIDATO. Valor medio del rango
#     [0.0024, 0.0063] derivado por analogia con Myotis brandtii (Metodo de
#     Seis Pasos, Paso 5 -- Operador R: temperatura corporal operativa).
#     eta_basal(tuatara) = eta_basal(Myotis) x W_met x Q10^(-DeltaT/10)
#     DeltaT = 37-13 = 24 grados C (mamifero referencia vs. tuatara operativo).
#     Q10 no medido para esta especie -- se uso el rango estandar de la
#     literatura bioquimica general (2.0-3.0). Valor medio usado aqui: 0.0044.
#     PENDIENTE: dato directo de tasa de daño oxidativo o error de replicacion
#     en Sphenodon punctatus. Mismo estatus que Hydra y H. glaber.
#
#   rho_rep = 0.72/año -- NIVEL 3, CANDIDATO. Mismo tratamiento que eta_basal
#     (mismo factor combinado W_met x Q10), preservando el Margen de Seguridad
#     de Myotis (MS=300) por construccion -- SUPUESTO DECLARADO: la reparacion
#     proporcional se asume igual de sensible a la temperatura que la
#     produccion de daño. Si la cinetica de reparacion es menos sensible al
#     frio (plausible en ectotermos con enzimas adaptadas), rho_rep real
#     podria ser mayor y el MS real mayor a 300.
#     PENDIENTE: dato directo de tasa de renovacion celular o aclaramiento en
#     Sphenodon punctatus.
#
#   Xc=3.0, a=0.0988, b=0.0960, mu=3.0, n=2.0, eps=0.180*(3.0/17) -- heredados
#     del molde Myotis brandtii, NIVEL 3 (no medidos para esta especie).
#
#   CV_eta = 0.30 -- NIVEL 3, SUPUESTO. El script/funcion de Myotis brandtii
#     disponible (funcion ilustrativa simular_B en TMCSA_AnexoA_graficas_v10.R)
#     no expone este valor -- usa un eta_i fijo=0.06 para graficar, no un
#     sorteo log-normal. Se tomo prestado el CV_eta=0.30 de S. purpuratus (la
#     especie de regimen B con script de calibracion completo mas cercana
#     disponible) como aproximacion. PENDIENTE: verificar si existe el script
#     de calibracion completo de Myotis brandtii con su propio CV_eta.
#
# VERIFICACION PYTHON PREVIA (Paso 1 formal, regla B.12), exploratoria:
#   X_eq = eta_basal/rho_rep = 0.0044/0.72 = 0.006111 (identico en ambos
#   extremos del rango de Q10, por construccion). MS = Dc*Xc/X_eq = 294.5
#   (practicamente identico al MS=300 de Myotis, como se buscaba).
#
#   PENDIENTE: este script no ha sido parseado en R real. Correr antes de
#   ejecutar:
#     Rscript -e 'invisible(parse(file="TMCSA_EM2_Tuatara_exploratorio_v1.R"))'
#
# ==============================================================================

Dc <- 0.60

# ------------------------------------------------------------------------------
# R(Sphenodon punctatus) -- Regimen B, exploratorio v1
# ------------------------------------------------------------------------------
R_Tuatara_B <- list(
  nombre        = "Sphenodon punctatus",
  regimen       = "B",
  t_cese_target = 91.0,
  a=0.0988, b=0.0960, mu=3.0, n=2.0,
  Xc=3.0,
  eta_basal = 0.0044,              # NIVEL 3 -- candidato, valor medio del rango [0.0024,0.0063]
  rho_rep   = 0.72,                # NIVEL 3 -- candidato, valor medio del rango [0.39,1.04]
  Gamma_ext = log(2)/91.0,         # NIVEL 1 -- derivado del target empirico (Nelson)
  eps = 0.180*(3.0/17), W_met = 0.55,
  CV_eta = 0.30                    # NIVEL 3 -- supuesto, prestado de S. purpuratus
)

# ------------------------------------------------------------------------------
# FUNCION DE SIMULACION -- identica en estructura al motor de Hydra v2
# (TMCSA_EM2_Hydra_canonico_v2.R), sin modificar la logica del motor.
# ------------------------------------------------------------------------------
simular_tuatara_regB <- function(R, N_sto=300, seed=42, T_max_det=NULL) {

  dt        <- 1/365
  dt_sem    <- 7/365
  eps_sem   <- R$eps * sqrt(7)
  gamma_ext <- R$Gamma_ext

  CV_eta   <- R$CV_eta
  sigma_ln <- sqrt(log(1 + CV_eta^2))
  mu_ln    <- log(R$eta_basal) - sigma_ln^2/2

  set.seed(seed)

  if(is.null(T_max_det)) T_max_det <- R$t_cese_target * 2.5
  N_t <- round(T_max_det / dt)
  V_d <- rep(1.0, N_t); D_d <- rep(0.0, N_t)
  t_d <- (0:(N_t-1)) * dt

  X_d <- 0.0
  for(i in 2:N_t){
    ta <- t_d[i-1]; Va <- V_d[i-1]; Da <- D_d[i-1]
    if(isTRUE(Va <= 0.001)){ V_d[i]=0; D_d[i]=Da; next }
    dX_d <- (R$eta_basal - R$rho_rep * X_d) * dt
    X_d  <- max(0, X_d + dX_d)
    Da_nueva <- X_d / R$Xc
    fric <- max(0, (1 - Da/Dc)^R$n)
    dV   <- R$a * Va * fric - R$b * Va^R$mu - max(0, Da - Dc)
    V_d[i] <- min(1.10, max(0, Va + dV*dt))
    D_d[i] <- Da_nueva
  }

  idx_vida <- which(t_d <= R$t_cese_target)
  D_max_interno <- max(D_d[idx_vida], na.rm=TRUE)
  X_eq <- R$eta_basal / R$rho_rep
  MS   <- (Dc * R$Xc) / X_eq

  T_sto <- R$t_cese_target * 4
  N_s   <- round(T_sto / dt_sem)
  t_cese_sto <- rep(NA, N_sto)
  causa_ext  <- rep(FALSE, N_sto)

  for(sim in 1:N_sto){
    eta_i <- rlnorm(1, mu_ln, sigma_ln)
    X <- 0.0; V_i <- 1.0; t_fin <- NA; fue_ext <- FALSE
    for(sem in 1:N_s){
      ta <- (sem - 1) * dt_sem
      xi <- rnorm(1, 0, 1)

      if(runif(1) < gamma_ext * dt_sem){
        t_fin <- ta; fue_ext <- TRUE; break
      }

      dX <- (eta_i - R$rho_rep * X) * dt_sem + eps_sem * xi
      X  <- max(0, X + dX)
      D_i <- X / R$Xc

      fric <- max(0, (1 - D_i/Dc)^R$n)
      dV  <- R$a * V_i * fric - R$b * V_i^R$mu - max(0, D_i - Dc)
      V_i <- min(1.10, max(0, V_i + dV*dt))

      if(isTRUE(X >= R$Xc) || isTRUE(V_i <= 0.01)){ t_fin <- ta; break }
    }
    t_cese_sto[sim] <- t_fin; causa_ext[sim] <- fue_ext
  }

  tv <- t_cese_sto[!is.na(t_cese_sto)]
  pct_ext <- mean(causa_ext[!is.na(t_cese_sto)]) * 100

  cat(sprintf("  %s [Regimen B, exploratorio v1]: mediana=%.1fa  target=%.1fa  error=%.1f%%\n",
      R$nombre, median(tv), R$t_cese_target,
      abs(median(tv)-R$t_cese_target)/R$t_cese_target*100))
  cat(sprintf("    eta_basal=%.4f, rho_rep=%.2f (ambos NIVEL 3, candidatos)\n",
      R$eta_basal, R$rho_rep))
  cat(sprintf("    X_eq = eta_basal/rho_rep = %.6f  <<  Dc*Xc = %.3f  (MS = %.0f)\n",
      X_eq, Dc*R$Xc, MS))
  cat(sprintf("    D_max interno (sin cese ext) = %.6f -- %s\n",
      D_max_interno,
      ifelse(D_max_interno < Dc,
             "NO cruza Dc por si solo (cese es externo, correcto)",
             "cruza Dc -- revisar eta_basal/rho_rep")))
  cat(sprintf("    cese por causa externa (Gamma_ext) = %.1f%% de los individuos\n",
      pct_ext))
  cat("    [ADVERTENCIA] eta_basal y rho_rep son NIVEL 3 (candidatos, por analogia).\n")
  cat("    Requiere dato directo de tasa de daño oxidativo/renovacion celular.\n")
  cat("    CV_eta=0.30 es un supuesto prestado de S. purpuratus, no medido.\n")

  list(R=R, t_d=t_d, V_d=V_d, D_d=D_d,
       t_cese_sto=t_cese_sto, mediana=median(tv),
       D_max_interno=D_max_interno, pct_ext=pct_ext,
       X_eq=X_eq, MS=MS)
}

# ------------------------------------------------------------------------------
# CORRER CON PROMEDIO DE 5 SEMILLAS (regla B.12 / SS6bis.4), N=300 igual que Hydra
# ------------------------------------------------------------------------------
cat("=== TMCSA [EM-2] -- Sphenodon punctatus -- Regimen B exploratorio v1 ===\n\n")

semillas <- c(1, 2, 3, 4, 5)
medianas <- numeric(length(semillas))
pcts_ext <- numeric(length(semillas))

for(i in seq_along(semillas)){
  res <- simular_tuatara_regB(R_Tuatara_B, N_sto=300, seed=semillas[i])
  medianas[i] <- res$mediana
  pcts_ext[i] <- res$pct_ext
}

cat(sprintf("\n--- PROMEDIO 5 SEMILLAS ---\n"))
cat(sprintf("Mediana promedio: %.1fa (sd=%.1f)\n", mean(medianas), sd(medianas)))
cat(sprintf("Target (Nelson, longevidad maxima registrada): %.1fa\n", R_Tuatara_B$t_cese_target))
cat(sprintf("Error promedio: %.1f%%\n",
    abs(mean(medianas)-R_Tuatara_B$t_cese_target)/R_Tuatara_B$t_cese_target*100))
cat(sprintf("Cese externo promedio: %.1f%%\n", mean(pcts_ext)))
cat("\nRECORDATORIO: eta_basal y rho_rep son candidatos NIVEL 3 (por analogia con\n")
cat("Myotis brandtii, corregidos por W_met y Q10 termico no medido). Este resultado\n")
cat("es un aporte CUALITATIVO -- confirma si el diseño del regimen B es plausible\n")
cat("para esta especie, NO una calibracion verificada. Igual estatus que Hydra\n")
cat("antes del dato de David & Campbell.\n")
