# ==============================================================================
# TMCSA [EM-2] -- Hydra vulgaris -- Regimen B (v2: rho_rep actualizado a dato
# directo de ciclo celular epitelial)
# ==============================================================================
# CAMBIO respecto a v1: rho_rep pasa de NIVEL 2 (turnover corporal completo,
# ~24 dias, derivado por analogia de tiempo de residencia) a NIVEL 1/2 (dato
# cinetico directo de division celular, David & Campbell 1972).
#
# ORIGEN DE CADA PARAMETRO:
#
#   Gamma_ext = 0.00497/año -- NIVEL 1 (sin cambio respecto a v1).
#     Fuente: Schaible et al. 2015, PNAS 112(51):15701-15706.
#
#   rho_rep = ln(2)/3.5 * 365 = 72.29/año -- NIVEL 1/2, ACTUALIZADO.
#     Fuente: David CN & Campbell RD 1972 -- ciclo celular de las celulas
#     epiteliales de Hydra = 3-4 dias (medido por marcaje [3H]-timidina).
#     Sintetizado en Bosch TCG & David CN 1984, 1987, 1991; Bosch et al. 2010
#     ("The Hydra polyp: nothing but an active stem cell community").
#     Justificacion (Metodo de Seis Pasos, Paso 5 -- Operador R): el tejido
#     epitelial ES la masa estructural dominante de Hydra -- todas las celulas
#     epiteliales de la region gastrica son celulas madre en division activa
#     (Bosch et al. 2010; Wittlieb et al. 2006), no una poblacion nicho minoritaria
#     como en vertebrados. Condicion que anularia este termino: bloqueo de la
#     division epitelial (ej. hidroxiurea) colapsaria la renovacion tisular.
#     Formula: cada division diluye el daño acumulado por celula a la mitad
#     -- tasa de dilucion poblacional = ln(2)/T_ciclo (vida media), no el
#     reciproco simple 1/T_ciclo. T_ciclo=3.5d es el punto medio del rango 3-4d.
#     NO se usa el ciclo intersticial (16-30h, mas rapido) ni el glandular
#     (64-83h) como termino principal -- ver nota metodologica abajo.
#
#   eta_basal = 0.05/año -- NIVEL 3, CANDIDATO, SIN CAMBIO. PENDIENTE.
#     La tasa de division celular (arriba) informa la DILUCION del daño, no
#     su PRODUCCION. No existe en la literatura revisada (Bosch, David,
#     Campbell, Schmidt) una medicion directa de tasa de daño oxidativo o de
#     errores de replicacion en Hydra. No se estima este valor a partir de
#     datos de proliferacion -- serian magnitudes biologicas distintas.
#     Requiere: medicion directa de tasa de daño oxidativo (ej. produccion de
#     ROS por unidad de tiempo) o de tasa de error de replicacion en Hydra.
#     Mismo estatus que H. glaber (Registro de Especies, ficha A.EM1-6).
#
#   Xc=4.0, a=0.0988, b=0.0960, mu=1.0, n=1.5, W_met=1.0, eps=0.01,
#   CV_eta=0.10 -- SIN CAMBIO respecto a v1.
#
#   t_cese_target = ln(2)/Gamma_ext = 139.5 años -- SIN CAMBIO (depende solo
#     de Gamma_ext, no de rho_rep).
#
# NOTA METODOLOGICA -- por que no se promedian los tres tipos celulares:
#   Interstitial (16-30h, Campbell & David 1974): linaje minoritario
#     (neuronas, nematocitos, gametos), no estructura corporal dominante.
#   Glandular (64-83h, Schmidt & David 1986): linaje minoritario, funcion
#     secretora, no estructural.
#   Ninguno de los dos se usa como termino principal porque el Paso 5 del
#   Metodo de Seis Pasos identifica el tejido que domina el daño acumulado
#   TOTAL del organismo -- ese es el epitelial, dado que constituye la masa
#   estructural del cuerpo. Quedan documentados aqui como referencia, no
#   como alternativa de calibracion.
#
# VERIFICACION PYTHON (Paso 1 formal, regla B.12), 5 semillas promediadas:
#   Fraccion viva simulada a los 139.5 años = 50.8% (teorico = 50.0%).
#   Cese por cruce de daño (X/Xc >= Dc) = 0.000% en 1500 simulaciones (300
#   individuos x 5 semillas), para eta_basal en {0.02, 0.05, 0.10, 0.20}.
#   X_eq = eta_basal/rho_rep cae entre 4.8x y 4.7x respecto a v1 (dado que
#   rho_rep subio de 15.21 a 72.29/año) -- MS (Margen de Seguridad) sube en
#   la misma proporcion: de 730 a 3470 con eta_basal=0.05.
#   PENDIENTE: este script no ha sido parseado -- Rscript no disponible en
#   este entorno. Correr antes de ejecutar:
#     Rscript -e 'parse(file="TMCSA_EM2_Hydra_canonico_v2.R")'
#
# ==============================================================================

Dc <- 0.60

# ------------------------------------------------------------------------------
# R(Hydra vulgaris) -- Regimen B, v2
# ------------------------------------------------------------------------------
R_Hydra_B <- list(
  nombre        = "Hydra vulgaris",
  regimen       = "B",
  t_cese_target = 139.5,
  a=0.0988, b=0.0960, mu=1.0, n=1.5,
  Xc=4.0,
  eta_basal = 0.05,               # NIVEL 3 -- candidato, pendiente dato directo
  rho_rep   = log(2)/3.5 * 365,   # NIVEL 1/2 -- 72.29/año, David & Campbell 1972
  Gamma_ext = 0.00497,            # NIVEL 1 -- Schaible et al. 2015
  eps=0.01, W_met=1.0,
  CV_eta=0.10
)

# ------------------------------------------------------------------------------
# FUNCION DE SIMULACION -- identica en estructura a v1 (motor corregido,
# produccion y reparacion escalan igual con dt_sem). Unico cambio: rho_rep.
# ------------------------------------------------------------------------------
simular_hydra_regB <- function(R, N_sto=300, seed=42, T_max_det=NULL) {

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

  cat(sprintf("  %s [Regimen B, v2]: mediana=%.1fa  target=%.1fa  error=%.1f%%\n",
      R$nombre, median(tv), R$t_cese_target,
      abs(median(tv)-R$t_cese_target)/R$t_cese_target*100))
  cat(sprintf("    rho_rep = %.2f/año (NIVEL 1/2, ln2/T_ciclo epitelial=3.5d)\n", R$rho_rep))
  cat(sprintf("    X_eq = eta_basal/rho_rep = %.6f  <<  Dc*Xc = %.3f  (MS = %.0f)\n",
      X_eq, Dc*R$Xc, MS))
  cat(sprintf("    D_max interno (sin cese ext) = %.6f -- %s\n",
      D_max_interno,
      ifelse(D_max_interno < Dc,
             "NO cruza Dc por si solo (cese es externo, correcto)",
             "cruza Dc -- revisar eta_basal/rho_rep")))
  cat(sprintf("    cese por causa externa (Gamma_ext) = %.1f%% de los individuos\n",
      pct_ext))
  cat(sprintf("    [ADVERTENCIA] eta_basal=%.2f sigue siendo NIVEL 3 (candidato).\n",
      R$eta_basal))
  cat("    Requiere dato directo de tasa de daño oxidativo o error de replicacion.\n")

  list(R=R, t_d=t_d, V_d=V_d, D_d=D_d,
       t_cese_sto=t_cese_sto, mediana=median(tv),
       D_max_interno=D_max_interno, pct_ext=pct_ext,
       X_eq=X_eq, MS=MS)
}

# ------------------------------------------------------------------------------
# CORRER CON PROMEDIO DE 5 SEMILLAS (regla B.12 / SS6bis.4)
# ------------------------------------------------------------------------------
cat("=== TMCSA [EM-2] -- Hydra vulgaris -- Regimen B v2 (rho_rep actualizado) ===\n\n")

semillas <- c(1, 2, 3, 4, 5)
medianas <- numeric(length(semillas))
pcts_ext <- numeric(length(semillas))

for(i in seq_along(semillas)){
  res <- simular_hydra_regB(R_Hydra_B, N_sto=300, seed=semillas[i])
  medianas[i] <- res$mediana
  pcts_ext[i] <- res$pct_ext
}

cat(sprintf("\n--- PROMEDIO 5 SEMILLAS ---\n"))
cat(sprintf("Mediana promedio: %.1fa (sd=%.1f)\n", mean(medianas), sd(medianas)))
cat(sprintf("Target (ln2/Gamma_ext): %.1fa\n", R_Hydra_B$t_cese_target))
cat(sprintf("Error promedio: %.1f%%\n",
    abs(mean(medianas)-R_Hydra_B$t_cese_target)/R_Hydra_B$t_cese_target*100))
cat(sprintf("Cese externo promedio: %.1f%%\n", mean(pcts_ext)))
cat("\nCOMPARAR este error contra v1 (4.7%, rho_rep=15.21) -- si baja, confirma\n")
cat("que el dato real de ciclo celular mejora la calibracion respecto a la\n")
cat("estimacion anterior por analogia de turnover corporal.\n")
