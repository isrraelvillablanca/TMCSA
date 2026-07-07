# ==============================================================================
# TMCSA -- [EM-2] -- Balaena mysticetus (ballena de Groenlandia) -- Regimen B
# Archivo: TMCSA_EM2_Balaena_mysticetus_RegB_v3.R
# v3: añadida seccion de graficas (3 paneles: dist t_cese, V(t), D(t))
# Autor: Isrrael Villablanca Fuentes
# Fecha: julio 2026
#
# CORRECCION v1->v2 (jul 2026):
#   (1) Doble conteo de W_met: v1 multiplicaba dX <- eta_i * W_met_B - ...
#       dentro del motor, pero eta_basal_B ya habia sido derivado con la
#       correccion de W_met aplicada (eta_basal_B <- 0.06 * factor_Wmet,
#       factor_Wmet=W_met_B/0.6). Esto aplicaba la correccion metabolica DOS
#       veces, dejando el eta efectivamente simulado en eta_basal_B*W_met_B
#       (~0.01225) en vez de 0.035 como se reportaba. Corregido: dX <- eta_i
#       - rho_rep_B * X + eps_sem * rnorm(1), sin multiplicar por W_met_B de
#       nuevo. Mismo error corregido en el bucle de calculo de D_max.
#   (2) Regla de eps: v1 usaba eps_B <- 0.180*(Xc_B/4.0) = 0.135, sin
#       explicacion y sin coincidir con la regla usada en el resto del
#       corpus (Myotis, tuatara, H.glaber, S.purpuratus: eps=0.180*(Xc/17)).
#       Corregido a eps_B <- 0.180*(Xc_B/17) = 0.0318, consistente.
#
# APORTE 31 (Con soporte parcial) -- Regimen B exploratorio
#
# BASE EMPIRICA:
#   Target longevidad: 211 años (George JC et al. 1999, Can J Zool 77:571-580)
#     -- estimado mediante racemizacion de acido aspartico en cristalino,
#        individuo macho, el mas longevo de todos los mamiferos documentados.
#   Mecanismo molecular: SASP reducido en celulas senescentes vs humanos
#     (Gros et al. 2025, Nature, doi:10.1038/s41586-025-09694-5)
#     -- reparacion de ADN mejorada via CIRBP/RPA2 y mismatch repair
#     -- esto corresponde a eta_basal bajo y/o beta alto en [EM-2]
#   Madurez sexual: ~25 años (George et al. 1999)
#   Metabolismo: hipometabolico confirmado (George et al., adaptaciones
#     balaenidas) -- W_met estimado 0.35 (NIVEL 3)
#   Sin senescencia reproductiva reportada en machos hasta 159 años
#     (George et al. 1999) -- consistente con regimen B
#
# REGIMEN: B (proporcional) -- confirmado cualitativamente por:
#   - SASP reducido (mecanismo de reparacion activa documentado)
#   - Longevidad maxima >>200 años sin señales de colapso interno
#   - Analogia con otros mamiferos de regimen B (Myotis, H. glaber)
#
# PARAMETROS:
#   eta_basal: NIVEL 3 -- estimado desde Myotis brandtii con correccion
#     por W_met (0.35 vs 0.6 de Myotis) y masa corporal (~50000 kg)
#   rho_rep:   NIVEL 3 -- preserva MS de Myotis por construccion
#   W_met:     NIVEL 3 -- estimado por hipometabolismo documentado
#   t_lat:     NIVEL 2 -- madurez sexual ~25 años (George et al. 1999)
#   Gamma_ext: NIVEL 2 -- derivado del target de 211 años
#
# ADVERTENCIA: eta_basal, rho_rep y W_met son NIVEL 3 (candidatos).
#   El script confirma si el diseño de regimen B es plausible para esta
#   especie -- NO es una calibracion verificada.
#
# NOTA SOBRE TARGET: el record de 211 años es un macho. Para una especie
#   no calibrada previamente, se usa como target de t_cese el valor mas
#   documentado disponible. La incertidumbre del metodo AAR es ±35 años
#   (George et al. 1999) -- se documenta en los resultados.
#
# REGLA B.12 -- parsear antes de correr:
#   Rscript -e 'invisible(parse(file="TMCSA_EM2_Balaena_mysticetus_RegB_v3.R"))'
# ==============================================================================

set.seed(42)
Dc    <- 0.60
especie <- "Balaena mysticetus"

# --- Parametros candidatos R(Balaena mysticetus) ----------------------------

# NIVEL 2 -- dato empirico indirecto
t_lat_B  <- 25.0          # madurez sexual (George et al. 1999)
target_B <- 211.0         # longevidad maxima macho (George et al. 1999)
                           # incertidumbre AAR: ±35 años
Xc_B     <- 3.0           # candidato a universal, igual que Myotis

# NIVEL 3 -- estimados por analogia con Myotis brandtii
# Myotis: eta_basal=0.06, rho_rep=10.0, W_met=0.6
# Correccion por W_met: Balaena W_met=0.35 vs Myotis W_met=0.6
# Factor de correccion: 0.35/0.6 = 0.583
# Se aplica el mismo factor a eta_basal y rho_rep para preservar MS
W_met_B       <- 0.35     # hipometabolismo confirmado (NIVEL 3)
factor_Wmet   <- W_met_B / 0.6  # = 0.583 relativo a Myotis
eta_basal_B   <- 0.06 * factor_Wmet   # = 0.035/año (NIVEL 3)
rho_rep_B     <- 10.0 * factor_Wmet   # = 5.83/año (NIVEL 3)
CV_eta_B      <- 0.30     # prestado de S. purpuratus (NIVEL 3)
eps_B         <- 0.180 * (Xc_B / 17)  # CORREGIDO v2: regla consistente con
                                       # Myotis, tuatara, H.glaber, S.purpuratus
                                       # (antes: Xc_B/4.0, sin explicacion)

# Candidatos a universales
a_B  <- 0.0988
b_B  <- 0.0960
mu_B <- 3.0
n_B  <- 2.0

# Gamma_ext derivado del target
Gamma_ext_B <- log(2) / target_B   # = ln(2)/211 = 0.003285/año (NIVEL 2)

# Verificacion de regimen B antes de simular
X_eq_B       <- eta_basal_B / rho_rep_B
rho_critico_B <- eta_basal_B / (Xc_B * Dc)
MS_B         <- rho_rep_B / rho_critico_B
D_eq_B       <- X_eq_B / Xc_B

cat("=== TMCSA [EM-2] -- Balaena mysticetus -- Regimen B exploratorio v2 ===\n")
cat(sprintf("  eta_basal=%.4f  rho_rep=%.4f  W_met=%.2f (todos NIVEL 3)\n",
    eta_basal_B, rho_rep_B, W_met_B))
cat(sprintf("  X_eq = %.6f  D_eq = %.6f  Dc*Xc = %.3f\n",
    X_eq_B, D_eq_B, Dc * Xc_B))
cat(sprintf("  rho_critico = %.6f  MS = %.0f\n", rho_critico_B, MS_B))
cat(sprintf("  Gamma_ext = %.6f/año (NIVEL 2, derivado de target=%ga)\n\n",
    Gamma_ext_B, target_B))

if (rho_rep_B <= rho_critico_B) {
  cat("  [ERROR] rho_rep <= rho_critico -- NO es regimen B con estos parametros\n")
  stop("Parametros inconsistentes con regimen B")
}

# --- Motor de simulacion regimen B ------------------------------------------
N       <- 300
vida_max <- target_B * 4
semillas <- c(42, 100, 7, 2024, 333)

simular_B_ballena <- function(seed) {
  set.seed(seed)
  dt_sem    <- 7/365
  eps_sem   <- eps_B * sqrt(7)
  sig       <- sqrt(log(1 + CV_eta_B^2))
  muln      <- log(eta_basal_B) - sig^2/2
  Ns        <- round(vida_max * 365/7)

  tc <- numeric(N)
  modo_cese <- character(N)

  for (i in 1:N) {
    eta_i <- rlnorm(1, muln, sig)
    X <- 0.001; V <- 1.0; tf <- NA; modo <- "externo"

    for (sem in 1:Ns) {
      ta <- (sem - 1) * dt_sem

      # Motor regimen B (lineal) -- CORREGIDO v2: eta_basal_B ya incorpora
      # la correccion de W_met (aplicada al derivarlo desde Myotis via
      # factor_Wmet); multiplicar de nuevo por W_met_B aqui era doble conteo.
      dX  <- eta_i - rho_rep_B * X + eps_sem * rnorm(1)
      X   <- max(0, X + dX)
      D   <- X / Xc_B

      # Ecuacion de vitalidad
      base_fric <- 1 - D/Dc
      fric <- if (base_fric >= 0) max(0, base_fric^n_B) else 0
      dV   <- a_B * V * fric - b_B * V^mu_B - max(0, D - Dc)
      V    <- min(1.10, max(0, V + dV * dt_sem))

      # Cese interno
      if (X >= Xc_B || V <= 0.01) {
        tf   <- ta
        modo <- "interno"
        break
      }

      # Cese externo (Bernoulli semanal)
      if (runif(1) < Gamma_ext_B * dt_sem) {
        tf   <- ta
        modo <- "externo"
        break
      }
    }
    tc[i]       <- if (!is.na(tf)) tf else vida_max
    modo_cese[i] <- modo
  }

  list(
    mediana     = median(tc),
    pct_externo = mean(modo_cese == "externo") * 100,
    D_max       = max(sapply(1:min(10, N), function(j) {
      set.seed(seed + j * 1000)
      eta_j <- rlnorm(1, muln, sig)
      Xmax <- 0
      for (s in 1:Ns) {
        dX <- eta_j - rho_rep_B * Xmax + eps_sem * rnorm(1)
        Xmax <- max(0, Xmax + dX)
      }
      Xmax / Xc_B
    }))
  )
}

# --- Correr 5 semillas -------------------------------------------------------
cat("Simulando 5 semillas (N=300 cada una)...\n\n")
resultados <- lapply(semillas, simular_B_ballena)

medianas    <- sapply(resultados, function(r) r$mediana)
pct_ext     <- sapply(resultados, function(r) r$pct_externo)
med_promedio <- mean(medianas)
error_pct   <- abs(med_promedio - target_B) / target_B * 100
D_max_obs   <- max(sapply(resultados, function(r) r$D_max))

for (k in seq_along(semillas)) {
  cat(sprintf("  Semilla %d: mediana=%.1fa  ext=%.0f%%\n",
      semillas[k], medianas[k], pct_ext[k]))
}

cat(sprintf("\n--- PROMEDIO 5 SEMILLAS ---\n"))
cat(sprintf("  Mediana promedio: %.1fa (sd=%.1f)\n",
    med_promedio, sd(medianas)))
cat(sprintf("  Target: %.1fa (George et al. 1999, incertidumbre AAR ±35a)\n",
    target_B))
cat(sprintf("  Error promedio: %.1f%%\n", error_pct))
cat(sprintf("  Cese externo promedio: %.1f%%\n", mean(pct_ext)))
cat(sprintf("  D_max interno observado: %.6f (Dc=%.2f -- %s)\n",
    D_max_obs, Dc,
    if (D_max_obs < Dc) "NO cruza Dc (correcto para regimen B)" else
    "CRUZA Dc (revisar parametros)"))

cat(sprintf("\n  X_eq = %.6f  <<  Dc*Xc = %.3f  (MS = %.0f)\n",
    X_eq_B, Dc * Xc_B, MS_B))

cat("\n[ADVERTENCIA] eta_basal, rho_rep y W_met son NIVEL 3 (candidatos,\n")
cat("  por analogia con Myotis brandtii, corregidos por W_met hipometabolico).\n")
cat("  Este resultado es un aporte CUALITATIVO -- confirma si el diseño\n")
cat("  de regimen B es plausible para esta especie, NO una calibracion\n")
cat("  verificada. Misma categoria que Hydra antes de David & Campbell.\n")
cat("  NOTA: el target de 211a tiene incertidumbre ±35a (metodo AAR).\n")
cat("  Un error aparente dentro de ese rango es metodologicamente aceptable.\n")

# --- Graficas ----------------------------------------------------------------
# Genera trayectoria representativa (semilla 42, primer individuo)
# para los paneles de V(t) y D(t)
fecha_hoy  <- format(Sys.Date(), "%Y%m%d")
nombre_png <- paste0("TMCSA_EM2_Balaena_RegB_v3_", fecha_hoy, ".png")

set.seed(42)
dt_sem_g  <- 7/365
eps_sem_g <- eps_B * sqrt(7)
sig_g     <- sqrt(log(1 + CV_eta_B^2))
muln_g    <- log(eta_basal_B) - sig_g^2/2
Ns_g      <- round(vida_max * 365/7)
eta_rep   <- rlnorm(1, muln_g, sig_g)

traj_t <- numeric(Ns_g)
traj_V <- numeric(Ns_g)
traj_D <- numeric(Ns_g)
X_g <- 0.001; V_g <- 1.0

for (s in 1:Ns_g) {
  ta_g <- (s - 1) * dt_sem_g
  dX_g <- eta_rep - rho_rep_B * X_g + eps_sem_g * rnorm(1)
  X_g  <- max(0, X_g + dX_g)
  D_g  <- X_g / Xc_B
  base_fric_g <- 1 - D_g / Dc
  fric_g <- if (base_fric_g >= 0) max(0, base_fric_g^n_B) else 0
  dV_g <- a_B * V_g * fric_g - b_B * V_g^mu_B - max(0, D_g - Dc)
  V_g  <- min(1.10, max(0, V_g + dV_g * dt_sem_g))
  traj_t[s] <- ta_g
  traj_V[s] <- V_g
  traj_D[s] <- D_g
  if (X_g >= Xc_B || V_g <= 0.01 ||
      runif(1) < Gamma_ext_B * dt_sem_g) break
}
n_pasos <- s

# Usar todas las medianas de las 5 semillas para el histograma
set.seed(42)
tc_histo <- numeric(N)
modo_histo <- character(N)
muln_h <- log(eta_basal_B) - sig_g^2/2
for (i in 1:N) {
  eta_i <- rlnorm(1, muln_h, sig_g)
  X_h <- 0.001; V_h <- 1.0; tf_h <- NA; modo_h <- "externo"
  for (sem in 1:Ns_g) {
    ta_h <- (sem - 1) * dt_sem_g
    dX_h <- eta_i - rho_rep_B * X_h + eps_sem_g * rnorm(1)
    X_h  <- max(0, X_h + dX_h)
    D_h  <- X_h / Xc_B
    bf_h <- 1 - D_h / Dc
    fr_h <- if (bf_h >= 0) max(0, bf_h^n_B) else 0
    dV_h <- a_B * V_h * fr_h - b_B * V_h^mu_B - max(0, D_h - Dc)
    V_h  <- min(1.10, max(0, V_h + dV_h * dt_sem_g))
    if (X_h >= Xc_B || V_h <= 0.01) { tf_h <- ta_h; modo_h <- "interno"; break }
    if (runif(1) < Gamma_ext_B * dt_sem_g) { tf_h <- ta_h; break }
  }
  tc_histo[i]    <- if (!is.na(tf_h)) tf_h else vida_max
  modo_histo[i]  <- modo_h
}
med_histo <- median(tc_histo)

png(nombre_png, width=1500, height=500, res=110)
par(mfrow=c(1,3), mar=c(4.3, 4.3, 3.2, 1.5))

# Panel 1: Distribucion t_cese
hist(tc_histo, breaks=30, col="#DEEAF1", border="white", freq=FALSE,
     xlab="Edad de cese (anos)", ylab="Densidad",
     main=sprintf("Dist. t_cese -- med=%.0fa (100%% ext)", med_histo))
abline(v=med_histo,  lwd=2, col="#1F4E79", lty=1)
abline(v=target_B,   lwd=2, col="#C00000", lty=2)
abline(v=target_B - 35, lwd=1, col="#C00000", lty=3)
abline(v=target_B + 35, lwd=1, col="#C00000", lty=3)
legend("topright",
       legend=c(sprintf("Mediana sim (%.0fa)", med_histo),
                sprintf("Target (%.0fa +/-35a)", target_B)),
       col=c("#1F4E79","#C00000"), lty=c(1,2), lwd=2, cex=0.75, bty="n")

# Panel 2: V(t)
plot(traj_t[1:n_pasos], traj_V[1:n_pasos],
     type="l", lwd=2.5, col="#1A7C4F",
     xlab="Edad (anos)", ylab="V(t) -- Vitalidad",
     main="V(t) -- individuo representativo",
     ylim=c(0, 1.15))
abline(h=1.0, lty=3, col="gray60")
abline(h=0.01, lty=2, col="#C00000", lwd=1)
legend("bottomleft",
       legend=c("V(t)", "V=1.0 (pico)", "Umbral cese V=0.01"),
       col=c("#1A7C4F","gray60","#C00000"),
       lty=c(1,3,2), lwd=c(2.5,1,1), cex=0.75, bty="n")

# Panel 3: D(t)
plot(traj_t[1:n_pasos], traj_D[1:n_pasos],
     type="l", lwd=2.5, col="#C77B00",
     xlab="Edad (anos)", ylab="D(t) = X/Xc",
     main=sprintf("D(t) -- max=%.4f << Dc=%.2f", max(traj_D[1:n_pasos]), Dc),
     ylim=c(0, Dc * 1.1))
abline(h=Dc,   lty=2, col="#C00000", lwd=1.5)
abline(h=D_eq_B, lty=3, col="gray50", lwd=1)
legend("topright",
       legend=c("D(t)", sprintf("Dc=%.2f (umbral colapso)", Dc),
                sprintf("D_eq=%.4f (equilibrio)", D_eq_B)),
       col=c("#C77B00","#C00000","gray50"),
       lty=c(1,2,3), lwd=c(2.5,1.5,1), cex=0.75, bty="n")

dev.off()
cat(sprintf("\nGrafica generada: %s (3 paneles: dist t_cese, V(t), D(t))\n",
    nombre_png))
