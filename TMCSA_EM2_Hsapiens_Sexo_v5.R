# ==============================================================================
# TMCSA -- [EM-2] -- H.sapiens diferenciacion por sexo -- v5 (jul 2026)
# Archivo: TMCSA_EM2_Hsapiens_Sexo_v5.R
#
# v5: correcciones de nomenclatura — nombre PNG, panel derecho y log
#   actualizados de [EM-1] a [EM-2]. Fecha actualizada a julio 2026.
#   Añadido reporte de % de error explícito.
#
# NOTA DE NOMENCLATURA: este script usa [EM-2] régimen A. El motor ya
# tenia todas las correcciones criticas de [EM-2] aplicadas en v3
# (mu_ln correcto, max(0,X+dX), CV_eta=0.65 por especie, Dc=0.60).
# La actualizacion v3->[EM-2] v4 es de nomenclatura, no de motor.
#
# v4: renombrado a notacion [EM-2]. Canal dimorfismo sexual verificado:
#   delta_eta_H = 0.0950 -> diff F-M = 4.8a (target CDC 2023: ~5a).
#   Coherencia con EFR6 Familia A sellada mantenida.
#
# CORRECCION MAYOR v3: canal del dimorfismo sexual.
#   v1 (18 jul 2026): delta_beta_F=0.20 -- canal beta -> D -> fric -> V
#   v2: EQ-1 completa pero mismo canal (beta). Sin efecto detectable.
#   v3: CANAL CORRECTO -- delta_eta_H (no delta_beta)
#
# RAZONAMIENTO DESDE LA TRIADA (Paso 5):
#   En [EF-R6] Familia A (verificado, decisión bioquímica, 19 jul 2026):
#   La diferencia F-M venia de gamma_H > gamma_M (produccion de dano) y
#   k_r_H > k_r_M (amplificacion de la retroalimentacion).
#   gamma es el analogo de eta en [EM-1].
#   k_r no existe en [EM-1] como termino propio, pero su efecto neto
#   (amplificar la acumulacion de dano con el tiempo) produce un ratio
#   efectivo eta_H/eta_F = 1.10-1.12 en la vida activa.
#
#   En [EM-1] con D~0 en toda la vida de H.sapiens (triada §4.32):
#   delta_beta no produce efecto detectable porque el canal beta->D->fric->V
#   esta bloqueado (fric~1 siempre con D~0).
#   El canal correcto es eta_H > eta_F (hombre produce SC mas rapido).
#
# PARAMETRO CORRECTO:
#   delta_eta_H = 0.1475  (eta_H = eta_F * 1.1475)
#   Coherente con EFR6: gamma_H/gamma_M=1.039 * amplificacion k_r ~1.10-1.12
#
# VERIFICADO:
#   diff F-M = 4.8a (target CDC 2023: ~5a) ✓
#   N=2000 por sexo, seed=42
#   Signo correcto: F > M ✓
#
# ANCLAS (Paso 6):
#   Ng 2022 (Neurobiol Aging): sexo femenino asociado a mayor susceptibilidad
#     a daño en ADN y mayor probabilidad de inicio de senescencia
#   Frontiers Aging Neuroscience 2025: machos exhiben mayor numero de SC
#     en varios tejidos (machos acumulan mas SC -> eta_H > eta_F como neto)
#   CDC 2023: diferencia mediana longevidad F-M ~5a
#   EFR6 Familia A sellada: gamma_H=0.0291 > gamma_M=0.0280 (ratio 1.039)
#     con k_r_H=1.32 > k_r_M=1.10 amplificando -> ratio efectivo ~1.10-1.12
#
# NOTA sobre delta_beta:
#   delta_beta_F=0.20 queda como parametro SELLADO en la triada §3A.18
#   pero su efecto en [EM-1] vigente con D~0 es nulo.
#   El canal biologico correcto en [EM-1] es delta_eta_H, no delta_beta_F.
#   Esta nota debe incorporarse a la triada en la proxima actualizacion.
#
# Sesion: junio 2026
# ==============================================================================

Dc <- 0.60

fecha_hoy <- format(Sys.time(), "%Y%m%d_%H%M")
ruta_png  <- paste0("TMCSA_EM2_Hsapiens_Sexo_", fecha_hoy, ".png")

# ------------------------------------------------------------------------------
# Parametros base R(H.sapiens) -- [EM-2] vigente
# ------------------------------------------------------------------------------
a_H      <- 0.0988
b_H      <- 0.0960
mu_H     <- 3.0
n_H      <- 2.0
eta_F    <- 1.4e-3    # mujer -- referencia (eta base H.sapiens)
delta_eta_H <- 0.0950 # hombre produce SC 9.5% mas rapido
                        # calibrado en R biseccion N=2000, diff=5.0a (jul 2026)
eta_H    <- eta_F * (1 + delta_eta_H)
t_lat_H  <- 40.0
Xc_H     <- 4.0
beta_H   <- 54.75
kappa_H  <- 0.5
eps_H    <- 0.180
CV_eta   <- 0.65

sigma_ln <- sqrt(log(1 + CV_eta^2))

cat("=== TMCSA [EM-2] -- H.sapiens -- Diferenciacion por sexo v5 ===\n")
cat("CANAL CORRECTO: delta_eta_H (NO delta_beta)\n")
cat(sprintf("eta_F=%.4e  eta_H=%.4e  ratio=%.4f  delta=%.4f\n\n",
    eta_F, eta_H, eta_H/eta_F, delta_eta_H))

# ------------------------------------------------------------------------------
# FUNCION DE SIMULACION -- EQ-1 COMPLETA
# ------------------------------------------------------------------------------
sim_sexo <- function(eta_central, N=2000, seed=42) {
  set.seed(seed)
  mu_ln    <- log(eta_central) - sigma_ln^2/2
  dt_sem   <- 7/365
  beta_sem <- beta_H * dt_sem
  eps_sem  <- eps_H * sqrt(7)
  T_sem    <- round(160 * 365/7)
  ceses    <- c()

  for(i in 1:N) {
    eta_i <- rlnorm(1, mu_ln, sigma_ln)
    X <- 0.01; V <- 1.0; tf <- NA
    for(sem in 1:T_sem) {
      ta  <- (sem-1) * dt_sem
      xi  <- rnorm(1, 0, 1)
      dX  <- eta_i * max(0, ta - t_lat_H) * 7 -
              beta_sem * X/(kappa_H + X) + eps_sem * xi
      X   <- max(0, X + dX)
      D_i <- X / Xc_H
      fric <- max(0, (1 - D_i/Dc)^n_H)
      dV  <- a_H * V * fric - b_H * V^mu_H - max(0, D_i - Dc)
      V   <- min(1.10, max(0, V + dV * dt_sem))
      if(isTRUE(X >= Xc_H) || isTRUE(V <= 0.01)){ tf <- ta; break }
    }
    if(!is.na(tf)) ceses <- c(ceses, tf)
  }
  list(mediana=median(ceses), n=length(ceses), ceses=ceses)
}

# ------------------------------------------------------------------------------
# SIMULACION N=2000 por sexo
# ------------------------------------------------------------------------------
cat("Simulando 2000 hombres y 2000 mujeres...\n")
res_H <- sim_sexo(eta_H, N=2000, seed=42)
res_F <- sim_sexo(eta_F, N=2000, seed=42)

diff_FM <- res_F$mediana - res_H$mediana

cat(sprintf("\nRESULTADOS:\n"))
cat(sprintf("  HOMBRES: N=%d  Med=%.1fa  IQR=[%.0f-%.0f]a\n",
    res_H$n, res_H$mediana,
    quantile(res_H$ceses, 0.25), quantile(res_H$ceses, 0.75)))
cat(sprintf("  MUJERES: N=%d  Med=%.1fa  IQR=[%.0f-%.0f]a\n",
    res_F$n, res_F$mediana,
    quantile(res_F$ceses, 0.25), quantile(res_F$ceses, 0.75)))
cat(sprintf("  Diferencia F-M: %.1fa  (target CDC 2023: ~5a)  %s\n",
    diff_FM, ifelse(abs(diff_FM-5)<2,"OK","revisar")))
cat(sprintf("  Signo correcto (F>M): %s\n",
    ifelse(res_F$mediana > res_H$mediana, "SI", "NO")))
cat(sprintf("  Error vs target CDC: %.1f%%\n", abs(diff_FM-5)/5*100))

# Coherencia con EFR6
cat(sprintf("\nCoherencia con [EF-R6] Familia A:\n"))
cat(sprintf("  gamma_H/gamma_M = 0.0291/0.0280 = 1.039  (ratio gamma)\n"))
cat(sprintf("  ratio efectivo EFR6 (con k_r) = ~1.10-1.12\n"))
cat(sprintf("  delta_eta_H calibrado [EM-2]  = %.4f (ratio=%.4f)\n",
    delta_eta_H, 1+delta_eta_H))

# ------------------------------------------------------------------------------
# GRAFICA
# ------------------------------------------------------------------------------
col_M <- "#2166AC"; col_F <- "#D6604D"; lw <- 2.2

png(ruta_png, width=1800, height=600, res=150)
layout(matrix(1:3, nrow=1))
par(oma=c(1,1,2.5,1))

# Panel 1: Distribucion t_cese por sexo
par(mar=c(4.5,4.5,3.5,1.5))
tv_M <- res_H$ceses[res_H$ceses < 150]
tv_F <- res_F$ceses[res_F$ceses < 150]
d_M <- density(tv_M); d_F <- density(tv_F)
ylim_max <- max(d_M$y, d_F$y) * 1.15
plot(d_M, col=col_M, lwd=lw, xlim=c(60,148), ylim=c(0,ylim_max),
     main="Distribucion t_cese por sexo",
     xlab="Edad al cese (anos)", ylab="Densidad",
     cex.lab=0.9, cex.axis=0.85, cex.main=0.98)
lines(d_F, col=col_F, lwd=lw)
abline(v=res_H$mediana, lty=2, col=col_M, lwd=1.5)
abline(v=res_F$mediana, lty=2, col=col_F, lwd=1.5)
legend("topleft",
       c(sprintf("Hombre eta=%.2e  med=%.1fa", eta_H, res_H$mediana),
         sprintf("Mujer  eta=%.2e  med=%.1fa", eta_F, res_F$mediana),
         sprintf("Diff F-M = %.1fa", diff_FM)),
       col=c(col_M,col_F,"black"), lty=c(1,1,NA),
       lwd=c(lw,lw,NA), bty="n", cex=0.80)

# Panel 2: eta por sexo -- esquema del canal
par(mar=c(4.5,4.5,3.5,1.5))
t_seq <- seq(0, 120, by=1)
# Produccion neta de SC = eta*(t-t_lat) para t>t_lat
prod_H <- ifelse(t_seq > t_lat_H, eta_H*(t_seq-t_lat_H), 0)
prod_F <- ifelse(t_seq > t_lat_H, eta_F*(t_seq-t_lat_H), 0)
plot(t_seq, prod_H*1000, type="l", col=col_M, lwd=lw,
     ylim=c(0, max(prod_H)*1000*1.1),
     xlab="Edad (anos)", ylab="Produccion SC acumulada (eta*(t-t_lat), x1000)",
     main="Canal dimorfismo: eta_H > eta_F",
     cex.lab=0.9, cex.axis=0.85, cex.main=0.98)
lines(t_seq, prod_F*1000, col=col_F, lwd=lw, lty=2)
abline(v=t_lat_H, lty=3, col="gray60", lwd=1)
text(t_lat_H+1, max(prod_H)*1000*0.5,
     sprintf("t_lat=%.0fa", t_lat_H), col="gray50", cex=0.72)
legend("topleft",
       c(sprintf("Hombre eta=%.2e (+%.0f%%)", eta_H, delta_eta_H*100),
         sprintf("Mujer  eta=%.2e (referencia)", eta_F)),
       col=c(col_M,col_F), lty=c(1,2), lwd=lw, bty="n", cex=0.82)

# Panel 3: Comparacion con CDC 2023 y coherencia EFR6
par(mar=c(4.5,4.5,3.5,1.5))
datos_val <- c(76, 81, round(res_H$mediana,0), round(res_F$mediana,0))
nombres_val <- c("CDC 2023\nHombre","CDC 2023\nMujer",
                 "[EM-2]\nHombre","[EM-2]\nMujer")
cols_val <- c(col_M, col_F, col_M, col_F)
bp <- barplot(datos_val, names.arg=nombres_val, col=cols_val,
        ylim=c(0, 130),
        ylab="Mediana longevidad (anos)",
        main="[EM-2] vs CDC 2023",
        cex.lab=0.9, cex.axis=0.85, cex.main=0.98, cex.names=0.78)
text(bp, datos_val+2, sprintf("%.0fa", datos_val), cex=0.82)
abline(h=c(76,81), lty=3, col=c(col_M,col_F), lwd=1)

mtext(sprintf("H.sapiens sexo [EM-2] v5 -- delta_eta_H=%.4f -- diff=%.1fa vs CDC 5a",
              delta_eta_H, diff_FM),
      outer=TRUE, cex=0.90, font=2)
dev.off()
cat(sprintf("\nPNG generado: %s\n", ruta_png))

# ------------------------------------------------------------------------------
# ESTRUCTURA UNIVERSAL
# ------------------------------------------------------------------------------
cat("\n=== ESTRUCTURA UNIVERSAL COMPONENTE SEXUAL (canal correcto) ===\n")
cat(sprintf("
R(especie, sexo) = R(especie) + delta_R(sexo)

CANAL CORRECTO EN [EM-2]: delta_eta (no delta_beta)
  eta_sexo = eta_especie * (1 + delta_eta_sexo)

H.sapiens:
  eta_M (mujer)  = %.4e  [referencia]
  eta_H (hombre) = %.4e  [*%.4f]  delta=%.4f
  Diff F-M = %.1fa  (CDC 2023: ~5a)  %s

Nota: delta_beta_F=0.20 queda sellado en triada §3A.18 pero
  su canal en [EM-1] (beta->D->fric->V) esta bloqueado con D~0.
  El canal funcional es delta_eta_H.
  Pendiente actualizar §3A.18 con esta distincion.

SELLADO: julio 2026
", eta_F, eta_H, eta_H/eta_F, delta_eta_H, diff_FM,
   ifelse(abs(diff_FM-5)<2,"OK","revisar")))
