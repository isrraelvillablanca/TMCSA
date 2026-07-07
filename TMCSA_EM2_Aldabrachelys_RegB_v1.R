# ==============================================================================
# TMCSA [EM-2] -- Aldabrachelys gigantea (tortuga gigante de Aldabra) -- Regimen B
# CALIBRACION EXPLORATORIA v1 -- aporte CUALITATIVO, NO VERIFICADO
# ==============================================================================
# Molde base: Myotis brandtii (B.EM2-12), mismo molde ya usado para tuatara
# (Aporte 30) y Balaena mysticetus (Aporte 31).
#
# ORIGEN DE CADA PARAMETRO:
#
#   t_cese_target = 192 anos -- NIVEL 2/3, VERIFICACION DE TIPO DISTINTO A
#     BALAENA. Fuente: "Jonathan", tortuga gigante de Santa Elena, edad
#     estimada 192 anos (bioRxiv, Feb 2025, "Epigenomic insights into extreme
#     longevity... Jonathan"). ADVERTENCIA EXPLICITA: esta verificacion es
#     DOCUMENTAL/HISTORICA (fotografias de 1882-1886 + registros del
#     gobernador que lo recibio ya adulto), NO bioquimica como el metodo de
#     racemizacion de acido aspartico (AAR) usado por George et al. (1999)
#     para Balaena mysticetus. AnAge/genomics.senescence.info declara
#     explicitamente que esta edad "has not been confirmed" en su fuente
#     primaria original -- el respaldo del preprint 2025 es corroboracion
#     circunstancial (documentos historicos), no una medicion de reloj
#     biologico. Nivel de evidencia mas debil que el target de Balaena.
#     Alternativa NO usada aqui: esperanza de vida poblacional via ZIMS +
#     analisis de supervivencia bayesiano, 80.47a (media) / 93.68a (limite
#     superior, +-5.64 SD) -- un numero demografico, no un record individual,
#     de naturaleza distinta a la usada para tuatara/Balaena. Se prefiere
#     Jonathan (192a) por consistencia de tipo de dato (record maximo
#     individual) con los otros Aportes de regimen B, pese a su verificacion
#     mas debil.
#
#   W_met = 0.85 -- NIVEL 3, CANDIDATO, CON TRES CAPAS DE INCERTIDUMBRE
#     APILADAS, TODAS DECLARADAS EXPLICITAMENTE:
#     (1) Formula fuente: Hughes, Gaymer, Moore, Woakes (1971), J Exp Biol
#         55(3):651-665 -- consumo de O2 medido en 9 ejemplares reales de
#         Testudo gigantea (=Aldabrachelys gigantea), animales inactivos:
#         O2(ml/h) = 45.47 * W^0.82. UNIDAD DE W NO CONFIRMADA por el texto
#         completo (de pago, no accesible) -- se infirio W en kg por
#         consistencia con la cifra de "82 ml/kg/h para animales de 100g"
#         que el propio abstract reporta para la escala de actividad.
#     (2) Curva de referencia: en vez de comparar contra mamifero (lo que
#         daria W_met~0.11, no comparable con el criterio usado para
#         tuatara), se comparo contra la razon general ectotermo/endotermo
#         de la literatura (~10-15%), dando W_met = 0.11 / [0.10,0.15] =
#         [0.73, 1.10]. Se uso el punto medio conservador, W_met=0.85.
#     (3) No se encontro ninguna comparacion publicada DIRECTA de
#         Aldabrachelys contra otro reptil/tortuga de referencia -- el
#         calculo depende de una heuristica general de la literatura, no de
#         un dato especifico de la especie como el 55% del tuatara (Wilson &
#         Lee 1970, medicion directa contra lagarto).
#
#   HALLAZGO IMPORTANTE, DECLARADO EXPLICITAMENTE -- a diferencia del
#     tuatara (W_met=0.55, hipometabolismo real y medido incluso frente a
#     otros ectotermos), el rango calculado aqui (0.73-1.10) practicamente
#     incluye 1.0 -- es decir, Aldabrachelys gigantea NO muestra evidencia de
#     metabolismo deprimido respecto a otros ectotermos de su tamano. Esto
#     significa que si el diseno de regimen B resulta plausible para esta
#     especie, SU LONGEVIDAD EXTREMA PROBABLEMENTE NO ESTA MEDIADA POR
#     HIPOMETABOLISMO, a diferencia del tuatara. Explicaciones alternativas
#     no verificadas aqui: (a) un mecanismo de SASP reducido similar al
#     sugerido para Balaena mysticetus (Gros et al. 2025), consistente con
#     los hallazgos reales de metilacion diferencial en el epigenoma de
#     Jonathan (bioRxiv 2025) -- aunque esos datos no se usaron para derivar
#     ningun valor numerico aqui; o (b) ausencia de depredadores naturales en
#     el atolon de Aldabra, un factor ecologico, no biologico-molecular. NO
#     SE DEBE ASUMIR MECANICAMENTE QUE EL MISMO SUPUESTO DE HIPOMETABOLISMO
#     DEL TUATARA APLICA AQUI -- es una diferencia real entre especies, no
#     un error de calculo.
#
#   eta_basal, rho_rep -- NIVEL 3, CANDIDATOS. Mismo procedimiento que
#     tuatara y Balaena: factor_Wmet = W_met_A/0.6 = 1.4167 (NOTESE: mayor a
#     1, a diferencia de tuatara y Balaena donde el factor fue menor a 1 --
#     consecuencia directa de que W_met_A=0.85 > W_met_Myotis=0.6). Aplicado
#     al mismo par de valores de Myotis (0.06, 10.0), preservando el Margen
#     de Seguridad por construccion:
#       eta_basal_A = 0.06 * 1.4167 = 0.0850  (MAYOR que el de Myotis, a
#         diferencia de tuatara/Balaena donde salio menor)
#       rho_rep_A   = 10.0 * 1.4167 = 14.1667 (tambien mayor, misma proporcion)
#       X_eq = eta_basal_A/rho_rep_A = 0.006000 (IDENTICO a Myotis, por
#         construccion)
#       MS = Dc*Xc/X_eq = 300 (IDENTICO a Myotis y a los otros dos Aportes,
#         confirmado antes de correr este script)
#
#   Xc=3.0, eps=0.180*(3.0/17)=0.03176 -- heredados del molde Myotis, misma
#     regla de eps ya corregida y consistente con tuatara/Balaena/S.purpuratus.
#   CV_eta=0.30 -- NIVEL 3, supuesto prestado de S. purpuratus (mismo criterio
#     que tuatara/Balaena, no medido para esta especie).
#   a=0.0988, b=0.0960, mu=3.0, n=2.0 -- heredados del molde Myotis (B.EM2-12).
#
# APORTE 32 (Con soporte parcial, esperado) -- misma categoria que Hydra
#   (Aporte 28), H. glaber (Aporte 29), tuatara (Aporte 30), Balaena
#   (Aporte 31).
#
# GRAFICAS Y % DE ERROR -- incluidos por convencion establecida (jul 2026):
#   los scripts de calibracion/exploracion deben generar graficas de
#   verificacion visual y reportar el % de error explicitamente, sin que
#   esto implique subir la imagen a ningun documento salvo que se requiera.
#
# REGLA B.12 -- parsear antes de correr:
#   Rscript -e 'invisible(parse(file="TMCSA_EM2_Aldabrachelys_RegB_v1.R"))'
# ==============================================================================

Dc <- 0.60

# ------------------------------------------------------------------------------
# R(Aldabrachelys gigantea) -- Regimen B, exploratorio v1
# ------------------------------------------------------------------------------
W_met_A     <- 0.85
factor_Wmet <- W_met_A / 0.6
eta_basal_A <- 0.06 * factor_Wmet
rho_rep_A   <- 10.0 * factor_Wmet
Xc_A        <- 3.0
eps_A       <- 0.180 * (Xc_A / 17)
CV_eta_A    <- 0.30
a_A <- 0.0988; b_A <- 0.0960; mu_A <- 3.0; n_A <- 2.0

target_A    <- 192.0   # Jonathan -- NIVEL 2/3, verificacion documental, no bioquimica
Gamma_ext_A <- log(2) / target_A

X_eq_A <- eta_basal_A / rho_rep_A
MS_A   <- Dc * Xc_A / X_eq_A

cat("=== TMCSA [EM-2] -- Aldabrachelys gigantea -- Regimen B exploratorio v1 ===\n")
cat(sprintf("  eta_basal=%.4f  rho_rep=%.4f  W_met=%.2f (todos NIVEL 3)\n",
            eta_basal_A, rho_rep_A, W_met_A))
cat(sprintf("  X_eq = %.6f  Dc*Xc = %.3f\n", X_eq_A, Dc*Xc_A))
cat(sprintf("  MS = %.0f\n", MS_A))
cat(sprintf("  Gamma_ext = %.6f/anio (NIVEL 2/3, derivado de target=%.0fa, verificacion documental)\n",
            Gamma_ext_A, target_A))

# ------------------------------------------------------------------------------
# MOTOR -- identico en estructura al de Balaena mysticetus v2 (regimen B
# lineal, Gamma_ext primero en el bucle, sin doble conteo de W_met -- ya
# absorbido en eta_basal_A/rho_rep_A).
# ------------------------------------------------------------------------------
simular_aldabrachelys <- function(N_sto=300, seed=42) {
  dt_sem   <- 7/365
  eps_sem  <- eps_A * sqrt(7)
  Gext_sem <- Gamma_ext_A * dt_sem
  sig      <- sqrt(log(1 + CV_eta_A^2))
  muln     <- log(eta_basal_A) - sig^2/2
  vida_max <- target_A * 4
  N_s      <- round(vida_max * 365/7)

  set.seed(seed)
  t_cese_sto <- rep(NA, N_sto)
  causa_ext  <- rep(FALSE, N_sto)
  D_max_obs  <- 0

  for (sim in 1:N_sto) {
    eta_i <- rlnorm(1, muln, sig)
    X <- 0.001; V_i <- 1.0; t_fin <- NA; fue_ext <- FALSE

    for (sem in 1:N_s) {
      ta <- (sem - 1) * dt_sem

      if (runif(1) < Gext_sem) { t_fin <- ta; fue_ext <- TRUE; break }

      dX  <- eta_i - rho_rep_A * X + eps_sem * rnorm(1)
      X   <- max(0, X + dX)
      D_i <- X / Xc_A
      D_max_obs <- max(D_max_obs, D_i)

      base_f <- 1 - D_i/Dc
      fric <- if (base_f > 0) base_f^n_A else 0
      dV  <- a_A * V_i * fric - b_A * V_i^mu_A - max(0, D_i - Dc)
      V_i <- min(1.10, max(0, V_i + dV * dt_sem))

      if (X >= Xc_A || V_i <= 0.01) { t_fin <- ta; fue_ext <- FALSE; break }
    }
    t_cese_sto[sim] <- if (!is.na(t_fin)) t_fin else vida_max
    causa_ext[sim]  <- fue_ext
  }
  list(t_cese=t_cese_sto, causa_ext=causa_ext, D_max_obs=D_max_obs)
}

cat("\nSimulando 5 semillas (N=300 cada una)...\n\n")
semillas  <- c(42, 100, 7, 2024, 333)
medianas  <- numeric(length(semillas))
pcts_ext  <- numeric(length(semillas))
D_max_global <- 0

for (i in seq_along(semillas)) {
  res <- simular_aldabrachelys(N_sto=300, seed=semillas[i])
  med <- median(res$t_cese)
  pext <- mean(res$causa_ext) * 100
  medianas[i] <- med; pcts_ext[i] <- pext
  D_max_global <- max(D_max_global, res$D_max_obs)
  cat(sprintf("  Semilla %d: mediana=%.1fa  ext=%.0f%%\n", semillas[i], med, pext))
}

prom_med <- mean(medianas)
sd_med   <- sd(medianas)
error_pct <- abs(prom_med - target_A) / target_A * 100

cat(sprintf("\n--- PROMEDIO 5 SEMILLAS ---\n"))
cat(sprintf("  Mediana promedio: %.1fa (sd=%.1f)\n", prom_med, sd_med))
cat(sprintf("  Target: %.1fa (Jonathan, verificacion documental 1882-1886 + registros de gobernador)\n", target_A))
cat(sprintf("  Error promedio: %.1f%%\n", error_pct))
cat(sprintf("  Cese externo promedio: %.1f%%\n", mean(pcts_ext)))
cat(sprintf("  D_max interno observado: %.6f (Dc=%.2f -- %s)\n",
            D_max_global, Dc,
            ifelse(D_max_global < Dc, "NO cruza Dc (correcto para regimen B)", "cruza Dc -- revisar")))
cat(sprintf("\n  X_eq = %.6f  <<  Dc*Xc = %.3f  (MS = %.0f)\n", X_eq_A, Dc*Xc_A, MS_A))

cat("\n[ADVERTENCIA] eta_basal, rho_rep y W_met son NIVEL 3 (candidatos,\n")
cat("  por analogia con Myotis brandtii, corregidos por W_met calculado por\n")
cat("  heuristica, con tres capas de incertidumbre declaradas en el encabezado).\n")
cat("  Este resultado es un aporte CUALITATIVO -- confirma si el diseño\n")
cat("  de regimen B es plausible para esta especie, NO una calibracion\n")
cat("  verificada. Misma categoria que Hydra, tuatara y Balaena antes de\n")
cat("  dato directo.\n")
cat("[NOTA] El target de 192a (Jonathan) tiene verificacion DOCUMENTAL/\n")
cat("  HISTORICA, no bioquimica -- mas debil que la de Balaena (AAR).\n")
cat("[NOTA] W_met~0.85-1.0 sugiere que, a diferencia del tuatara, la\n")
cat("  longevidad de esta especie NO esta mediada por hipometabolismo --\n")
cat("  el mecanismo real (SASP reducido, ausencia de depredadores, u otro)\n")
cat("  permanece sin identificar directamente.\n")

# ------------------------------------------------------------------------------
# GRAFICAS -- convencion establecida: siempre incluir verificacion visual.
# ------------------------------------------------------------------------------
res_graf <- simular_aldabrachelys(N_sto=300, seed=42)

png("TMCSA_EM2_Aldabrachelys_RegB_v1_20260707.png", width=1500, height=500)
par(mfrow=c(1,3))

hist(res_graf$t_cese, breaks=20, col="lightblue",
     main=sprintf("Dist. t_cese -- med=%.0fa (%.0f%% ext)", median(res_graf$t_cese), mean(res_graf$causa_ext)*100),
     xlab="Edad de cese (anos)", ylab="Densidad", freq=FALSE)
abline(v=median(res_graf$t_cese), col="blue", lwd=2)
abline(v=target_A, col="red", lwd=2, lty=2)
legend("topright", legend=c(sprintf("Mediana sim (%.0fa)", median(res_graf$t_cese)),
                             sprintf("Target (%.0fa, doc.)", target_A)),
       col=c("blue","red"), lty=c(1,2), lwd=2, cex=0.8)

# Trayectoria individual representativa
dt_sem <- 7/365; eps_sem <- eps_A*sqrt(7); Gext_sem <- Gamma_ext_A*dt_sem
sig <- sqrt(log(1+CV_eta_A^2)); muln <- log(eta_basal_A)-sig^2/2
set.seed(7)
eta_i <- rlnorm(1, muln, sig)
X <- 0.001; V_i <- 1.0; N_max <- round((target_A*4)*365/7)
Vs <- numeric(N_max); Ds <- numeric(N_max); tf_idx <- N_max
for (sem in 1:N_max) {
  ta <- (sem-1)*dt_sem
  if (runif(1) < Gext_sem) { tf_idx <- sem; break }
  dX <- eta_i - rho_rep_A*X + eps_sem*rnorm(1)
  X <- max(0, X+dX); D_i <- X/Xc_A
  base_f <- 1-D_i/Dc; fric <- if(base_f>0) base_f^n_A else 0
  dV <- a_A*V_i*fric - b_A*V_i^mu_A - max(0,D_i-Dc)
  V_i <- min(1.10, max(0, V_i+dV*dt_sem))
  Vs[sem] <- V_i; Ds[sem] <- D_i
  if (X>=Xc_A || V_i<=0.01) { tf_idx <- sem; break }
}
t_ejes <- (1:tf_idx)*dt_sem
plot(t_ejes, Vs[1:tf_idx], type="l", col="darkgreen", lwd=2,
     main="V(t) -- individuo representativo", xlab="Edad (anos)", ylab="V(t) -- Vitalidad",
     ylim=c(0,1.1))
abline(h=1.0, col="gray", lty=3); abline(h=0.01, col="red", lty=2)
legend("bottomleft", legend=c("V(t)","V=1.0 (pico)","Umbral cese V=0.01"),
       col=c("darkgreen","gray","red"), lty=c(1,3,2), cex=0.8)

plot(t_ejes, Ds[1:tf_idx], type="l", col="darkorange", lwd=2,
     main=sprintf("D(t) -- max=%.4f << Dc=%.2f", max(Ds[1:tf_idx]), Dc),
     xlab="Edad (anos)", ylab="D(t) = X/Xc", ylim=c(0, Dc*1.1))
abline(h=Dc, col="black", lty=2); abline(h=X_eq_A, col="gray", lty=3)
legend("topright", legend=c("D(t)", sprintf("Dc=%.2f (umbral colapso)",Dc),
                             sprintf("D_eq=%.4f (equilibrio)",X_eq_A)),
       col=c("darkorange","black","gray"), lty=c(1,2,3), cex=0.8)

dev.off()
cat("\nGrafica generada: TMCSA_EM2_Aldabrachelys_RegB_v1_20260707.png (3 paneles: dist t_cese, V(t), D(t))\n")
