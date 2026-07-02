# ==============================================================================
# TMCSA -- [EM-2] -- Script canonico 10A -- v7 -- once especies verificadas
# Archivo: TMCSA_canonico_10A_EM2_v7.R
#
# CORRECCION v7 (jul 2026) -- asimetria beta en rama DETERMINISTA:
#   Hallazgo por identificabilidad estructural (analisis simbolico previo a
#   sensibilidad), confirmado en R real con diagnostico determinista vs
#   promedio de 300 trayectorias estocasticas (H. sapiens): la rama
#   ESTOCASTICA (la que decide t_cese, la unica que se compara contra el
#   target) opera en X crudo con beta*X/(kappa_X+X) -- correcta, sin cambios.
#   La rama DETERMINISTA (solo usada para graficar paneles D(t)/V(t) y para
#   calcular D_max_interno en simular_olivo) operaba en D=X/Xc con
#   R$beta*Da/(kappa_D+Da), SIN dividir beta por Xc -- inconsistente con la
#   sustitucion algebraica real (D=X/Xc exige beta_D = beta_X/Xc). Verificado
#   con diagnostico grafico: la curva D_d determinista subestimaba el dano
#   real por casi un orden de magnitud frente al promedio de 300 trayectorias
#   estocasticas (H. sapiens, D_d~0.05 vs D_promedio~0.20-0.25 a los 140-160a).
#   CORREGIDO: dos lineas (simular_especie linea ~101, simular_olivo linea
#   ~506) ahora usan (R$beta/R$Xc) y (beta_ef_det/R$Xc) respectivamente.
#   NO CAMBIA: ningun t_cese, ningun error%, ningun parametro de R(especie)
#   de ninguna especie -- el error nunca estuvo en la rama que calibra.
#   PENDIENTE DE RE-VERIFICAR: D_max_interno de Olea, Fitzroya, H.glaber y
#   Chelonoidis (calculado via simular_olivo) probablemente estaba
#   subestimado por el mismo motivo. Con el margen previo tan amplio frente
#   a Dc=0.60 (peor caso Fitzroya: 0.111), es improbable que la conclusion
#   cualitativa "cese externo, no interno" cambie -- pero debe confirmarse
#   corriendo este script v7 y mirando los cuatro valores de D_max_interno,
#   no asumirse.
#
# NOTA DE NOMENCLATURA: toda especie de este script fue calibrada con el
# metodo de [EM-2]. Las once confirmaron REGIMEN A (reparacion saturante)
# -- coinciden numericamente con la fase historica llamada "[EM-1]" porque
# su biologia real tiene reparacion con techo, no porque usen una ecuacion
# distinta. Ver Manual operativo de [EM-2] SS1 y Protocolo de calibracion
# canonica para la distincion completa entre regimen y version de ecuacion.
#
# Especies: H.sapiens, Mus, Somniosus, Hydra, Olea, Heterocephalus,
#            Fitzroya, Turritopsis, E.coli, Pan troglodytes, Chelonoidis
#
# ECUACION MADRE (forma universal -- no cambia entre especies):
#   dV/dt     = a*V*(1-D/Dc)^n - b*V^mu - max(0, D-Dc)
#   dD/dt_bio = eta*max(0,t_bio-t_lat) - beta*D/(kappa+D) + eps*xi(t_bio)
#   D         = X / Xc
#   t_bio     = t_cal * W_met
#
# Unico parametro universal verificado: Dc = 0.60
# Ancla molecular: umbral NF-kappaB/RAGE en eje AGE-RAGE
# Todos los demas parametros son de especie -- entregados por R(especie)
#
# CORRECCIONES v6:
#   1. CV_eta leido desde R$CV_eta (no hardcodeado)
#   2. kappa_X = 0.5 por defecto (pendiente verificacion por especie)
#   3. W_met aplicado al tiempo antes de max(0, t_bio - t_lat)
#
# Fuentes: TMCSA_Cap4_v1.docx | TMCSA_Fichas_Variables_EM1_completo.md
#          TMCSA_Acta_Sesion_21jun2026.md | sesion jun 2026
# ==============================================================================
# ==============================================================================
# AUDITORIA DE HERENCIAS v4 (24jun2026) -- 10A
# CV_eta PASO 5 por especie:
#   H.sapiens 0.65 = FUENTE propia (no herencia). Mus 0.20, Hydra 0.10 = propios.
#   Somniosus/Olea/Fitzroya/Pan/Chelonoidis 0.65 = coincide-DIVERSIDAD real.
#   E.coli 0.65 = coincide pero por heterogeneidad polo viejo/nuevo (asimetria),
#     NO por diversidad genetica (Stewart 2005). Razon distinta.
#   Turritopsis 0.10 = clonal laboratorio + rejuvenece como Hydra (NIVEL3).
#   H.glaber 0.20 = colonia laboratorio (ya corregido).
# eps: H.glaber escalado 0.180*(Xc/17)=0.0684 (antes heredado 0.180).
#   HALLAZGO: eps NO escala universal con Xc. Somniosus se DEJA en eps=0.180
#   (verificado 4.0%); al escalarlo a 0.021 el error salta a 280% (casi no
#   colapsa). El escalado eps=0.180*(Xc/17) funciona en unas especies pero
#   NO en todas -- es coincidencia por especie, no ley universal. Cada eps
#   se verifica corriendo, no se asume por formula.
# a,b=0.0988/0.0960: COINCIDENCIA ratio conservado, pendiente metabolismo. NO universal.
# Dc=0.60: geometria del 60% del Xc PROPIO, NO herencia.
# ==============================================================================


graphics.off()

# PNG con fecha en directorio de trabajo -- no requiere ruta absoluta
fecha_hoy  <- format(Sys.time(), "%Y%m%d_%H%M")
ruta_png   <- paste0("TMCSA_EM1_10A_", fecha_hoy, ".png")

# ------------------------------------------------------------------------------
# PARAMETRO UNIVERSAL
# ------------------------------------------------------------------------------
Dc <- 0.60  # Universal verificado: umbral NF-kappaB/RAGE

# ------------------------------------------------------------------------------
# FUNCION DE SIMULACION -- recibe R(especie)
# ------------------------------------------------------------------------------
simular_especie <- function(R, N_sto=500, seed=42, T_max_det=NULL) {

  # kappa_X: pendiente verificacion directa por especie
  # Karin 2019 midio kappa=0.5 en raton (Xc=17%)
  # Se usa 0.5 para todas las especies hasta tener medicion directa
  kappa_X <- if(!is.null(R$kappa)) R$kappa else 0.5
  kappa_D <- kappa_X / R$Xc

  dt_sem   <- 7/365
  beta_sem <- R$beta * 7/365
  eps_sem  <- R$eps  * sqrt(7)
  eta_det  <- R$eta  * 365 / R$Xc

  stopifnot("P1 violada: a > b requerido" = R$a > R$b)

  CV_eta   <- R$CV_eta
  sigma_ln <- sqrt(log(1 + CV_eta^2))
  # mu_ln: convencion media=eta (calibrada en v5, error <2% en todas las especies)
  mu_ln    <- log(R$eta) - sigma_ln^2/2

  set.seed(seed)

  # Simulacion determinista
  if(is.null(T_max_det)) T_max_det <- R$t_cese_target * 2.5
  dt  <- 1/365
  N_t <- round(T_max_det/dt)
  t_d <- seq(0, T_max_det, length.out=N_t)
  V_d <- rep(1.0, N_t); D_d <- rep(0.0, N_t); dVdt_d <- rep(0.0, N_t)

  for(i in 2:N_t){
    ta <- t_d[i-1]; Va <- V_d[i-1]; Da <- D_d[i-1]
    if(isTRUE(Va <= 0.001)){ V_d[i]=0; D_d[i]=Da; next }
    t_bio_a <- ta * R$W_met  # CORRECCION v6: W_met aplicado al tiempo
    fric <- max(0, (1 - Da/Dc)^R$n)
    dD   <- eta_det * max(0, t_bio_a - R$t_lat) - (R$beta/R$Xc) * Da/(kappa_D + Da)  # CORREGIDO v7: beta/Xc (ver nota v7 en cabecera)
    dV   <- R$a * Va * fric - R$b * Va^R$mu - max(0, Da - Dc)
    dVdt_d[i] <- dV
    V_d[i] <- min(1.10, max(0, Va + dV*dt))
    D_d[i] <- max(Da, Da + dD*dt)  # D monotonicamente no decreciente
  }

  # Simulacion estocastica
  T_sto <- R$t_cese_target * 4
  N_s   <- round(T_sto * 365/7)
  t_s   <- seq(0, T_sto, length.out=N_s)
  t_cese_sto <- rep(NA, N_sto); T_activo <- rep(NA, N_sto)
  N_tray <- 30; V_tray <- matrix(NA, N_tray, N_s)

  for(sim in 1:N_sto){
    eta_i <- rlnorm(1, mu_ln, sigma_ln)
    X <- 0.01; V_i <- 1.0; V_acum <- 0; t_fin <- NA
    Vs <- if(sim <= N_tray) rep(NA, N_s) else NULL
    for(sem in 1:N_s){
      ta       <- t_s[sem]
      xi       <- rnorm(1, 0, 1)
      t_bio_sem <- ta * R$W_met  # CORRECCION v6
      dX <- eta_i * max(0, t_bio_sem - R$t_lat) * 7 -
            beta_sem * X/(kappa_X + X) + eps_sem * xi
      X  <- max(0, X + dX)
      D_i  <- X / R$Xc
      fric <- max(0, (1 - D_i/Dc)^R$n)
      dV   <- R$a * V_i * fric - R$b * V_i^R$mu - max(0, D_i - Dc)
      V_i  <- min(1.10, max(0, V_i + dV*dt_sem))
      V_acum <- V_acum + V_i*dt_sem
      if(!is.null(Vs)) Vs[sem] <- V_i
      if(isTRUE(X >= R$Xc) || isTRUE(V_i <= 0.01)){
        t_fin <- ta; break
      }
    }
    t_cese_sto[sim] <- t_fin; T_activo[sim] <- V_acum
    if(sim <= N_tray) V_tray[sim,] <- Vs
  }

  tv <- t_cese_sto[!is.na(t_cese_sto)]

  # Reporte diferenciado segun tipo de especie
  if(!is.null(R$verificar_sin_colapso) && R$verificar_sin_colapso) {
    # Negligible senescence: criterio es ausencia de colapso en t_obs anos
    t_obs <- R$t_obs_verificado
    n_ok  <- sum(t_cese_sto >= t_obs | is.na(t_cese_sto), na.rm=TRUE)
    cat(sprintf("  %s: sin colapso en %.0fa = %d/%d (%.1f%%) -- negligible senescence verificado\n",
        R$nombre, t_obs, n_ok, N_sto, 100*n_ok/N_sto))
  } else {
    cat(sprintf("  %s: muertes=%d/%d  mediana=%.2fa  target=%.1fa  error=%.1f%%\n",
        R$nombre, length(tv), N_sto,
        median(tv), R$t_cese_target,
        abs(median(tv)-R$t_cese_target)/R$t_cese_target*100))
  }

  list(R=R, t_d=t_d, V_d=V_d, D_d=D_d, dVdt_d=dVdt_d,
       t_s=t_s, V_tray=V_tray, T_activo=T_activo,
       t_cese_sto=t_cese_sto, mediana=median(tv))
}

# ------------------------------------------------------------------------------
# R(H. sapiens)
# Fuentes: Karin 2019, Idda 2020, Liu 2009, Calment 122.5a
# ------------------------------------------------------------------------------
R_Hsapiens <- list(
  nombre        = "H. sapiens",
  t_cese_target = 122.0,
  a=0.0988, b=0.0960, mu=3.0, n=2.0,
  eta=1.4e-3, t_lat=40.0, Xc=4.0, beta=54.75, eps=0.180,
  W_met=1.0, delta_beta=0.0, rho=0.0,
  CV_eta=0.65  # Idda 2020 -- poblacion humana heterogenea
)

# ------------------------------------------------------------------------------
# R(Mus musculus) -- cepa C57BL/6J
# Fuentes: Karin 2019, Schultz 2020
# CV_eta=0.20 derivado de IQR real B6 -- correccion principal v6
# ------------------------------------------------------------------------------
R_Musculus <- list(
  nombre        = "Mus musculus",
  t_cese_target = 4.0,
  a=0.0988, b=0.0960, mu=3.0, n=2.0,
  eta=0.084, t_lat=1.31, Xc=17.0, beta=54.75, eps=0.16,
  W_met=1.0, delta_beta=0.0, rho=0.0,
  CV_eta=0.20  # Schultz 2020 -- cepa consanguinea, variabilidad real B6
)

# ------------------------------------------------------------------------------
# R(Somniosus microcephalus) -- Tiburon de Groenlandia
# Modo 1 K basico extremo
#
# ORIGEN DE CADA PARAMETRO:
# t_cese_target=392a: Nielsen 2016 (Science), radiocarbono cristalino ocular.
#   Nivel NIVEL 1 -- dato real verificado.
# t_lat=156a: madurez sexual Nielsen 2016 -- proxy de t_lat.
#   Razon t_lat/t_vida=0.40, comparable a H.sapiens (0.33).
# Xc=2.0281%: formula alometrica EQ-3: 17*(2.5/392)^0.4206
#   Ancla: Karin 2019 (raton Xc=17%) + Idda 2020 (humano Xc=4%).
#   Interpretacion: umbral mas estricto que H.sapiens -- longevidad extrema.
# eta=2.7813e-05: calibracion inversa Monte Carlo N=400, error 1.4% (jun 2026).
#   Nota de escalamiento: eta_H/eta_S=50.3x para longevidad 3.2x mayor.
#   No es escala lineal ni cuadratica -- emerge de la interaccion eta y Xc.
#   Xc=2.03% (umbral estricto) requiere eta muy bajo para no colapsar antes.
#   Pendiente: medicion directa de SC en tejido de Somniosus.
# mu=3.0, n=2.0: aproximacion desde mamiferos.
#   NOTA: Somniosus es condrictio (cartilago, no hueso). Pendiente verificacion.
# beta=54.75, eps=0.180, CV_eta=0.65: referencia H.sapiens. Pendiente medicion.
# ------------------------------------------------------------------------------
R_Somniosus <- list(
  nombre        = "Somniosus microcephalus",
  t_cese_target = 392.0,
  a=0.0988, b=0.0960, mu=3.0, n=2.0,
  eta=2.7813e-05, t_lat=156.0, Xc=2.0281, beta=54.75, eps=0.180,
  W_met=1.0, delta_beta=0.0, rho=0.0,
  CV_eta=0.65  # PASO5: salvaje artico diverso -- coincide-diversidad
)

# ------------------------------------------------------------------------------
# R(Hydra vulgaris)
# Negligible senescence -- mecanismo de renovacion celular continua
#
# MECANISMO EN [EM-1] -- leer antes de modificar:
#   Hydra renueva todos sus tipos celulares en ~20 dias via tres poblaciones
#   de celulas madre autorrenovables (Bradshaw et al. PMC9770530).
#   No existe soporte molecular para acumulacion de SC.
#   Implementacion: t_lat=1000a >> t_observacion.
#   El interruptor max(0, t_bio - t_lat) = 0 siempre.
#   El valor de eta es irrelevante -- t_lat lo neutraliza.
#   NO implementar como eta~0 -- el mecanismo es t_lat, no eta.
#
# ORIGEN DE CADA PARAMETRO:
# t_cese_target=1e4: simbolico. Criterio real: sin colapso en 41a (Martinez 1998).
# t_lat=1000a: efectivamente infinito -- interruptor nunca activo.
#   Verificado: 0 colapsos en N=300 en 200a (jun 2026).
# eta=1e-6: simbolico -- no determina el resultado.
# eps=0.01: muy bajo -- Hydra de laboratorio es clonal.
# mu=1.0: sin ECM densa -- diblastico, sin tejido conectivo organizado.
# n=1.5: redundancia moderada -- pendiente verificacion directa.
# Xc=4.0: no relevante con t_lat=1000a. Referencia H.sapiens.
# CV_eta=0.10: linea clonal -- baja variabilidad genetica.
# verificar_sin_colapso=TRUE: criterio de verificacion es Martinez 1998 (41a).
# ------------------------------------------------------------------------------
R_Hydra <- list(
  nombre        = "Hydra vulgaris",
  t_cese_target = 1e4,
  a=0.0988, b=0.0960, mu=1.0, n=1.5,
  eta=1e-6, t_lat=1000.0, Xc=4.0, beta=54.75, eps=0.01,
  W_met=1.0, delta_beta=0.0, rho=0.0,
  CV_eta=0.10,
  verificar_sin_colapso = TRUE,
  t_obs_verificado      = 41.0
)

# ------------------------------------------------------------------------------
# CORRER LAS CUATRO ESPECIES
# ------------------------------------------------------------------------------
cat("=== TMCSA [EM-1] -- Canonico 10A -- once especies verificadas ===\n")
cat("Correcciones v6: CV_eta por especie | W_met en tiempo | kappa_X pendiente\n\n")

res_H <- simular_especie(R_Hsapiens,  N_sto=500, T_max_det=160)
res_M <- simular_especie(R_Musculus,  N_sto=500, T_max_det=12)
res_S <- simular_especie(R_Somniosus, N_sto=300, T_max_det=600)
res_Y <- simular_especie(R_Hydra,     N_sto=30,  T_max_det=45)

# ------------------------------------------------------------------------------
# GRAFICAS -- H.sapiens y Mus (las dos especies con graficas completas)
# ------------------------------------------------------------------------------
col_H <- "#1D9E75"; col_M <- "#E85D24"; lw <- 2.2

# Grafica PNG -- se guarda en el directorio de trabajo de R
png(ruta_png, width=2100, height=1500, res=150)
  layout(matrix(1:6, nrow=2, byrow=TRUE))
  par(oma=c(1,1,2.5,1))

  # Panel 1: V(t) H.sapiens
  par(mar=c(4.5,4.5,3.5,1.5))
  idx_max_H <- which.max(res_H$V_d)
  plot(res_H$t_d, res_H$V_d, type="l", col=col_H, lwd=lw,
       xlim=c(0,160), ylim=c(0.997, res_H$V_d[idx_max_H]+0.004),
       xlab="Edad (anos)", ylab="V(t)",
       main="V(t) H. sapiens -- tres etapas",
       cex.lab=0.9, cex.axis=0.85, cex.main=0.98)
  abline(h=1.0, lty=2, col="gray60", lwd=1)
  abline(v=R_Hsapiens$t_lat, lty=3, col="gray60", lwd=1)
  mtext("construccion", side=3, at=18,  line=-2.5, cex=0.72, col="darkgreen")
  mtext("plateau",      side=3, at=72,  line=-2.5, cex=0.72, col=col_H)
  mtext("declive",      side=3, at=128, line=-2.5, cex=0.72, col=col_M)
  mtext(sprintf("mediana=%.1fa", res_H$mediana),
        side=3, at=100, line=0.2, cex=0.75, col=col_H)

  # Panel 2: V(t) Mus trayectorias
  par(mar=c(4.5,4.5,3.5,1.5))
  plot(NA, NA, xlim=c(0, R_Musculus$t_cese_target*3), ylim=c(0,1.05),
       xlab="Edad (anos)", ylab="V(t)",
       main="V(t) Mus musculus -- trayectorias (N=30)",
       cex.lab=0.9, cex.axis=0.85, cex.main=0.98)
  abline(h=1.0, lty=2, col="gray80", lwd=1)
  for(sim in 1:30){
    v_i <- res_M$V_tray[sim,]; idx_v <- which(!is.na(v_i))
    if(length(idx_v)>2)
      lines(res_M$t_s[idx_v], v_i[idx_v], col=adjustcolor(col_M,0.35), lwd=0.8)
  }
  V_mean_M <- apply(res_M$V_tray, 2, function(x) mean(x, na.rm=TRUE))
  idx_m <- which(!is.na(V_mean_M) & V_mean_M > 0)
  lines(res_M$t_s[idx_m], V_mean_M[idx_m], col=col_M, lwd=lw)
  mtext(sprintf("mediana=%.2fa  target=%.0fa", res_M$mediana, R_Musculus$t_cese_target),
        side=3, at=R_Musculus$t_cese_target, line=0.2, cex=0.72, col=col_M)

  # Panel 3: D(t) H.sapiens
  par(mar=c(4.5,4.5,3.5,1.5))
  D_H <- res_H$D_d[res_H$D_d>1e-8 & res_H$t_d>R_Hsapiens$t_lat]
  t_H <- res_H$t_d[res_H$D_d>1e-8 & res_H$t_d>R_Hsapiens$t_lat]
  plot(t_H, log10(D_H), type="l", col=col_H, lwd=lw,
       xlim=c(R_Hsapiens$t_lat,160),
       xlab="Edad (anos)", ylab="log10[ D(t) ]",
       main="D(t) H. sapiens -- acumulacion SC",
       cex.lab=0.9, cex.axis=0.85, cex.main=0.98)
  abline(h=log10(Dc), lty=2, col="gray40", lwd=1.5)
  text(50, log10(Dc)+0.08, "Dc=0.60", col="gray40", cex=0.78)

  # Panel 4: D(t) Mus
  par(mar=c(4.5,4.5,3.5,1.5))
  D_M <- res_M$D_d[res_M$D_d>1e-8 & res_M$t_d>R_Musculus$t_lat]
  t_M <- res_M$t_d[res_M$D_d>1e-8 & res_M$t_d>R_Musculus$t_lat]
  if(length(D_M)>2){
    plot(t_M, log10(D_M), type="l", col=col_M, lwd=lw,
         xlim=c(R_Musculus$t_lat, R_Musculus$t_cese_target*2.5),
         xlab="Edad (anos)", ylab="log10[ D(t) ]",
         main="D(t) Mus musculus -- acumulacion SC",
         cex.lab=0.9, cex.axis=0.85, cex.main=0.98)
    abline(h=log10(Dc), lty=2, col="gray40", lwd=1.5)
  }

  # Panel 5: Distribucion t_cese normalizada
  par(mar=c(4.5,4.5,3.5,1.5))
  tv_H <- res_H$t_cese_sto[!is.na(res_H$t_cese_sto)]
  tv_M <- res_M$t_cese_sto[!is.na(res_M$t_cese_sto)]
  d_H  <- density(tv_H/R_Hsapiens$t_cese_target)
  d_M  <- density(tv_M/R_Musculus$t_cese_target)
  ylim_d <- c(0, max(d_H$y, d_M$y)*1.1)
  plot(d_H, col=col_H, lwd=lw, ylim=ylim_d,
       main="Distribucion t_cese (normalizada al target)",
       xlab="t_cese / t_target", ylab="Densidad",
       cex.lab=0.9, cex.axis=0.85, cex.main=0.98)
  lines(d_M, col=col_M, lwd=lw)
  abline(v=1.0, lty=2, col="gray40", lwd=1.5)
  legend("topright",
         c(sprintf("H.sapiens  med=%.1fa", res_H$mediana),
           sprintf("Mus musc.  med=%.2fa", res_M$mediana)),
         col=c(col_H,col_M), lty=1, lwd=lw, bty="n", cex=0.82)

  # Panel 6: Espacio de fase V-D H.sapiens
  par(mar=c(4.5,4.5,3.5,1.5))
  idx_ok <- which(res_H$D_d>1e-8 & res_H$V_d>0.001)
  xl <- range(log10(res_H$D_d[idx_ok]))
  yl <- range(res_H$V_d[idx_ok])
  plot(log10(res_H$D_d[idx_ok]), res_H$V_d[idx_ok],
       type="l", col=col_H, lwd=lw,
       xlim=xl, ylim=c(yl[1]*0.99, yl[2]*1.001),
       xlab="log10[ D(t) ]", ylab="V(t)",
       main="Espacio de fase V-D -- H. sapiens",
       cex.lab=0.9, cex.axis=0.85, cex.main=0.98)
  abline(v=log10(Dc), lty=2, col="gray40", lwd=1.5)
  text(log10(Dc)+0.05, yl[1]*0.991, "Dc", col="gray40", cex=0.8, adj=0)
  mtext(sprintf("Somniosus: med=%.0fa (error %.1f%%) | Hydra: sin colapso en 41a",
        res_S$mediana,
        abs(res_S$mediana - R_Somniosus$t_cese_target)/R_Somniosus$t_cese_target*100),
        side=1, line=3.5, cex=0.70, col="gray30")

  mtext("[EM-1] v6 -- cuatro especies verificadas -- TMCSA jun 2026",
        outer=TRUE, cex=0.95, font=2)
dev.off()
cat(sprintf("\nPNG generado: %s\n", ruta_png))

# ==============================================================================
# R(Olea europaea) -- Olivo mediterraneo
# ESPECIE SINGULAR: (1) veceria (Gamma_vec) + (2) cese externo/estructural
# ------------------------------------------------------------------------------
# COHERENCIA CON LA TRIADA -- Cap5 §5.2i (Araucaria araucana):
#   La Triada establece TRES mecanismos de cese en vegetal lenoso longevo,
#   NINGUNO de los cuales es saturacion interna de SC (X>=Xc):
#     (1) E0 degradado externamente (megasequia/DFA)
#     (2) limite hidraulico geometrico por cavitacion
#     (3) perdida de meristemos apicales por incendio -- colapso de delta
#   El olivo sigue el MISMO patron: el meristemo cambial se renueva, no hay
#   senescencia programada que lleve D a cruzar Dc en la ventana de vida.
#
# CORRECCION v2 (coherencia interna + Triada):
#   v1 calibraba eta para que el colapso INTERNO (X>=Xc) diera t_cese=1000a.
#   Eso CONTRADECIA su propio comentario y la Triada (§5.2i): el cese del
#   olivo NO es por D->Dc. Calibrar eta al colapso interno forzaba un
#   mecanismo que el concepto (Paso 1) dice que no opera.
#
#   v2 corrige siguiendo el patron de Hydra y araucaria:
#     - El proceso interno (D) corre con eta REAL bajo, SIN cruzar Dc en
#       la ventana de vida observada. t_lat alto neutraliza acumulacion neta
#       (meristemo persistente, igual que Hydra neutraliza con t_lat=1000a).
#     - El cese lo da un mecanismo EXTERNO/ESTRUCTURAL (Gamma_ext), no X>=Xc.
#       Gamma_ext es la tasa anual de cese por causa externa: megasequia,
#       incendio de copa, cavitacion, tala. Analogo a gamma_M de hominidos.
#     - Criterio de verificacion: como Hydra, NO es "error vs target de
#       colapso interno" sino "el proceso interno NO colapsa solo; el cese
#       viene de Gamma_ext". El t_cese observado emerge de Gamma_ext.
#
# Esto respeta Paso 1 (concepto de cese = externo, no interno) y la Triada.
#
# t_lat = 5.0a
#   Inicio de actividad reproductiva en olivo (~3-5 anos en condiciones normales).
#   Proxy de t_lat -- pendiente ancla directa de acumulacion neta de SC.
#
# Xc = 1.3678%
#   Formula alometrica EQ-3: 17*(2.5/1000)^0.4206 = 1.368%
#   Umbral mas estricto que H.sapiens (4%) y Somniosus (2.03%).
#   Interpretacion: organismos con longevidad milenaria tienen umbral
#   critico de SC muy bajo -- su sistema de aclaramiento esta calibrado
#   para tolerancia minima al dano.
#
# eps = 0.01448
#   HALLAZGO NUEVO (sesion jun 2026): eps escala con Xc entre especies.
#   eps_O = eps_Mus * (Xc_O / Xc_Mus) = 0.180 * (1.368/17) = 0.01448
#   Razon: con Xc pequeno el sistema es mas fragil al ruido. Si eps no
#   escala con Xc, cualquier fluctuacion cruza el umbral en anos.
#   Verificado: con eps=0.180 (valor de raton) el olivo colapsa en ~4a
#   independientemente de eta. Con eps escalado, la dinamica es correcta.
#   PENDIENTE: verificar si eps escala con Xc tambien en Somniosus.
#
# eta = 7.6996e-05
#   Calibracion inversa biseccion Monte Carlo N=200, error 2.5% (jun 2026).
#   Procedimiento: dado t_cese=1000a, t_lat=5a, Xc=1.368%, Gamma_vec=0.3,
#   se busco el eta que minimiza |mediana_simulada - 1000|.
#   Pendiente: medicion directa de SC en tejido de olivo.
#
# Gamma_vec = 0.3 -- NUEVO PARAMETRO DE CAPA 3
#   Amplitud de la reduccion de beta_ef en ano ON de la veceria.
#   beta_ef(t) = beta * (1 - Gamma_vec * f_ON(t))
#   f_ON(t) = 1 si floor(t) es impar (ano ON), 0 si par (ano OFF)
#   Ancla conceptual: veceria moderada -- reduccion del 30% de la capacidad
#   de remocion de SC en ano de cosecha completa, porque los recursos
#   nutritivos compiten entre reproduccion y mantenimiento membranal.
#   Fuente: Delta Trees 2024 -- veceria parcial: tasas 60-70% ON vs 30-40% OFF.
#   PMC3564680 -- diferencia ON/OFF: 5-30 t/ha.
#   Dominio: Gamma_vec in [0,1). Neutro: 0 (todas las especies anteriores).
#   Otras especies con veceria (manzano, pistache, mango): sus propios valores.
#   OPCION B PENDIENTE: si Gamma_vec=0.3 no verifica bien con datos demograficos
#   directos de olivos milenarios, revisar con funcion coseno (transicion suave).
#
# mu = 2.0 -- mu>1 requerido para P2 (punto fijo estable en dV/dt)
#   Vegetal le oso sin ECM animal. Declive lineal si ocurriera.
#   Misma razon que Hydra pero por arquitectura vegetal, no animal.
#   Pendiente verificacion directa en vegetales le osos.
#
# n = 2.0
#   Pendiente verificacion directa en vegetales. Heredado de mamiferos.
#
# CV_eta = 0.65
#   Referencia H.sapiens -- pendiente datos de variabilidad en olivos.
#   Nota: el IQR simulado [~700-1400a] es coherente con la variabilidad
#   real observada en longevidad de olivos milenarios.
# ==============================================================================

# Funcion de simulacion para vegetal lenoso longevo con:
#   - veceria (Gamma_vec): beta_ef(t) oscila bienalmente
#   - cese externo/estructural (Gamma_ext): tasa anual de cese externo
#     El proceso interno (D) NO colapsa solo en la ventana; el cese
#     viene de Gamma_ext (megasequia/incendio/cavitacion/tala), coherente
#     con araucaria §5.2i. Patron analogo a gamma_M de hominidos.
simular_olivo <- function(R, N_sto=200, seed=42, T_max_det=NULL) {

  kappa_X   <- if(!is.null(R$kappa)) R$kappa else 0.5
  kappa_D   <- kappa_X / R$Xc
  dt        <- 1/365
  dt_sem    <- 7/365
  beta_sem  <- R$beta * dt_sem
  eps_sem   <- R$eps  * sqrt(7)
  eta_det   <- R$eta  * 365 / R$Xc
  gamma_vec <- if(!is.null(R$Gamma_vec)) R$Gamma_vec else 0.0
  gamma_ext <- if(!is.null(R$Gamma_ext)) R$Gamma_ext else 0.0  # cese externo

  stopifnot("P1 violada: a > b requerido" = R$a > R$b)

  CV_eta   <- R$CV_eta
  sigma_ln <- sqrt(log(1 + CV_eta^2))
  mu_ln    <- log(R$eta) - sigma_ln^2/2

  set.seed(seed)

  # Simulacion determinista del proceso interno (sin cese externo)
  # Muestra que D NO cruza Dc por si solo DENTRO DE LA VENTANA DE VIDA REAL.
  # IMPORTANTE: D_max se mide hasta t_cese_target (vida real del organismo),
  # NO sobre T_max_det extendido. Medir D sobre una ventana arbitraria
  # responde la pregunta equivocada: con t_lat alto la acumulacion es lenta
  # pero no nula, asi que dado tiempo infinito siempre cruza Dc. La pregunta
  # correcta es si colapsa DENTRO de su vida real (la que produce Gamma_ext),
  # y a esa edad ningun individuo llega porque Gamma_ext ya lo ceso antes.
  if(is.null(T_max_det)) T_max_det <- R$t_cese_target * 2.5
  N_t <- round(T_max_det / dt)
  V_d <- rep(1.0, N_t); D_d <- rep(0.0, N_t)
  t_d <- (0:(N_t-1)) * dt

  for(i in 2:N_t){
    ta <- t_d[i-1]; Va <- V_d[i-1]; Da <- D_d[i-1]
    if(isTRUE(Va <= 0.001)){ V_d[i]=0; D_d[i]=Da; next }
    f_ON        <- if(floor(ta) %% 2 == 1) 1 else 0
    beta_ef_det <- R$beta * (1.0 - gamma_vec * f_ON)
    fric <- max(0, (1 - Da/Dc)^R$n)
    dD   <- eta_det * max(0, ta * R$W_met - R$t_lat) - (beta_ef_det/R$Xc) * Da/(kappa_D + Da)  # CORREGIDO v7: beta_ef_det/Xc (ver nota v7 en cabecera)
    dV   <- R$a * Va * fric - R$b * Va^R$mu - max(0, Da - Dc)
    V_d[i] <- min(1.10, max(0, Va + dV*dt))
    D_d[i] <- max(Da, Da + dD*dt)  # D monotonicamente no decreciente
  }

  # D_max medido sobre la ventana de vida REAL (hasta t_cese_target),
  # no sobre toda la simulacion extendida
  idx_vida <- which(t_d <= R$t_cese_target)
  D_max_interno <- max(D_d[idx_vida], na.rm=TRUE)

  # Simulacion estocastica: cese por Gamma_ext (externo), no por X>=Xc
  T_sto <- R$t_cese_target * 4
  N_s   <- round(T_sto / dt_sem)
  t_cese_sto <- rep(NA, N_sto)
  causa_ext  <- rep(FALSE, N_sto)

  for(sim in 1:N_sto){
    eta_i <- rlnorm(1, mu_ln, sigma_ln)
    X <- 0.01; V_i <- 1.0; t_fin <- NA; fue_ext <- FALSE
    for(sem in 1:N_s){
      ta          <- (sem - 1) * dt_sem
      xi          <- rnorm(1, 0, 1)
      # Cese externo: tasa anual gamma_ext -> probabilidad por semana
      if(runif(1) < gamma_ext * dt_sem){
        t_fin <- ta; fue_ext <- TRUE; break
      }
      f_ON        <- if(floor(ta) %% 2 == 1) 1 else 0
      beta_ef_sem <- beta_sem * (1.0 - gamma_vec * f_ON)
      dX <- eta_i * max(0, ta * R$W_met - R$t_lat) * 7 -
            beta_ef_sem * X/(kappa_X + X) + eps_sem * xi
      X   <- max(0, X + dX)
      D_i <- X / R$Xc
      fric <- max(0, (1 - D_i/Dc)^R$n)
      dV  <- R$a * V_i * fric - R$b * V_i^R$mu - max(0, D_i - Dc)
      V_i <- min(1.10, max(0, V_i + dV*dt_sem))
      # El cese interno SIGUE disponible (condicion canonica), pero con
      # eta bajo + t_lat alto NO deberia activarse antes que Gamma_ext
      if(isTRUE(X >= R$Xc) || isTRUE(V_i <= 0.01)){ t_fin <- ta; break }
    }
    t_cese_sto[sim] <- t_fin; causa_ext[sim] <- fue_ext
  }

  tv <- t_cese_sto[!is.na(t_cese_sto)]
  pct_ext <- mean(causa_ext[!is.na(t_cese_sto)]) * 100

  cat(sprintf("  %s: mediana=%.0fa  target=%.0fa  error=%.1f%%\n",
      R$nombre, median(tv), R$t_cese_target,
      abs(median(tv)-R$t_cese_target)/R$t_cese_target*100))
  cat(sprintf("    D_max interno (sin cese ext) = %.3f -- %s\n",
      D_max_interno,
      ifelse(D_max_interno < Dc,
             "NO cruza Dc por si solo (cese es externo, correcto)",
             "cruza Dc -- revisar eta/t_lat")))
  cat(sprintf("    cese por causa externa (Gamma_ext) = %.0f%% de los individuos\n",
      pct_ext))

  list(R=R, t_d=t_d, V_d=V_d, D_d=D_d,
       t_cese_sto=t_cese_sto, mediana=median(tv),
       D_max_interno=D_max_interno, pct_ext=pct_ext)
}

R_Olivo <- list(
  nombre        = "Olea europaea",
  t_cese_target = 350.0,   # MEDIANA REAL DEL RAMET (Camarero et al. 2024)
                           # Distribucion C14: mayoria <=300a, clase 400-650a,
                           # outlier 1161a. Ehrlich 2017: duramen rara vez >300a.
                           # 350a = centro de la distribucion; su p90 exponencial
                           # (1163a) reproduce el outlier real de 1161a.
                           # CORRECCION: el target anterior (1000a) era del GENET
                           # o del outlier, no la mediana del RAMET (tronco vivo).

  # Capa 1 -- Viabilidad
  a  = 0.0988,
  b  = 0.0960,
  mu = 2.0,    # mu>1 requerido para punto fijo estable (P2)
  n  = 2.0,    # pendiente verificacion en vegetales

  # Capa 2 -- Dano (SC)
  eta   = 7.6996e-05,  # eta REAL bajo -- NO calibrado al colapso interno
  t_lat = 150.0,       # meristemo cambial persistente: acumulacion neta
                       # retrasada. El proceso interno determinista NO cruza
                       # Dc dentro de la vida real del ramet -- cese externo.
                       # Respaldo: meristemo persistente (mismo mecanismo que
                       # Hydra t_lat=1000a). Concepto, no ajuste numerico.
  Xc    = 1.3678,      # EQ-3 con t_activo del genet (referencia estructural)
  beta  = 54.75,
  eps   = 0.01448,     # escalado por Xc: 0.180*(1.3678/17)

  # Capa 3 -- Modificadores
  W_met      = 1.0,
  delta_beta = 0.0,
  rho        = 0.0,
  Gamma_vec  = 0.3,    # veceria moderada -- ancla Delta Trees 2024
  # Gamma_ext -- cese externo/estructural (coherente araucaria §5.2i)
  # CALIBRADO desde dato demografico real: distribucion de edades del ramet
  # (Camarero 2024). Gamma_ext = ln(2)/mediana = ln(2)/350 = 0.00198/a.
  # La exponencial resultante reproduce: mediana 350a, clase intermedia
  # 400-650a, y p90=1163a == outlier real 1161a. El cese es por megasequia/
  # incendio/cavitacion/tala. Analogo a gamma_M de hominidos.
  Gamma_ext  = 0.00198,  # ln(2)/350 -- dato demografico real Camarero 2024.
  # HALLAZGO jun2026: el valor 0.00198 es CORRECTO. El error ~10% que aparece
  # con seed=42 es VARIANZA DE MUESTREO, no un parametro mal calibrado: la
  # distribucion del olivo es de COLA MUY LARGA (individuos hasta 1800+a), y la
  # mediana muestral con N finito es inestable entre semillas (seed=42->323a err
  # 10.4%; seed=100->358a err 2.3%). El PROMEDIO de medianas cae cerca de 350.
  # NO ajustar Gamma_ext. Reportar error como promedio de semillas o subir N.
  # VERIFICADO: error por semilla varia de 0.5% a 8.2% (seed7=352a/0.5%,
  # seed100=341a/2.4%, seed42=321a/8.2%). ERROR REAL PROMEDIADO ~3-4%, NO 10.4%.
  # El 10.4% es artefacto de la semilla 42 (la mas desfavorable). El olivo esta
  # CORRECTO; su error real es ~3-4%. (Ver metodo Parte XII: varianza muestreo.)
  # NOTA: el metodo inverso arrojo 0.00175 como valor que daba 1.1% CON seed=42,
  # pero ese valor DA 20.7% con seed=100 -> INESTABLE, DESCARTADO. NO reintroducir.
  # No es una "tendencia con mejores datos": es un parche a una semilla. El valor
  # correcto es 0.00198, que con N suficiente / semillas promediadas tiende a ~3-4%.

  # CV_eta: PENDIENTE PASO 5 -- NO heredar de H.sapiens sin verificar.
  # Pregunta Paso 5 abierta: ¿que condicion biologica determina la
  # variabilidad individual de eta en olivos? (clonal vs semilla, sitio).
  CV_eta=0.65  # PASO5: monumental cultivares mixtos (Muzzalupo,Belaj) -- coincide-diversidad SELLADO  # PROVISIONAL -- pendiente Paso 5 para olivo
)

cat("\n=== Olea europaea (veceria + cese externo Gamma_ext) ===\n")
cat("Target = mediana real del ramet (Camarero 2024), no del genet\n")
res_O <- simular_olivo(R_Olivo, N_sto=200, T_max_det=1400)

# ==============================================================================
# R(Heterocephalus glaber) -- Rata topo desnuda
# CORRECCION MAYOR v2 (jun 2026): modo de cese corregido.
#   v1: K basico con beta alto (beta_apoptotico), error 10.2% vs target 30a.
#   v2: negligible senescence por TASA CONSTANTE. Error 0.7% vs half-life 19a.
#
#   DATO QUE OBLIGA EL CAMBIO (NIVEL 1):
#   Ruby et al. 2018 (eLife) + GeroScience 2024 (N>3000):
#   H.glaber NO sigue Gompertz. Hazard ~0.0365/a a CUALQUIER edad.
#   Consistencia: 1/10000/dia = 0.0365/a -> half-life=19.0a (exacto).
#
#   El proceso interno (D) corre con beta=400 (SCD) y NO colapsa
#   (D_max=0.209<Dc). El cese viene de Gamma_ext=0.0365/a.
#   Mismo patron que olivo/Fitzroya pero Gamma_ext no es ambiental:
#   es mortalidad intrinseca de tasa constante de especie no-senescente.
#
# REGLA DERIVADA: ante error alto, verificar el MODO antes de ajustar
# parametros. El dato que decide el modo es la curva de mortalidad por edad.
#
# PARAMETROS v2:
#   Gamma_ext=0.0365/a: NIVEL 1 -- Ruby 2018 (1/10000/dia, half-life 19a).
#   t_lat=13a:          NIVEL 1 -- Emmrich 2021.
#   beta=400:           NIVEL 2 -- beta_ef SCD, Kawaguchi 2020.
#   a=0.0692,b=0.0672:  NIVEL 2 -- Buffenstein & Yahav 1991.
#   eta=0.050:          NIVEL 3 -- estimado alometrico.
#   Xc=6.46:            NIVEL 2 -- alometrico 17*(2.5/25)^0.42.
#   CV_eta=0.20:        NIVEL 2 -- colonia laboratorio.
# ------------------------------------------------------------------------------
R_Heterocephalus <- list(
  nombre        = "Heterocephalus glaber",
  t_cese_target = 19.0,   # half-life (mediana cese), Ruby 2018
  a=0.0692, b=0.0672, mu=3.0, n=2.0,
  eta   = 0.050,
  t_lat = 13.0,
  Xc    = 6.46,
  beta  = 400.0,
  eps=0.0684,  # AUDIT: escalado 0.180*(Xc/17), antes 0.180 heredado
  W_met=1.0, delta_beta=0.0, rho=0.0, CV_eta=0.20,
  Gamma_ext = 0.0365
)

cat("\n=== Heterocephalus glaber (v2 -- negligible senescence, tasa constante) ===\n")
cat("Modo corregido: Gamma_ext=0.0365/a (Ruby 2018). NO senescencia Gompertziana.\n")
res_HG <- simular_olivo(R_Heterocephalus, N_sto=400, T_max_det=40)


# ==============================================================================
# R(Fitzroya cupressoides) -- Alerce patagonico
# Modo: vegetal leñoso longevo con cese externo dominante
# Mismo patron que Olea europaea (Triadad S5.2i)
# ------------------------------------------------------------------------------
# ORIGEN DE CADA PARAMETRO:
#
# Gamma_ext = 0.00170/a
#   Fuente: estimador Nelson-Aalen con censura a la derecha sobre CHIL015
#   (Villalba 1996, ITRDB/NOAA, Tiuchue Isla de Chiloe, n=33 arboles).
#   12 eventos (arboles que terminaron antes de 1987) + 21 censurados.
#   Mediana ramet = ln(2)/0.00170 = 408a.
#   Coherencia: p99.9 exponencial = 4063a vs maximo verificado 3622a (Lara 1993).
#   Nivel NIVEL 1 -- datos demograficos reales propios.
#
# t_cese_target = 408a
#   Mediana ramet derivada de Gamma_ext. NO es el maximo de 3622a (ese es p99.9).
#
# Xc = 0.7960%
#   EQ-3: 17*(2.5/3622)^0.4206 usando t_activo del genet como referencia
#   estructural. Xc mas estricto que olivo (1.37%) -- genet mas longevo.
#
# eps = 0.008429
#   Escalado por Xc: 0.180*(0.7960/17).
#
# t_lat = 175a
#   Meristemo cambial persistente -- mismo concepto que olivo (t_lat=150a).
#   Escalado por longevidad: 150*(408/350) = 175a.
#   D_max=0.111 < Dc=0.60 dentro de la vida real del ramet (verificado).
#
# eta = 7.6996e-05
#   Mismo orden que Olea europaea -- vegetal leñoso longevo.
#   Con Gamma_ext dominante, eta tiene poco efecto en rango bajo.
#   Pendiente medicion directa de SC en tejido de alerce.
#
# mu = 2.0
#   Vegetal leñoso -- mu>1 requerido para P2. Igual Olea europaea.
#
# Gamma_vec = 0.0
#   Conifera -- sin veceria bienal verificada.
#
# Error: 6.0% vs target=408a. Aceptable dado n=33 en CHIL015.
# D_max=0.111 < Dc=0.60 -- cese externo correcto (86% de los ceses).
# ==============================================================================
R_Fitzroya <- list(
  nombre        = "Fitzroya cupressoides",
  t_cese_target = 408.0,
  a=0.0988, b=0.0960, mu=2.0, n=2.0,
  eta   = 7.6996e-05,
  t_lat = 175.0,
  Xc    = 0.7960,
  beta  = 54.75,
  eps   = 0.008429,
  W_met=1.0, delta_beta=0.0, rho=0.0,
  Gamma_vec = 0.0,
  Gamma_ext = 0.00170,
  CV_eta=0.65  # PASO5: salvaje bosque nativo diverso -- coincide-diversidad
)

# ==============================================================================
# R(Turritopsis dohrnii) -- Medusa inmortal
# Modo: unico caso verificado de rho>0 en [EM-1]
# ------------------------------------------------------------------------------
# MECANISMO -- leer antes de modificar:
#   La transdegeneracion es un evento DISCRETO (no continuo).
#   Se activa cuando D > umbral_stress con tasa p_stress por unidad de tiempo.
#   Al activarse: D->0, V->1 (reseteo completo, rho=1.0).
#   El individuo reinicia el ciclo como polipo juvenil.
#   Fuente: Piraino 1996; Miglietta et al. 2009; Schmich et al. 2007.
#
# DOS REGIMENES:
#   Laboratorio (Gamma_ext=0): sin colapso interno en 4a = 500/500 (100%)
#   Natural (Gamma_ext=2.0/a): 100% depredacion, mediana=0.31a (~4 meses)
#
# ORIGEN DE CADA PARAMETRO:
#
# rho = 1.0
#   Reseteo completo: D->0, V->1.
#   Unico caso verificado en [EM-1]. Para todas las demas especies: rho=0.
#
# p_stress = 20.0/a
#   Tasa de transdegeneracion bajo estres.
#   Calibrado: ~0.65 resets/individuo en 4a de laboratorio.
#   Sin medicion directa de frecuencia. Pendiente.
#
# umbral_stress = 0.10
#   D/Dc donde se activa la posibilidad de transdegeneracion.
#   Sin medicion directa. Pendiente.
#
# Gamma_ext = 2.0/a
#   Depredacion natural. Mediana silvestre ~0.31a (~4 meses).
#   "la mayoria sucumbe a la depredacion sin revertir" (Animalia.bio).
#   Sin medicion directa de tasa de depredacion. Estimacion cualitativa.
#
# eta = 1.0e-02
#   Calibrado para que rho sea necesario (D sube hasta umbral) y suficiente
#   (0 colapsos en 4a con p_stress=20).
#   Pendiente medicion directa de SC en tejido de T. dohrnii.
#
# Xc = 24.9932%
#   EQ-3: 17*(2.5/1.0)^0.4206 con t_activo=1a (vida media medusa).
#
# mu = 1.0
#   Cnidario sin ECM densa -- igual que Hydra. Justificado.
#
# Criterio de verificacion: sin colapso INTERNO en 4a (laboratorio).
# Analogo a Hydra: t_obs_verificado=4.0a.
# ==============================================================================
R_Turritopsis <- list(
  nombre        = "Turritopsis dohrnii",
  t_cese_target = 4.0,
  a=0.0988, b=0.0960, mu=1.0, n=1.5,
  eta   = 1.0e-02,
  t_lat = 0.1,
  Xc    = 24.9932,
  beta  = 54.75,
  eps   = 0.26463,
  W_met=1.0, delta_beta=0.0,
  rho          = 1.0,
  p_stress     = 20.0,
  umbral_stress= 0.10,
  Gamma_ext    = 2.0,
  verificar_sin_colapso = TRUE,
  t_obs_verificado      = 4.0,
  CV_eta=0.1  # PASO5: clonal laboratorio (1500 medusas 1 clon, DNA Research 2023) + rejuvenece como Hydra -> CV bajo NIVEL3
)

# Funcion de simulacion para T. dohrnii (rho como evento discreto)
simular_turritopsis <- function(R, N_sto=500, seed=42, con_gamma_ext=FALSE) {
  kappa_X  <- 0.5
  dt_sem   <- 7/365
  beta_sem <- R$beta * dt_sem
  eps_sem  <- R$eps  * sqrt(7)
  sigma_ln <- sqrt(log(1 + R$CV_eta^2))
  mu_ln    <- log(R$eta) - sigma_ln^2/2
  gamma_ext <- if(con_gamma_ext) R$Gamma_ext else 0.0
  set.seed(seed)
  T_sto <- R$t_obs_verificado * 1.2
  N_s   <- round(T_sto / dt_sem)
  t_cese_sto <- rep(NA, N_sto); n_resets <- rep(0, N_sto)
  causa_ext  <- rep(FALSE, N_sto)
  for(sim in 1:N_sto){
    eta_i <- rlnorm(1, mu_ln, sigma_ln)
    X <- 0.01; V_i <- 1.0; t_fin <- NA; fue_ext <- FALSE; nr <- 0
    for(sem in 1:N_s){
      ta <- (sem-1)*dt_sem
      if(con_gamma_ext && runif(1) < gamma_ext*dt_sem){
        t_fin <- ta; fue_ext <- TRUE; break
      }
      D_i <- X/R$Xc
      if(D_i > R$umbral_stress && runif(1) < R$p_stress*dt_sem){
        X <- 0.001; V_i <- 1.0; nr <- nr+1; next
      }
      xi  <- rnorm(1,0,1)
      dX  <- eta_i*max(0,ta-R$t_lat)*7 - beta_sem*X/(kappa_X+X) + eps_sem*xi
      X   <- max(0, X+dX)
      D_i <- X/R$Xc
      fric <- max(0,(1-min(D_i/Dc,1.0))^R$n)
      dV  <- R$a*V_i*fric - R$b*V_i^R$mu - max(0,D_i-Dc)
      V_i <- min(1.10, max(0,V_i+dV*dt_sem))
      if(isTRUE(X>=R$Xc)||isTRUE(V_i<=0.01)){ t_fin<-ta; break }
    }
    t_cese_sto[sim]<-t_fin; n_resets[sim]<-nr; causa_ext[sim]<-fue_ext
  }
  t_obs <- R$t_obs_verificado
  n_ok  <- sum(t_cese_sto >= t_obs | is.na(t_cese_sto), na.rm=TRUE)
  regimen <- if(con_gamma_ext) "natural" else "laboratorio"
  cat(sprintf("  %s [%s]: sin colapso interno en %.0fa = %d/%d (%.0f%%)  resets/ind=%.2f\n",
      R$nombre, regimen, t_obs, n_ok, N_sto, 100*n_ok/N_sto, mean(n_resets)))
  list(R=R, t_cese_sto=t_cese_sto, n_resets=n_resets, n_ok=n_ok)
}

cat("\n=== Fitzroya cupressoides (cese externo -- Nelson-Aalen CHIL015) ===\n")
res_FZ <- simular_olivo(R_Fitzroya, N_sto=300, T_max_det=1200)

cat("\n=== Turritopsis dohrnii (rho=1.0 -- transdegeneracion) ===\n")
res_TD_lab <- simular_turritopsis(R_Turritopsis, N_sto=500, con_gamma_ext=FALSE)
res_TD_nat <- simular_turritopsis(R_Turritopsis, N_sto=500, con_gamma_ext=TRUE)

# ==============================================================================
# R(Escherichia coli) -- bacteria gram-negativa
# Modo: asimetria de polo -- primer caso de bifurcacion del vector vida en [EM-1]
# ------------------------------------------------------------------------------
# CONCEPTO FUNDAMENTAL:
#   En E. coli el "cese" del individuo no es muerte -- es bifurcacion.
#   El proceso vida se ramifica en dos nodos: polo viejo (envejece) y
#   polo nuevo (rejuvenecido). Se simula la linea del polo viejo.
#
# UNIDADES DE TIEMPO: GENERACIONES (no anos ni minutos).
#   Razon: Stewart 2005 mide en generaciones. EQ-3 con t_activo=100gen
#   da Xc=3.6% comparable a H.sapiens (4%) -- invariante biologica.
#
# MECANISMO DE DIVISION (correccion conceptual):
#   1. Se integra dV/dt y dD/dt durante la generacion
#   2. Division ANTES del cese: X_post = alpha * X
#   3. Cese evaluado sobre D_post (no D_pre) -- respeta P10
#
# VERIFICACION CIENTIIFICA:
#   Stewart et al. 2005 (PLoS Biol, n=35049 celulas):
#   cells inheriting old poles had reduced growth rate and increased death.
#   Lindner et al. 2008 (PNAS): ratio tasas bajo estres = 0.83 -> alpha=0.70
#
# ORIGEN DE CADA PARAMETRO:
# t_cese_target=100 gen: Kirkwood 2005. Vida maxima de la celula madre.
# t_lat=33 gen: 33% de t_vida -- mismo ratio que H.sapiens (40/122=0.33).
# Xc=3.6027%: EQ-3: 17*(2.5/100)^0.4206. Casi identico a H.sap (4.0%) --
#   invariante biologica del umbral critico de dano entre organismos.
# eps=0.03815: escalado por Xc: 0.180*(3.6027/17).
# eta=0.04641 gen^-1: calibracion inversa MC N=500, error 2.0% (jun 2026).
#   Logica corregida: division antes del cese.
# alpha=0.70: fraccion de X que hereda la madre (polo viejo).
#   Verificado: Stewart 2005, Lindner 2008.
# beta=54.75 gen^-1: referencia H.sapiens en unidades de generacion.
#   Pendiente calibracion directa para E. coli.
# mu=2.0: sin ECM densa -- mu>1 requerido para P2.
# HALLAZGO: alpha es el primer parametro de BIFURCACION en [EM-1].
#   No es un parametro de cese sino de transmision del vector vida.
# ==============================================================================
R_Ecoli <- list(
  nombre        = "Escherichia coli (linea polo viejo)",
  t_cese_target = 100.0,
  a=0.0988, b=0.0960, mu=2.0, n=2.0,
  eta   = 0.04641,
  t_lat = 33.0,
  Xc    = 3.6027,
  beta  = 54.75,
  eps   = 0.03815,
  W_met=1.0, delta_beta=0.0, rho=0.0,
  alpha = 0.70,
  CV_eta=0.65  # PASO5: clonal PERO heterogeneidad polo viejo/nuevo por division asimetrica (Stewart 2005, Ackermann) -- coincide por mecanismo, NO diversidad genetica
)

# Funcion de simulacion con asimetria de polo
simular_ecoli <- function(R, N_sto=500, seed=42) {
  kappa_X  <- 0.5
  sigma_ln <- sqrt(log(1 + R$CV_eta^2))
  mu_ln    <- log(R$eta) - sigma_ln^2/2
  alpha    <- R$alpha
  set.seed(seed)
  T_max <- R$t_cese_target * 5
  t_cese_sto <- rep(NA, N_sto)
  for(sim in 1:N_sto){
    eta_i <- rlnorm(1, mu_ln, sigma_ln)
    X <- 0.01; V_i <- 1.0; t_fin <- NA
    for(gen in 1:T_max){
      xi  <- rnorm(1,0,1)
      dX  <- eta_i*max(0,gen-R$t_lat) - R$beta*X/(kappa_X+X) + R$eps*xi
      X   <- max(0, X+dX)
      D_i <- X/R$Xc
      fric <- max(0,(1-min(D_i/Dc,1.0))^R$n)
      dV  <- R$a*V_i*fric - R$b*V_i^R$mu - max(0,D_i-Dc)
      V_i <- min(1.10,max(0,V_i+dV))
      X   <- alpha * X  # division primero: D_post = alpha*D_pre
      if(isTRUE(X >= R$Xc) || isTRUE(V_i <= 0.01)){ t_fin<-gen; break }
    }
    t_cese_sto[sim] <- t_fin
  }
  tv <- t_cese_sto[!is.na(t_cese_sto)]
  cat(sprintf("  %s: muertes=%d/%d  mediana=%.0f gen  target=%.0f gen  error=%.1f%%\n",
      R$nombre, length(tv), N_sto, median(tv), R$t_cese_target,
      abs(median(tv)-R$t_cese_target)/R$t_cese_target*100))
  list(R=R, t_cese_sto=t_cese_sto, mediana=median(tv))
}

cat("\n=== Escherichia coli (asimetria de polo alpha=0.70) ===\n")
cat("Unidades de tiempo: GENERACIONES\n")
res_EC <- simular_ecoli(R_Ecoli, N_sto=500)

# ==============================================================================
# R(Pan troglodytes) -- Chimpance comun
# Especie 10 -- cierra el canonico 10A
# Modo 1 K basico -- comparacion directa con H.sapiens (98.8% ADN compartido)
# ------------------------------------------------------------------------------
# ORIGEN DE CADA PARAMETRO:
# t_cese_target=42.4a: MLE hembras, Havercamp 2023 (Am J Primatol, n=2349).
#   Se usa MLE hembras: mayor muestra, mas comparable con H.sapiens.
# t_lat=21.74a: formula triada S3A.17: 40*(42.4/78).
#   Mismo mecanismo que H.sapiens: involucion timica.
# Xc=5.1684%: EQ-3: 17*(2.5/42.4)^0.4206.
# eps=0.05472: 0.180*(5.1684/17).
# eta=8.125e-3: calibracion inversa MC N=500, error 0.1%.
#   eta_PT/eta_H = 5.8x para longevidad 2.9x menor.
#   Patron emergente entre mamiferos K basico:
#   H.sap(1.4e-3,122a), Pan(8.1e-3,42a), H.glaber(0.050,19a), Mus(0.084,4a)
#   R2=0.865, exponente -1.15 (eta ~ t_vida^-1.15).
# ==============================================================================
R_Pan <- list(
  nombre        = "Pan troglodytes",
  t_cese_target = 42.4,
  a=0.0988, b=0.0960, mu=3.0, n=2.0,
  eta   = 8.125e-3,
  t_lat = 21.74,
  Xc    = 5.1684,
  beta  = 54.75,
  eps   = 0.05472,
  W_met=1.0, delta_beta=0.0, rho=0.0,
  Gamma_vec=0.0, Gamma_ext=0.0,
  CV_eta=0.65  # PASO5: salvaje/cautiverio diverso -- coincide-diversidad
)

# ==============================================================================
# R(Chelonoidis sp.) -- Tortuga gigante de Galapagos
# Modo: negligible senescence por tasa constante
# PRIMERA APLICACION CORRECTA DE LA REGLA DEL MODO DESDE EL PASO 1:
#   modo verificado con dato demografico ANTES de calibrar.
# ------------------------------------------------------------------------------
# DATO QUE DECIDE EL MODO (NIVEL 1):
#   Da Silva et al. 2022 (Science): negligible senescence verificada
#   en 160 tortugas de Galapagos (ZIMS). Supervivencia adulta ~98%/anio.
#
# ORIGEN DE CADA PARAMETRO:
# Gamma_ext=0.0202/a: -ln(0.98). Supervivencia adulta 98%/anio (da Silva 2022).
#   Coherencia: 1/Gamma_ext=49.5a vida adulta media + t_lat=20a = 69.5a total
#   coincide con esperanza de vida adulta media 69.2a (da Silva 2022).
# t_lat=20a: madurez sexual, da Silva 2022.
# Xc=4.2060%: EQ-3 con t_activo=69.2a.
#   INVARIANTE BIOLOGICA: H.sap(4.0%), E.coli(3.6%), Chelonoidis(4.21%)
#   -- tres reinos distintos, misma escala de dano critico.
# W_met=0.80: ectotermo ~80% actividad metabolica. NIVEL 2 estimado.
# eta=5e-4: NIVEL 3. Gamma_ext domina -- eta sin efecto en rango bajo.
# D_max=0.013 << Dc=0.60 -- proceso interno no colapsa. Correcto.
# ==============================================================================
R_Chelonoidis <- list(
  nombre        = "Chelonoidis sp.",
  t_cese_target = 34.3,   # mediana adulta = ln(2)/Gamma_ext
  a=0.0988, b=0.0960, mu=2.0, n=2.0,
  eta   = 5e-4,
  t_lat = 20.0,
  Xc    = 4.2060,
  beta  = 54.75,
  eps   = 0.04453,
  W_met=1.0,  # PASO 5: sin condicion que lo modifique. No hiberna.
              # W_met<1 empeora ajuste y contaria doble con Gamma_ext+eta.
              # Verificado: W_met=1.0->err 1.1%; W_met=0.80->err 2.2%.
  delta_beta=0.0, rho=0.0,
  Gamma_vec=0.0,
  Gamma_ext=0.0202,
  CV_eta=0.65  # PASO5: zoo multi-origen ZIMS (da Silva 2022) -- coincide-diversidad
)

cat("\n=== Pan troglodytes (K basico -- especie 10 del canonico 10A) ===\n")
res_PT <- simular_especie(R_Pan, N_sto=500)

cat("\n=== Chelonoidis sp. (negligible senescence -- regla modo correcta) ===\n")
res_CH <- simular_olivo(R_Chelonoidis, N_sto=500, T_max_det=180)

# ==============================================================================
# GRAFICO DE PARIDAD de conjunto (10A) -- target vs mediana, log-log
# Alimentado por las corridas REALES de arriba (res_*). NO numeros a mano.
# Aporte: muestra que [EM-1] predice longevidad sobre ordenes de magnitud
# (Mus 4a -> Fitzroya 408a) con la misma ecuacion, dos modos (color).
# Hydra y Turritopsis se excluyen (negligible/sin colapso, sin t_cese finito).
# ==============================================================================
fecha_par <- format(Sys.time(), "%Y%m%d_%H%M")

# Recolectar resultados reales con t_cese finito
res_list <- list(
  list(nombre="H.sapiens",  res=res_H,  modo="sen"),
  list(nombre="Mus",        res=res_M,  modo="sen"),
  list(nombre="Somniosus",  res=res_S,  modo="sen"),
  list(nombre="Olea",       res=res_O,  modo="ext"),
  list(nombre="H.glaber",   res=res_HG, modo="ext"),
  list(nombre="Fitzroya",   res=res_FZ, modo="ext"),
  list(nombre="E.coli",     res=res_EC, modo="sen"),
  list(nombre="Pan",        res=res_PT, modo="sen"),
  list(nombre="Chelonoidis",res=res_CH, modo="ext")
)

targets <- sapply(res_list, function(x) x$res$R$t_cese_target)
medianas <- sapply(res_list, function(x){
  tv <- x$res$t_cese_sto[!is.na(x$res$t_cese_sto)]
  median(tv)
})
nombres <- sapply(res_list, function(x) x$nombre)
modos   <- sapply(res_list, function(x) x$modo)
errs    <- abs(medianas - targets)/targets*100

ruta_par_10A <- paste0("TMCSA_EM1_10A_paridad_", fecha_par, ".png")
png(ruta_par_10A, width=2400, height=2400, res=300)
par(mar=c(5,5,4,2))
cols <- ifelse(modos=="sen","#2166AC","#B2182B")
rango <- range(c(targets, medianas))
plot(targets, medianas, log="xy", pch=19, col=cols, cex=1.8,
     xlim=rango, ylim=rango,
     xlab="Target real (anos, escala log)",
     ylab="Mediana simulada [EM-1] (anos, log)",
     main="10A: paridad [EM-1] sobre ordenes de magnitud")
abline(0,1,lty=2,col="gray40",lwd=1.5)
lines(rango,rango*1.1,lty=3,col="gray70"); lines(rango,rango*0.9,lty=3,col="gray70")
text(targets, medianas, labels=nombres, pos=4, cex=0.62, col="gray20")
legend("topleft",
       legend=c("Senescencia","Cese externo","Prediccion y=x","Banda +-10%"),
       col=c("#2166AC","#B2182B","gray40","gray70"),
       pch=c(19,19,NA,NA), lty=c(NA,NA,2,3), lwd=c(NA,NA,1.5,1), bty="n", cex=0.7)
mtext(sprintf("9 especies con t_cese finito (Hydra/Turritopsis negligible). Error medio=%.1f%%",
      mean(errs)), side=1, line=3.8, cex=0.58, font=3)
dev.off()
cat(sprintf("\nPNG paridad 10A: %s\n", ruta_par_10A))
cat(sprintf("Ruta: %s/%s\n", getwd(), ruta_par_10A))
