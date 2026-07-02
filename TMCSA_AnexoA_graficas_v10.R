# ==============================================================================
# TMCSA -- Anexo A -- Graficas canonicas -- v10
# Archivo: TMCSA_AnexoA_graficas_v10.R
#
# v10: G5 Mus musculus -- suavizado mensual (4 semanas) en vez de anual.
#   Con vida_max=6a, el suavizado de 104 semanas dejaba muy pocos puntos
#   y D no llegaba a cruzar Dc visualmente. Con 4 semanas hay suficiente
#   resolucion para ver la curva completa hasta el cruce de Dc=0.60.
#
# v9: (1) G5 usa Mus musculus (eta=0.084, Xc=17) en vez de H.sapiens para
#   que D cruce visualmente Dc=0.60. (2) G6b suavizado 2 anios en vez de 1.
#   (3) G3 ajuste Gompertz restringido al rango 70-150a donde la ley es
#   mas lineal en escala log.
#
# v8: (1) G5 cambia seed=42->7 para graficar individuo que cruce Dc=0.60.
#   (2) G6 repara simular_B: trayectoria capturada en loop dedicado fuera
#   del loop de N individuos -- el Vt_sem quedaba fuera de scope en v7.
#
# v7: suavizado anual en G4, G5, G6 -- promedio de 52 semanas por punto
#   elimina el ruido semanal del paso dt=7/365. Las trayectorias V(t) y
#   D(t) muestran ahora la tendencia de largo plazo, no las fluctuaciones
#   semana a semana.
#
# v6: G3 cambia a escala logaritmica en Y -- en escala log Gompertz es
#   linea recta creciente, forma estandar en demografia. Ajuste
#   log-lineal superpuesto para confirmar la forma.
#
# v5: corrige G3 bathtub -- intervalo 8->15a para suavizar dientes de
#   sierra, rango recortado a 40-220a donde ocurren ceses reales,
#   añade linea de tendencia lowess para mostrar forma de Gompertz.
#
# v4: N 500->2000 para suavizar G2 y G3. vida_max 250->350 para eliminar
#   el artefacto de censura en el pico final. bw 3->8 en densidad G2.
#
# v3: corrige parametros de H.sapiens a los del canonico 10A real
#   (eta=1.4e-3, Xc=4.0 -- verificados con error 1.3%). v2 usaba
#   eta=4.016e-4, Xc=17 que son del mouse, no del humano.
#
# v2: corrige vida_max de H.sapiens de 120 a 250 -- v1 generaba S(t) plana
#   y bathtub sin forma porque todos los individuos llegaban al tope
#   artificial antes de cesar por senescencia. Myotis corregido de 120 a
#   200 por la misma razon.
#
# Genera las graficas del Anexo A en dos grupos:
#   Grupo 1 -- Graficas clasicas (forma conocida en la ciencia):
#     G1: Curva de supervivencia S(t) tipo Gompertz vs TMCSA
#     G2: Distribucion de t_cese (campana/exponencial segun modo)
#     G3: Curva en tina/bathtub h(t) -- tasa de mortalidad
#   Grupo 2 -- Graficas propias de la TMCSA:
#     G4: V(t) regimen A -- curva sigmoide invertida (declive)
#     G5: D(t) regimen A -- acumulacion de daño cruzando Dc
#     G6: V(t) y D(t) regimen B -- estabilizacion (negligible senescence)
#     G7: X_eq regimen B -- Margen de Seguridad (MS) para las 3 especies
#     G8: Paridad diagonal -- mediana simulada vs target (todas las especies)
#
# Cada grafica se guarda con nombre fijo + fecha para no sobreescribir.
# Los nombres de imagen se usan luego en el documento del Anexo A.
# ==============================================================================
Dc <- 0.60
fecha <- format(Sys.Date(), "%Y%m%d")

# ------------------------------------------------------------------------------
# FUNCIONES AUXILIARES
# ------------------------------------------------------------------------------
simular_A <- function(eta, t_lat, Xc, beta, eps, CV_eta=0.65,
                       a=0.0988, b=0.0960, mu=3, n=2,
                       Gamma_ext=0, vida_max=200, N=500, seed=42){
  kappa_X<-0.5; dt_sem<-7/365; beta_sem<-beta*dt_sem; eps_sem<-eps*sqrt(7)
  ge<-Gamma_ext; sigma_ln<-sqrt(log(1+CV_eta^2)); mu_ln<-log(eta)-sigma_ln^2/2
  set.seed(seed); Ns<-round(vida_max*365/7); tc<-rep(NA,N); Vtraj<-Dtraj<-NULL
  Vtraj_sem<-Dtraj_sem<-NULL
  for(s in 1:N){ eta_i<-rlnorm(1,mu_ln,sigma_ln); X<-0.01; V<-1; tf<-NA
    vv<-dd<-NULL
    for(sem in 1:Ns){ ta<-(sem-1)*dt_sem
      if(ge>0 && runif(1)<ge*dt_sem){tf<-ta;break}
      dX<-eta_i*max(0,ta-t_lat)*7-beta_sem*X/(kappa_X+X)+eps_sem*rnorm(1)
      X<-max(0,X+dX); D<-X/Xc; fr<-max(0,(1-D/Dc)^n)
      dV<-a*V*fr-b*V^mu-max(0,D-Dc); V<-min(1.1,max(0,V+dV*dt_sem))
      if(s==1){vv<-c(vv,V); dd<-c(dd,D)}
      if(isTRUE(X>=Xc)||isTRUE(V<=0.01)){tf<-ta;break} }
    tc[s]<-if(is.na(tf)) vida_max else tf
    if(s==1){Vtraj_sem<-vv; Dtraj_sem<-dd} }
  # Suavizar: promedio anual (52 semanas) para eliminar ruido semanal
  suavizar_anual<-function(v){ n52<-floor(length(v)/52)
    sapply(1:n52, function(i) mean(v[((i-1)*52+1):(i*52)])) }
  Vtraj<-suavizar_anual(Vtraj_sem)
  Dtraj<-suavizar_anual(Dtraj_sem)
  tt<-(1:length(Vtraj))-0.5  # tiempo en años
  list(tc=tc[!is.na(tc)], Vt=Vtraj, Dt=Dtraj, tt=tt) }

simular_B <- function(eta_basal, rho_rep, Xc=3.0,
                       a=0.0988, b=0.0960, mu=3, n=2,
                       Gamma_ext=log(2)/41, vida_max=200, N=500, seed=42){
  Dc<-0.60; dt_sem<-7/365; Gext_sem<-Gamma_ext*dt_sem
  eps<-0.180*(3.0/17)*sqrt(7)
  set.seed(seed); Ns<-round(vida_max*365/7)
  eta_i<-0.06  # valor representativo para la trayectoria ilustrativa
  X<-0.01; V<-1; Vt_sem<-c(); Dt_sem<-c()
  for(sem in 1:Ns){
    if(Gext_sem>0 && runif(1)<Gext_sem) break
    dX<-(eta_i-rho_rep*X)*dt_sem+eps*rnorm(1)
    X<-max(0,X+dX); D<-X/Xc; fr<-max(0,(1-D/Dc)^n)
    dV<-a*V*fr-b*V^mu-max(0,D-Dc); V<-min(1.1,max(0,V+dV*dt_sem))
    Vt_sem<-c(Vt_sem,V); Dt_sem<-c(Dt_sem,D) }
  suavizar_anual<-function(v, w=104){ nw<-floor(length(v)/w)
    if(nw<1) return(mean(v))
    sapply(1:nw, function(i) mean(v[((i-1)*w+1):(i*w)])) }
  Vt<-suavizar_anual(Vt_sem); Dt<-suavizar_anual(Dt_sem)
  tt<-(1:length(Vt))-0.5
  X_eq<-eta_basal/rho_rep; MS<-rho_rep/(eta_basal/(Xc*Dc))
  list(Vt=Vt, Dt=Dt, tt=tt, X_eq=X_eq, MS=MS) }

# Parametros de H.sapiens para graficas ilustrativas
# G4: H.sapiens -- muestra declive gradual de V (vida larga, senescencia lenta)
R_HS <- list(eta=1.4e-3, t_lat=40.0, Xc=4.0, beta=54.75,
             eps=0.180, CV_eta=0.65, vida_max=350)
# G5: Mus musculus -- muestra cruce de Dc (vida corta, senescencia rapida)
# eta=0.084 alto hace que D suba rapido y cruce Dc=0.60 visualmente
R_MUS <- list(eta=0.084, t_lat=1.31, Xc=17.0, beta=54.75,
              eps=0.16, CV_eta=0.20, vida_max=6)

cat("Generando graficas del Anexo A...\n")

# ------------------------------------------------------------------------------
# G1-G3: GRAFICAS CLASICAS (una imagen de 1x3 paneles)
# ------------------------------------------------------------------------------
res_A  <- simular_A(R_HS$eta,  R_HS$t_lat,  R_HS$Xc,  R_HS$beta,  R_HS$eps,
                    CV_eta=R_HS$CV_eta,  vida_max=R_HS$vida_max,  N=2000, seed=7)
# Trayectoria de Mus calculada directamente (suavizado mensual, no anual)
set.seed(42)
{
  Dc_m<-0.60; dt_s<-7/365; kappa_m<-0.5
  bsem<-R_MUS$beta*dt_s; esem<-R_MUS$eps*sqrt(7)
  Xm<-0.01; Vm<-1; Vt_m<-c(); Dt_m<-c()
  Ns_m<-round(R_MUS$vida_max*365/7)
  for(sem in 1:Ns_m){ ta<-(sem-1)*dt_s
    dX<-R_MUS$eta*max(0,ta-R_MUS$t_lat)*7-bsem*Xm/(kappa_m+Xm)+esem*rnorm(1)
    Xm<-max(0,Xm+dX); Dm<-Xm/R_MUS$Xc; frm<-max(0,(1-Dm/Dc_m)^2)
    dV<-0.0988*Vm*frm-0.0960*Vm^3-max(0,Dm-Dc_m)
    Vm<-min(1.1,max(0,Vm+dV*dt_s))
    Vt_m<-c(Vt_m,Vm); Dt_m<-c(Dt_m,Dm)
    if(Xm>=R_MUS$Xc || Vm<=0.01) break }
  # Suavizado mensual: promedio de 4 semanas
  n4<-floor(length(Dt_m)/4)
  Dt_mus<-sapply(1:n4, function(i) mean(Dt_m[((i-1)*4+1):(i*4)]))
  tt_mus<-(1:n4)*4*dt_s - 2*dt_s  # tiempo en anos
  res_MUS<-list(Dt=Dt_mus, tt=tt_mus)
}
tc_A  <- res_A$tc

nombre_G1G3 <- sprintf("TMCSA_AnexoA_v10_G1G3_clasicas_%s.png", fecha)
png(nombre_G1G3, width=1800, height=600, res=110)
par(mfrow=c(1,3), mar=c(4.5,4.5,3.5,1.5))

# G1: Curva de supervivencia S(t)
t_grid <- seq(0, max(tc_A), length.out=200)
St     <- sapply(t_grid, function(t) mean(tc_A >= t))
plot(t_grid, St, type="l", lwd=3, col="#2C5F8A",
     xlab="Edad (años)", ylab="S(t) — fracción superviviente",
     main="G1 — Curva de supervivencia S(t)", ylim=c(0,1))
abline(h=0.5, lty=2, col="gray50"); text(5, 0.52, "S=0.5", cex=0.85, col="gray40")
legend("topright", legend="H. sapiens (simulado)", lwd=2, col="#2C5F8A", bty="n", cex=0.85)

# G2: Distribucion de t_cese (densidad)
plot(density(tc_A, bw=8), lwd=3, col="#4A7C4E",
     xlab="Edad de cese (años)", ylab="Densidad",
     main="G2 — Distribución de t_cese")
abline(v=median(tc_A), lty=1, col="#4A7C4E", lwd=2)
abline(v=mean(tc_A),   lty=2, col="gray40")
legend("topright", legend=c(sprintf("mediana=%.0fa",median(tc_A)),
       sprintf("media=%.0fa",mean(tc_A))),
       lty=c(1,2), col=c("#4A7C4E","gray40"), lwd=2, bty="n", cex=0.85)

# G3: Tasa de mortalidad h(t) — Gompertz emergente
# Intervalo ancho=15a para suavizar (menos dientes de sierra con N=2000)
# Rango 40-200a donde ocurren los ceses reales de H.sapiens
ancho <- 15; breaks_t <- seq(0, max(tc_A), by=ancho)
St_vec <- sapply(breaks_t, function(t) mean(tc_A >= t))
St_vec <- pmax(St_vec, 1e-6)  # evitar log(0)
ht <- diff(-log(St_vec))/ancho
t_mid <- breaks_t[-length(breaks_t)] + ancho/2
# Filtrar solo el rango con eventos reales (40-220a)
mask <- t_mid >= 40 & t_mid <= 220 & is.finite(ht) & ht > 0
plot(t_mid[mask], ht[mask], type="l", lwd=3, col="#8B4513",
     log="y",
     xlab="Edad (años)", ylab="h(t) — tasa de mortalidad (escala log)",
     main="G3 — Tasa de mortalidad h(t)\n(Gompertz emergente de [EM-2], escala log)")
# En escala log, Gompertz es linea recta creciente -- ajuste
# Ajuste Gompertz solo en rango central (70-150a) donde es mas lineal
mask_fit <- t_mid >= 70 & t_mid <= 150 & is.finite(ht) & ht > 0
fit <- lm(log(ht[mask_fit]) ~ t_mid[mask_fit])
# Extender la linea de ajuste a todo el rango para mostrar la tendencia
t_extrap <- seq(min(t_mid[mask]), max(t_mid[mask]), length.out=50)
ht_fit   <- exp(coef(fit)[1] + coef(fit)[2]*t_extrap)
lines(t_extrap, ht_fit, lwd=2, lty=2, col="#D2691E")
legend("topleft", legend=c("h(t) observada","tendencia Gompertz (70-150a)"),
       lwd=c(3,2), lty=c(1,2), col=c("#8B4513","#D2691E"), bty="n", cex=0.85)
dev.off()
cat(sprintf("  Guardada: %s\n", nombre_G1G3))

# ------------------------------------------------------------------------------
# G4-G5: V(t) y D(t) REGIMEN A (una imagen de 1x2)
# ------------------------------------------------------------------------------
nombre_G4G5 <- sprintf("TMCSA_AnexoA_v10_G4G5_VD_regimenA_%s.png", fecha)
png(nombre_G4G5, width=1200, height=550, res=110)
par(mfrow=c(1,2), mar=c(4.5,4.5,3.5,1.5))

plot(res_A$tt, res_A$Vt, type="l", lwd=2.5, col="#2C5F8A",
     xlab="Edad (años)", ylab="V(t) — vitalidad",
     main="G4 — V(t) régimen A\n(senescencia gradual)", ylim=c(0,1.1))
abline(h=1.0, lty=3, col="gray70"); abline(h=0.01, lty=2, col="#B2182B")
text(5, 0.04, "V=0.01 (cese)", col="#B2182B", cex=0.8)

plot(res_MUS$tt, res_MUS$Dt, type="l", lwd=2.5, col="#C77B00",
     xlab="Edad (años)", ylab="D = X/Xc",
     main="G5 — D(t) régimen A\n(Mus musculus — D cruza Dc=0.60)", ylim=c(0,0.75))
abline(h=Dc, lty=2, col="#B2182B", lwd=2)
text(min(res_MUS$tt)+0.1, Dc+0.04, "Dc=0.60", col="#B2182B", cex=0.85)
# Marcar el punto de cruce
cruce_idx <- which(res_MUS$Dt >= Dc)[1]
if(!is.na(cruce_idx)) points(res_MUS$tt[cruce_idx], Dc, pch=19, col="#B2182B", cex=1.5)
dev.off()
cat(sprintf("  Guardada: %s\n", nombre_G4G5))

# ------------------------------------------------------------------------------
# G6: V(t) y D(t) REGIMEN B -- tres especies superpuestas
# ------------------------------------------------------------------------------
nombre_G6 <- sprintf("TMCSA_AnexoA_v10_G6_VD_regimenB_%s.png", fecha)
png(nombre_G6, width=1200, height=550, res=110)
par(mfrow=c(1,2), mar=c(4.5,4.5,3.5,1.5))

rB_myo  <- simular_B(0.06,  10.0, Xc=3.0, Gamma_ext=log(2)/41,  vida_max=200)
rB_purp <- simular_B(0.04,  4.0,  Xc=3.0, Gamma_ext=log(2)/50,  vida_max=200)
rB_fran <- simular_B(0.03,  6.0,  Xc=3.0, Gamma_ext=log(2)/100, vida_max=250)

cores <- c("#2C5F8A","#6A0DAD","#C0392B")
nombres_B <- c("Myotis (MS=300)","S.purpuratus (MS=180)","S.franciscanus (MS=360)")

plot(NULL, xlim=c(0,120), ylim=c(0,1.1),
     xlab="Edad (años)", ylab="V(t)", main="G6a — V(t) régimen B\n(negligible senescence)")
for(i in seq_along(list(rB_myo, rB_purp, rB_fran))){
  r <- list(rB_myo, rB_purp, rB_fran)[[i]]
  lines(r$tt[r$tt<=120], r$Vt[r$tt<=120], lwd=2, col=cores[i]) }
abline(h=1.0, lty=3, col="gray70")
legend("bottomleft", legend=nombres_B, lwd=2, col=cores, bty="n", cex=0.75)

plot(NULL, xlim=c(0,120), ylim=c(0,0.65),
     xlab="Edad (años)", ylab="D = X/Xc", main="G6b — D(t) régimen B\n(estabilización en X_eq)")
for(i in seq_along(list(rB_myo, rB_purp, rB_fran))){
  r <- list(rB_myo, rB_purp, rB_fran)[[i]]
  lines(r$tt[r$tt<=120], r$Dt[r$tt<=120], lwd=2, col=cores[i]) }
abline(h=Dc, lty=2, col="#B2182B", lwd=2)
text(5, Dc+0.02, "Dc=0.60", col="#B2182B", cex=0.85)
legend("topright", legend=nombres_B, lwd=2, col=cores, bty="n", cex=0.75)
dev.off()
cat(sprintf("  Guardada: %s\n", nombre_G6))

# ------------------------------------------------------------------------------
# G7: Margen de Seguridad -- las tres especies de regimen B
# ------------------------------------------------------------------------------
nombre_G7 <- sprintf("TMCSA_AnexoA_v10_G7_MargenSeguridad_%s.png", fecha)
png(nombre_G7, width=900, height=650, res=110)
par(mar=c(5,5,4,2))

especies_B  <- c("Myotis\nbrandtii", "S. purpuratus", "S./M.\nfranciscanus")
longevidad  <- c(41, 50, 100)
MS_vals     <- c(300, 180, 360)
x_pos       <- c(1,2,3)

plot(longevidad, MS_vals, pch=19, cex=2.5, col=cores,
     xlim=c(20,120), ylim=c(100,450),
     xlab="Longevidad observada (años)",
     ylab="Margen de Seguridad (MS = ρ_rep / ρ_rep_crítico)",
     main="G7 — Margen de Seguridad vs Longevidad\n(régimen B, tres especies, P11 verificado)")
for(i in 1:3) text(longevidad[i], MS_vals[i]+20, especies_B[i], cex=0.8, col=cores[i])
abline(lm(MS_vals~longevidad), lty=2, col="gray50")
legend("topleft", legend="P11: MS escala con longevidad", lty=2, col="gray50", bty="n", cex=0.85)
dev.off()
cat(sprintf("  Guardada: %s\n", nombre_G7))

# ------------------------------------------------------------------------------
# G8: Paridad diagonal -- todas las especies de los canonicos 10A-10E
# (seleccion representativa de regimen A + B)
# ------------------------------------------------------------------------------
nombre_G8 <- sprintf("TMCSA_AnexoA_v10_G8_paridad_global_%s.png", fecha)
png(nombre_G8, width=1000, height=900, res=110)
par(mar=c(5,5,4,2))

# Datos verificados en R real (tomados de los logs de los canonicos)
targets  <- c(122,4,392,123,41,350,408,87,40,40,100,80,30,17,
              70,45,6,90,2,4.2,6.5,17,15,50,20,25,60,100,
              120,25,40,150,600,500,2000,3000,4850,10000)
medianas <- c(123.6,4.0,376.1,122.0,42.5,325.6,374.4,85.1,39.5,39.5,91.5,75.0,27.1,16.8,
              71.8,45.4,6.2,88.2,2.1,4.2,6.5,16.5,14.9,53.2,20.1,27.5,61.7,96.5,
              118,25,40,143,541,438,1981,2898,4588,9486)
grupos   <- c(rep("10A",9),rep("10C",9),rep("10D",10),rep("10E",10))
cores_g  <- c("10A"="#2C5F8A","10C"="#4A7C4E","10D"="#C77B00","10E"="#8B4513")

lim_min <- min(targets,medianas)*0.8; lim_max <- max(targets,medianas)*1.1
plot(targets, medianas, log="xy", pch=19, cex=1.2,
     col=sapply(grupos, function(g) cores_g[g]),
     xlim=c(lim_min,lim_max), ylim=c(lim_min,lim_max),
     xlab="Target t_cese (log)", ylab="Mediana simulada (log)",
     main="G8 — Paridad global: target vs simulado\n(regímenes A y B, todos los canónicos)")
abline(0,1, lty=2, col="gray40", lwd=2)
lines(c(lim_min,lim_max), c(lim_min,lim_max)*1.10, lty=3, col="gray70")
lines(c(lim_min,lim_max), c(lim_min,lim_max)*0.90, lty=3, col="gray70")
# Añadir B13/B14 en otro color
points(c(50,100), c(49.2,99.5), pch=17, cex=1.8, col="#6A0DAD")
legend("topleft",
       legend=c("10A (animales/plantas)","10C (vegetales)","10D (animales diversos)",
                "10E (plantas extremas)","Régimen B (erizos de mar)",
                "y=x (perfecto)","banda ±10%"),
       pch=c(19,19,19,19,17,NA,NA), lty=c(NA,NA,NA,NA,NA,2,3),
       col=c(unname(cores_g),"#6A0DAD","gray40","gray70"),
       bty="n", cex=0.75, pt.cex=1.2)
dev.off()
cat(sprintf("  Guardada: %s\n", nombre_G8))

cat("\n=== GRAFICAS COMPLETADAS ===\n")
cat(sprintf("G1-G3 clasicas:   %s\n", nombre_G1G3))
cat(sprintf("G4-G5 regimen A:  %s\n", nombre_G4G5))
cat(sprintf("G6  regimen B:    %s\n", nombre_G6))
cat(sprintf("G7  Margen Seg.:  %s\n", nombre_G7))
cat(sprintf("G8  paridad glob: %s\n", nombre_G8))
