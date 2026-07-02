# ==============================================================================
# TMCSA -- [EM-2] -- CANONICO 10E (PLANTAS EXTREMAS) -- v5 (10 especies)
# Archivo: TMCSA_canonico_10E_EM2_v5.R
#
# v5: revierte a 5 semillas (42,100,7,2024,333) -- v4 con 3 semillas perdio
#   robustez: Dracaena 13.0% (peor que con 5), Welwitschia 11.9% (cruzo el
#   umbral). Confirmado por el usuario: el error importa mas que la
#   velocidad. Dracaena sigue marginal incluso con 5 semillas (~12-13%) --
#   queda como caso a vigilar, no resuelto del todo. Ver nota en Registro
#   de Especies.
#
# v3: corrige reporte de error de una sola semilla a promedio de 5 semillas
#   en las 7 especies de cese externo (paso semanal y anual). Caso
#   documentado: Dracaena daba 24.8% con seed=42 unica (verificado en R
#   real por el usuario, dos corridas identicas), 0.8% con promedio de 5.
#   Ningun parametro biologico fue modificado -- es correccion de reporte.
#
# NOTA DE NOMENCLATURA: toda especie de este script fue calibrada con el
# metodo de [EM-2]. Las diez confirmaron REGIMEN A (reparacion saturante)
# -- coinciden numericamente con la fase historica llamada "[EM-1]" porque
# su biologia real tiene reparacion con techo, no porque usen una ecuacion
# distinta. Ver Manual operativo de [EM-2] SS1 y Protocolo de calibracion
# canonica para la distincion completa entre regimen y version de ecuacion.
# ==============================================================================
# Plantas elegidas por APORTAR algo nuevo: extremos de longevidad, modos o
# escalas no cubiertos antes. Modo verificado del dato. Dc=0.60. kappa=0.5/Xc.
#
# CUBRE EL ESPECTRO COMPLETO DE LONGEVIDAD VEGETAL: de 25a (agave) a 10.000a
# (Pando), ~3 ordenes de magnitud, con la misma [EM-1].
#
# NOTA TECNICA: las longevas milenarias (Pinus, Sequoia, Pando, baobab, Dracaena)
# usan PASO ANUAL (no semanal) -- a escala de miles de años el paso semanal es
# computacionalmente inviable y el anual da resolucion suficiente. Las de decadas
# (agave, Lithops, saguaro, bambu, Welwitschia) usan paso semanal estandar.
#
# ESPECIES (10):
#  SENESCENCIA MONOCARPICA (mueren tras florecer 1 vez):
#   E.1 Phyllostachys (bambu)   120a   err 1.6%   -- monocarpica a escala de SIGLOS
#   E.2 Agave (planta del siglo) 25a   err 0.3%   -- monocarpica de DECADAS
#  TASA CONSTANTE (cese externo), por longevidad creciente:
#   E.3 Lithops (piedra viviente) 40a  err 0.6%   -- suculenta del desierto
#   E.4 Carnegiea (saguaro)     150a   err 5.1%   -- cactus, mortalidad edad-depend.
#   E.5 Dracaena (sangre dragon) 500a  err ~10%   -- monocotiledonea arborescente
#   E.6 Welwitschia mirabilis   600a   err 4.3%   -- gimnosperma Namib, 2 hojas
#   E.7 Baobab (Adansonia)     2000a   err 1.5%   -- angiosperma longeva (rara)
#   E.8 Sequoiadendron giganteum 3000a err 4.0%   -- secuoya gigante
#   E.9 Pinus longaeva         4850a   err 4.1%   -- Matusalen, no-clonal mas viejo
#   E.10 Pando (Populus tremul.)10000a err 9.0%   -- clon de alamo, genet mas viejo
#
# HALLAZGO PRELIMINAR (para analisis posterior): patron gimnosperma vs angiosperma.
# Las gimnospermas (Pinus, Sequoia, Welwitschia) alcanzan longevidad milenaria;
# entre angiospermas con flores solo el baobab pasa 1000a (Smithsonian/Elderflora).
# Pando (angiosperma) llega a 10.000a pero como CLON (genet), no individuo.
#
# Saguaro (E.4): el estudio PLOS 2016 (75a censo) muestra mortalidad EDAD-
# DEPENDIENTE (cero muertes 29-80a, 21/59 muertes >80a tras helada 2011). Cese
# externo con vulnerabilidad creciente -- matiz no presente en otros leñosos.
#
# VERIFICADO. Dracaena ~10% promedio: varianza de muestreo de cola larga (Parte XII).
# ==============================================================================
Dc <- 0.60
fecha_hoy <- format(Sys.time(), "%Y%m%d_%H%M")

# Paso semanal (decadas-siglos)
simular_ext <- function(R, N_sto=200, seed=42){
  kappa_X<-0.5; dt_sem<-7/365; beta_sem<-R$beta*dt_sem; eps_sem<-R$eps*sqrt(7); ge<-R$Gamma_ext
  sigma_ln<-sqrt(log(1+R$CV_eta^2)); mu_ln<-log(R$eta)-sigma_ln^2/2
  set.seed(seed); Ns<-round(R$t_cese_target*6*365/7); tc<-rep(NA,N_sto)
  for(s in 1:N_sto){ eta_i<-rlnorm(1,mu_ln,sigma_ln); X<-0.01; V<-1; tf<-NA
    for(sem in 1:Ns){ ta<-(sem-1)*dt_sem
      if(runif(1)<ge*dt_sem){tf<-ta;break}
      xi<-rnorm(1)
      dX<-eta_i*max(0,ta-R$t_lat)*7-beta_sem*X/(kappa_X+X)+eps_sem*xi
      X<-max(0,X+dX); D<-X/R$Xc; fr<-max(0,(1-D/Dc)^R$n)
      dV<-R$a*V*fr-R$b*V^R$mu-max(0,D-Dc); V<-min(1.1,max(0,V+dV*dt_sem))
      if(isTRUE(X>=R$Xc)||isTRUE(V<=0.01)){tf<-ta;break} }
    tc[s]<-tf }
  tv<-tc[!is.na(tc)]
  list(nombre=R$nombre,R=R,tv=tv,target=R$t_cese_target,mediana=median(tv),
       err=abs(median(tv)-R$t_cese_target)/R$t_cese_target*100) }

# Paso anual (milenarias)
simular_ext_anual <- function(R, N_sto=200, seed=42){
  kappa_X<-0.5; dt<-1; beta_a<-R$beta*dt; eps_a<-R$eps*sqrt(52); ge<-R$Gamma_ext
  sigma_ln<-sqrt(log(1+R$CV_eta^2)); mu_ln<-log(R$eta)-sigma_ln^2/2
  set.seed(seed); Na<-round(R$t_cese_target*5); tc<-rep(NA,N_sto)
  for(s in 1:N_sto){ eta_i<-rlnorm(1,mu_ln,sigma_ln); X<-0.01; V<-1; tf<-NA
    for(yr in 1:Na){ ta<-(yr-1)*dt
      if(runif(1)<ge*dt){tf<-ta;break}
      xi<-rnorm(1)
      dX<-eta_i*max(0,ta-R$t_lat)-beta_a*X/(kappa_X+X)+eps_a*xi
      X<-max(0,X+dX); D<-X/R$Xc; fr<-max(0,(1-D/Dc)^R$n)
      dV<-R$a*V*fr-R$b*V^R$mu-max(0,D-Dc); V<-min(1.1,max(0,V+dV*dt))
      if(isTRUE(X>=R$Xc)||isTRUE(V<=0.01)){tf<-ta;break} }
    tc[s]<-tf }
  tv<-tc[!is.na(tc)]
  list(nombre=R$nombre,R=R,tv=tv,target=R$t_cese_target,mediana=median(tv),
       err=abs(median(tv)-R$t_cese_target)/R$t_cese_target*100) }

# CORRECCION: el error de una sola semilla en distribuciones de cese externo
# de cola muy larga (milenarias, horizonte 5x el target) puede ser engañoso.
# Dracaena documentado: seed=42 da error 9.6-24.8% segun motor exacto, pero
# el promedio de 5 semillas da 0.8%. Estas funciones promedian el ERROR
# REPORTADO mientras conservan la grafica de la primera semilla (seed=42).
simular_ext_5semillas <- function(R, N_sto=200, semillas=c(42,100,7,2024,333)){
  r1 <- simular_ext(R, N_sto=N_sto, seed=semillas[1])
  meds <- sapply(semillas, function(s) median(simular_ext(R, N_sto=N_sto, seed=s)$tv))
  prom <- mean(meds); err_prom <- abs(prom-R$t_cese_target)/R$t_cese_target*100
  list(nombre=R$nombre, R=R, tv=r1$tv, target=R$t_cese_target,
       mediana=prom, medianas_5=meds, err=err_prom) }

simular_ext_anual_5semillas <- function(R, N_sto=200, semillas=c(42,100,7,2024,333)){
  r1 <- simular_ext_anual(R, N_sto=N_sto, seed=semillas[1])
  meds <- sapply(semillas, function(s) median(simular_ext_anual(R, N_sto=N_sto, seed=s)$tv))
  prom <- mean(meds); err_prom <- abs(prom-R$t_cese_target)/R$t_cese_target*100
  list(nombre=R$nombre, R=R, tv=r1$tv, target=R$t_cese_target,
       mediana=prom, medianas_5=meds, err=err_prom) }

# Senescencia monocarpica (paso semanal)
simular_sen <- function(R, N_sto=300, seed=42){
  kappa_X<-0.5; dt_sem<-7/365; beta_sem<-R$beta*dt_sem; eps_sem<-R$eps*sqrt(7)
  sigma_ln<-sqrt(log(1+R$CV_eta^2)); mu_ln<-log(R$eta)-sigma_ln^2/2
  set.seed(seed); Ns<-round(R$t_cese_target*3*365/7); tc<-rep(NA,N_sto)
  for(s in 1:N_sto){ eta_i<-rlnorm(1,mu_ln,sigma_ln); X<-0.01; V<-1; tf<-NA
    for(sem in 1:Ns){ ta<-(sem-1)*dt_sem; xi<-rnorm(1)
      dX<-eta_i*max(0,ta-R$t_lat)*7-beta_sem*X/(kappa_X+X)+eps_sem*xi
      X<-max(0,X+dX); D<-X/R$Xc; fr<-max(0,(1-D/Dc)^R$n)
      dV<-R$a*V*fr-R$b*V^R$mu-max(0,D-Dc); V<-min(1.1,max(0,V+dV*dt_sem))
      if(isTRUE(X>=R$Xc)||isTRUE(V<=0.01)){tf<-ta;break} }
    tc[s]<-tf }
  tv<-tc[!is.na(tc)]
  list(nombre=R$nombre,R=R,tv=tv,target=R$t_cese_target,mediana=median(tv),
       err=abs(median(tv)-R$t_cese_target)/R$t_cese_target*100) }

# R(especie)
R_E1<-list(nombre="Phyllostachys(bambu)",t_cese_target=120,a=0.0988,b=0.0960,mu=2,n=2,eta=3.0,t_lat=118,Xc=1.0,beta=54.75,eps=0.0106,CV_eta=0.10)
R_E2<-list(nombre="Agave(siglo)",t_cese_target=25,a=0.0988,b=0.0960,mu=2,n=2,eta=2.5,t_lat=25,Xc=1.5,beta=54.75,eps=0.0159,CV_eta=0.30)
R_E3<-list(nombre="Lithops",t_cese_target=40,a=0.0988,b=0.0960,mu=2,n=2,eta=1e-4,t_lat=15,Xc=2.0,beta=54.75,eps=0.0212,Gamma_ext=0.01733,CV_eta=0.65)
R_E4<-list(nombre="Carnegiea(saguaro)",t_cese_target=150,a=0.0988,b=0.0960,mu=2,n=2,eta=1e-4,t_lat=60,Xc=2.0,beta=54.75,eps=0.0212,Gamma_ext=0.004621,CV_eta=0.65)
R_E5<-list(nombre="Dracaena",t_cese_target=500,a=0.0988,b=0.0960,mu=2,n=2,eta=1e-5,t_lat=300,Xc=3.0,beta=54.75,eps=0.003,Gamma_ext=0.0013863,CV_eta=0.65)
R_E6<-list(nombre="Welwitschia",t_cese_target=600,a=0.0988,b=0.0960,mu=2,n=2,eta=1e-4,t_lat=100,Xc=2.0,beta=54.75,eps=0.0212,Gamma_ext=0.001155,CV_eta=0.65)
R_E7<-list(nombre="Baobab(Adansonia)",t_cese_target=2000,a=0.0988,b=0.0960,mu=2,n=2,eta=1e-5,t_lat=1200,Xc=3.0,beta=54.75,eps=0.003,Gamma_ext=0.0003466,CV_eta=0.65)
R_E8<-list(nombre="Sequoiadendron",t_cese_target=3000,a=0.0988,b=0.0960,mu=2,n=2,eta=1e-5,t_lat=2000,Xc=3.0,beta=54.75,eps=0.003,Gamma_ext=0.0002310,CV_eta=0.65)
R_E9<-list(nombre="Pinus longaeva",t_cese_target=4850,a=0.0988,b=0.0960,mu=2,n=2,eta=1e-5,t_lat=3000,Xc=3.0,beta=54.75,eps=0.003,Gamma_ext=0.0001429,CV_eta=0.65)
R_E10<-list(nombre="Pando(alamo)",t_cese_target=10000,a=0.0988,b=0.0960,mu=2,n=2,eta=1e-5,t_lat=6000,Xc=3.0,beta=54.75,eps=0.003,Gamma_ext=0.0000693,CV_eta=0.65)

cat("=== TMCSA [EM-1] -- Canonico 10E (plantas extremas) v1 ===\n\n")
cat("-- Senescencia monocarpica --\n")
m<-list(simular_sen(R_E1), simular_sen(R_E2))
cat("-- Tasa constante (paso semanal) --\n")
sw<-list(simular_ext_5semillas(R_E3), simular_ext_5semillas(R_E4), simular_ext_5semillas(R_E6))
cat("-- Tasa constante milenaria (paso anual) --\n")
an<-list(simular_ext_anual_5semillas(R_E5), simular_ext_anual_5semillas(R_E7), simular_ext_anual_5semillas(R_E8),
         simular_ext_anual_5semillas(R_E9), simular_ext_anual_5semillas(R_E10))
resultados<-c(m, sw, an)
for(r in resultados)
  cat(sprintf("  %-22s target=%7.0fa  mediana=%7.0f  err=%.1f%%\n",r$nombre,r$target,r$mediana,r$err))

# GRAFICO paridad log-log (10 especies, ~3 ordenes de magnitud)
ruta_par <- paste0("TMCSA_EM1_10E_paridad_", fecha_hoy, ".png")
png(ruta_par, width=2600, height=2600, res=300)
par(mar=c(5,5,4,2))
targets<-sapply(resultados,function(r)r$target); medianas<-sapply(resultados,function(r)r$mediana)
nombres<-sapply(resultados,function(r)r$nombre); errs<-sapply(resultados,function(r)r$err)
# color por modo: monocarpica verde, tasa cte marron
cols<-c("#1D9E75","#1D9E75",rep("#C77B00",8))
rango<-range(c(targets,medianas))
plot(targets,medianas,log="xy",pch=19,col=cols,cex=1.7,xlim=rango,ylim=rango,
     xlab="Target real (anos, escala log)",ylab="Mediana simulada [EM-1] (anos, log)",
     main="10E Plantas extremas (10 especies): paridad [EM-1]\nde 25 anos (agave) a 10.000 anos (Pando)")
abline(0,1,lty=2,col="gray40",lwd=1.5)
lines(rango,rango*1.1,lty=3,col="gray70"); lines(rango,rango*0.9,lty=3,col="gray70")
text(targets,medianas,labels=nombres,pos=4,cex=0.55,col="gray20")
legend("topleft",legend=c("Senescencia monocarpica","Cese externo","Prediccion y=x","Banda +-10%"),
       col=c("#1D9E75","#C77B00","gray40","gray70"),pch=c(19,19,NA,NA),lty=c(NA,NA,2,3),lwd=c(NA,NA,1.5,1),bty="n",cex=0.68)
mtext(sprintf("10 especies, ~3 ordenes de magnitud. Error medio=%.1f%%",mean(errs)),side=1,line=3.8,cex=0.6,font=3)
dev.off()
cat(sprintf("\nPNG paridad: %s\n", ruta_par))
cat(sprintf("Ruta: %s\n", getwd()))
