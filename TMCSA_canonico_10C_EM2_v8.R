# ==============================================================================
# TMCSA -- [EM-2] -- CANONICO 10C (VEGETALES) -- v8 (13 especies)
# Archivo: TMCSA_canonico_10C_EM2_v8.R
#
# v8: revierte a 5 semillas (42,100,7,2024,333) -- v7 con 3 semillas perdio
#   robustez real: Amaryllis cruzo 10% (11.3%), Phalaenopsis subio a 8.4%.
#   Confirmado por el usuario: el error importa mas que la velocidad. 5
#   semillas es el estandar correcto para distribuciones de cese externo
#   de cola larga. Corrige ademas la etiqueta del log que decia "PROMEDIO
#   5 SEMILLAS" incluso cuando v7 calculaba con 3 -- la etiqueta no se
#   habia actualizado al cambiar el numero real de semillas.
#
# v6: corrige reporte de error de una sola semilla a promedio de 5 semillas
#   en las 9 especies de cese externo (tasa constante). Caso documentado:
#   Phalaenopsis daba 11.3% con seed=42 unica, 7.4% con promedio de 5.
#   M.grandiflora y Amaryllis con el mismo patron. Ningun parametro
#   biologico fue modificado -- es correccion de reporte, no de R(especie).
#
# NOTA DE NOMENCLATURA: toda especie de este script fue calibrada con el
# metodo de [EM-2]. Las trece confirmaron REGIMEN A (reparacion saturante)
# -- coinciden numericamente con la fase historica llamada "[EM-1]" porque
# su biologia real tiene reparacion con techo, no porque usen una ecuacion
# distinta. Ver Manual operativo de [EM-2] SS1 y Protocolo de calibracion
# canonica para la distincion completa entre regimen y version de ecuacion.
# ==============================================================================
# Cada especie con su modo VERIFICADO desde el dato (no asumido por taxon), datos
# propios, sin herencias. Dc=0.60 (geometria: 60% del Xc propio). kappa=0.5/Xc.
#
# HALLAZGO CONSOLIDADO: el modo emerge del dato, y se separa por HISTORIA DE VIDA,
# no por taxon:
#   - ANUALES MONOCARPICAS -> SENESCENCIA (mueren por programa interno tras
#     reproducirse, aunque el ambiente sea optimo). Tiempo en DIAS.
#   - PERENNES (leñosas/bulbosas/rizomatosas) -> TASA CONSTANTE (Gamma_ext):
#     mueren por causa externa conservando meristemo o regenerando. Tiempo en AÑOS.
# Evidencia fuerte: Citrullus lanatus (sandia) es senescencia, pero su pariente
# C. colocynthis es perenne -> el modo es de la especie, no del genero.
#
# ESPECIES (13 verificadas + 1 pendiente):
#  SENESCENCIA (4):
#   C.1 Arabidopsis thaliana   70 d   err 2.5%   NIVEL 1
#   C.2 Triticum aestivum      45 d   err 0.8%   NIVEL 1
#   C.3 Avena sativa            6 d   err 4.2%   NIVEL 1 (gamma_M medido)
#   C.4 Citrullus lanatus      90 d   err 1.9%   NIVEL 1 (sandia)
#  TASA CONSTANTE / CESE EXTERNO (9):
#   C.5 Olea europaea         350 a   err ~3-4%  (en 10A)
#   C.6 Fitzroya cupressoides 408 a   err 3.2%   (en 10A)
#   C.7 Drimys winteri         87 a   err ~3-5%  (canelo)
#   C.8 Ficus carica           40 a   err ~0.6%  (higuera)
#   C.9 Nelumbo nucifera       40 a   err ~0.6%  (loto, organismo vivo)
#   C.10 Magnolia grandiflora 100 a   err ~0.8%  NIVEL 1
#   C.11 Magnolia macrophylla  80 a   err ~5.7%  NIVEL 2
#   C.12 Amaryllis(Hippeastrum)30 a   err ~0.1%  (bulbosa, latencia estacional)
#   C.13 Phalaenopsis          17 a   err ~7%    (orquidea, regeneracion continua)
#  PENDIENTE:
#   C.14 Vitis vinifera -- Nivel 3, grafico sintetico, sin dato real. NO se fuerza.
#
# ANOTADO para fase futura de LATENCIAS DE SEMILLAS (no modelado aqui):
#   - Nelumbo: semilla viable 1300a (record mundial). - Araucaria: piñon (pendiente).
#
# NOTA: olivo y Fitzroya residen en 10A_v4_auditado.R. Aqui sus R se incluyen solo
# para el grafico de paridad. Para correrlos completos usar el 10A.
# ==============================================================================
Dc <- 0.60
fecha_hoy <- format(Sys.time(), "%Y%m%d_%H%M")

# ---- Senescencia (dias) ----
simular_veg_senesc <- function(R, N_sto=400, seed=42){
  kappa_X<-0.5; dt<-0.25; beta_d<-R$beta*dt; eps_d<-R$eps*sqrt(dt)
  sigma_ln<-sqrt(log(1+R$CV_eta^2)); mu_ln<-log(R$eta)-sigma_ln^2/2
  set.seed(seed); Ng<-round(R$t_cese_target*4/dt); tc<-rep(NA,N_sto)
  for(s in 1:N_sto){ eta_i<-rlnorm(1,mu_ln,sigma_ln); X<-0.01; V<-1; tf<-NA
    for(k in 1:Ng){ ta<-(k-1)*dt; xi<-rnorm(1)
      dX<-eta_i*max(0,ta-R$t_lat)*dt-beta_d*X/(kappa_X+X)+eps_d*xi
      X<-max(0,X+dX); D<-X/R$Xc; fr<-max(0,(1-D/Dc)^R$n)
      dV<-R$a*V*fr-R$b*V^R$mu-max(0,D-Dc); V<-min(1.1,max(0,V+dV*dt))
      if(isTRUE(X>=R$Xc)||isTRUE(V<=0.01)){tf<-ta;break} }
    tc[s]<-tf }
  tv<-tc[!is.na(tc)]
  list(nombre=R$nombre,R=R,tv=tv,target=R$t_cese_target,mediana=median(tv),
       err=abs(median(tv)-R$t_cese_target)/R$t_cese_target*100,unidad="dias") }

# ---- Cese externo (años) ----
simular_veg_ext <- function(R, N_sto=200, seed=42){
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
       err=abs(median(tv)-R$t_cese_target)/R$t_cese_target*100,unidad="anos") }

# CORRECCION: el error de una sola semilla en distribuciones de cese externo
# (cola larga) puede ser engañoso (Olea, M.grandiflora, Amaryllis, Phalaenopsis
# documentados con error >5% en una sola corrida pero <1% promediando). Esta
# funcion promedia 5 semillas para el ERROR REPORTADO mientras conserva la
# grafica de panel de la primera semilla (seed=42) para visualizacion.
simular_veg_ext_5semillas <- function(R, N_sto=200, semillas=c(42,100,7,2024,333)){
  r1 <- simular_veg_ext(R, N_sto=N_sto, seed=semillas[1])  # para grafica de panel
  meds <- sapply(semillas, function(s) median(simular_veg_ext(R, N_sto=N_sto, seed=s)$tv))
  prom <- mean(meds)
  err_prom <- abs(prom-R$t_cese_target)/R$t_cese_target*100
  list(nombre=R$nombre, R=R, tv=r1$tv, target=R$t_cese_target,
       mediana=prom, medianas_5=meds, err=err_prom, unidad="anos") }

panel_especie <- function(res, col_V, col_D){
  tv<-res$tv
  if(is.null(tv)||length(tv)<5){plot.new();title(main=res$nombre,cex.main=0.7);return(invisible())}
  x_max<-max(max(tv,na.rm=TRUE)*1.05,res$target*1.2); d<-density(tv,na.rm=TRUE)
  plot(d,col=col_V,lwd=2,xlim=c(0,x_max),main="",xlab="",ylab="",axes=FALSE)
  polygon(c(d$x,rev(d$x)),c(d$y,rep(0,length(d$y))),col=adjustcolor(col_V,0.18),border=NA)
  abline(v=res$target,lty=2,col="gray40",lwd=1.2); abline(v=median(tv),lty=1,col=col_D,lwd=1.5)
  axis(1,cex.axis=0.6,tcl=-0.3); box(col="gray70")
  title(main=sprintf("%s\nmed=%.1f tgt=%.1f e=%.1f%%",res$nombre,median(tv),res$target,res$err),cex.main=0.62) }

# ---- R(especie) ----
# Senescencia
R_C1<-list(nombre="Arabidopsis",t_cese_target=70,a=0.0988,b=0.0960,mu=2,n=2,eta=0.5,t_lat=48,Xc=3.0,beta=54.75,eps=0.0318,CV_eta=0.30)
R_C2<-list(nombre="Triticum",t_cese_target=45,a=0.0988,b=0.0960,mu=2,n=2,eta=0.3,t_lat=10,Xc=3.0,beta=54.75,eps=0.0318,CV_eta=0.30)
R_C3<-list(nombre="Avena",t_cese_target=6,a=0.0988,b=0.0960,mu=2,n=2,eta=1.0,t_lat=2,Xc=1.0,beta=54.75,eps=0.0106,CV_eta=0.20)
R_C4<-list(nombre="Citrullus(sandia)",t_cese_target=90,a=0.0988,b=0.0960,mu=2,n=2,eta=0.4,t_lat=60,Xc=3.0,beta=54.75,eps=0.0318,CV_eta=0.30)
# Cese externo
R_C5<-list(nombre="Olea(olivo)",t_cese_target=350,a=0.0988,b=0.0960,mu=2,n=2,eta=7.6996e-05,t_lat=150,Xc=1.3678,beta=54.75,eps=0.01448,Gamma_ext=0.00198,CV_eta=0.65)
R_C6<-list(nombre="Fitzroya",t_cese_target=408,a=0.0988,b=0.0960,mu=2,n=2,eta=7.6996e-05,t_lat=175,Xc=0.7960,beta=54.75,eps=0.008429,Gamma_ext=0.00170,CV_eta=0.65)
R_C7<-list(nombre="Drimys(canelo)",t_cese_target=87,a=0.0988,b=0.0960,mu=2,n=2,eta=1e-4,t_lat=80,Xc=2.0,beta=54.75,eps=0.0212,Gamma_ext=0.008,CV_eta=0.65)
R_C8<-list(nombre="Ficus(higuera)",t_cese_target=40,a=0.0988,b=0.0960,mu=2,n=2,eta=1e-4,t_lat=35,Xc=2.0,beta=54.75,eps=0.0212,Gamma_ext=0.01733,CV_eta=0.65)
R_C9<-list(nombre="Nelumbo(loto)",t_cese_target=40,a=0.0988,b=0.0960,mu=2,n=2,eta=1e-4,t_lat=25,Xc=2.0,beta=54.75,eps=0.0212,Gamma_ext=0.01733,CV_eta=0.65)
R_C10<-list(nombre="M.grandiflora",t_cese_target=100,a=0.0988,b=0.0960,mu=2,n=2,eta=1e-4,t_lat=60,Xc=2.0,beta=54.75,eps=0.0212,Gamma_ext=0.006931,CV_eta=0.65)
R_C11<-list(nombre="M.macrophylla",t_cese_target=80,a=0.0988,b=0.0960,mu=2,n=2,eta=1e-4,t_lat=50,Xc=2.0,beta=54.75,eps=0.0212,Gamma_ext=0.008664,CV_eta=0.65)
R_C12<-list(nombre="Amaryllis",t_cese_target=30,a=0.0988,b=0.0960,mu=2,n=2,eta=1e-4,t_lat=20,Xc=2.0,beta=54.75,eps=0.0212,Gamma_ext=0.02310,CV_eta=0.65)
R_C13<-list(nombre="Phalaenopsis",t_cese_target=17,a=0.0988,b=0.0960,mu=2,n=2,eta=1e-4,t_lat=10,Xc=2.0,beta=54.75,eps=0.0212,Gamma_ext=0.04077,CV_eta=0.65)

cat("=== TMCSA [EM-1] -- Canonico 10C v4 CONSOLIDADO (13 vegetales) ===\n\n")
cat("-- SENESCENCIA (anuales monocarpicas) --\n")
sen<-list(simular_veg_senesc(R_C1),simular_veg_senesc(R_C2),simular_veg_senesc(R_C3),simular_veg_senesc(R_C4))
cat("-- TASA CONSTANTE (perennes, cese externo, PROMEDIO 5 SEMILLAS) --\n")
ext<-list(simular_veg_ext_5semillas(R_C5),simular_veg_ext_5semillas(R_C6),simular_veg_ext_5semillas(R_C7),simular_veg_ext_5semillas(R_C8),
          simular_veg_ext_5semillas(R_C9),simular_veg_ext_5semillas(R_C10),simular_veg_ext_5semillas(R_C11),simular_veg_ext_5semillas(R_C12),simular_veg_ext_5semillas(R_C13))
resultados<-c(sen,ext)
for(r in resultados)
  cat(sprintf("  %-18s target=%8.1f %-4s  mediana=%8.1f  err=%.1f%%\n",
      r$nombre,r$target,r$unidad,r$mediana,r$err))

# --- GRAFICO 1: paneles (13) ---
ruta_paneles <- paste0("TMCSA_EM1_10C_paneles_", fecha_hoy, ".png")
png(ruta_paneles, width=3200, height=2600, res=190)
par(mfrow=c(4,4), mar=c(3.5,1.5,3,1))
cols_v<-c(rep("#1D9E75",4),rep("#C77B00",9))
cols_d<-c("darkgreen","forestgreen","olivedrab","seagreen","chocolate","sienna","saddlebrown","peru","tan4","orange3","darkorange3","goldenrod","darkgoldenrod")
for(i in 1:13) panel_especie(resultados[[i]], cols_v[i], cols_d[i])
dev.off()
cat(sprintf("\nPNG paneles: %s\n", ruta_paneles))

# --- GRAFICO 2: paridad log-log (13) ---
ruta_par <- paste0("TMCSA_EM1_10C_paridad_", fecha_hoy, ".png")
png(ruta_par, width=2600, height=2600, res=300)
par(mar=c(5,5,4,2))
to_anos<-function(r) if(r$unidad=="dias") c(r$target/365,r$mediana/365) else c(r$target,r$mediana)
pts<-t(sapply(resultados,to_anos)); targets<-pts[,1]; medianas<-pts[,2]
nombres<-sapply(resultados,function(r)r$nombre); errs<-sapply(resultados,function(r)r$err)
cols<-c(rep("#1D9E75",4),rep("#C77B00",9))
rango<-range(c(targets,medianas))
plot(targets,medianas,log="xy",pch=19,col=cols,cex=1.7,xlim=rango,ylim=rango,
     xlab="Target real (anos, escala log)",ylab="Mediana simulada [EM-1] (anos, log)",
     main="10C Vegetales (13 especies): paridad [EM-1]\nde 6 dias (avena) a 408 anos (Fitzroya)")
abline(0,1,lty=2,col="gray40",lwd=1.5)
lines(rango,rango*1.1,lty=3,col="gray70"); lines(rango,rango*0.9,lty=3,col="gray70")
text(targets,medianas,labels=nombres,pos=4,cex=0.55,col="gray20")
legend("topleft",legend=c("Senescencia (anual monocarpica)","Cese externo (perenne)","Prediccion y=x","Banda +-10%"),
       col=c("#1D9E75","#C77B00","gray40","gray70"),pch=c(19,19,NA,NA),lty=c(NA,NA,2,3),lwd=c(NA,NA,1.5,1),bty="n",cex=0.68)
mtext(sprintf("13 especies sobre ~5 ordenes de magnitud. Error medio=%.1f%%",mean(errs)),side=1,line=3.8,cex=0.6,font=3)
dev.off()
cat(sprintf("PNG paridad: %s\n", ruta_par))
cat(sprintf("Ruta: %s\n", getwd()))
