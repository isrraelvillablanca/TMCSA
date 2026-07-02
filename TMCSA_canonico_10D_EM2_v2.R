# ==============================================================================
# TMCSA -- [EM-2] -- CANONICO 10D (ANIMALES DIVERSOS) -- v2 (10 especies)
# Archivo: TMCSA_canonico_10D_EM2_v2.R
#
# NOTA DE NOMENCLATURA: toda especie de este script fue calibrada con el
# metodo de [EM-2]. Las diez confirmaron REGIMEN A (reparacion saturante)
# -- coinciden numericamente con la fase historica llamada "[EM-1]" porque
# su biologia real tiene reparacion con techo, no porque usen una ecuacion
# distinta. Ver Manual operativo de [EM-2] SS1 y Protocolo de calibracion
# canonica para la distincion completa entre regimen y version de ecuacion.
# ==============================================================================
# Animales de taxones y longevidades muy diversas, del Anexo (verificados) +
# nuevos. Modo K basico (senescencia Gompertziana) salvo notas. Sin herencias.
# Dc=0.60 (geometria). kappa=0.5/Xc por especie.
#
# Demuestra que [EM-1] captura desde una abeja reina (2a) hasta una tortuga de
# Galapagos (100a) -- ~2 ordenes de magnitud, invertebrados y vertebrados,
# con la misma ecuacion. El modo se determina del dato en cada caso.
#
# ESPECIES (10):
#  INVERTEBRADOS (Anexo, NIVEL 1):
#   D.1 Apis mellifera (reina)   2.0a   err 5.5%  -- K basico, casta reina
#   D.2 Eisenia fetida (lombriz) 4.25a  err 0.3%  -- K basico, Gompertz verificado
#   D.3 Homarus (langosta)       6.5a   err 0.3%  -- K basico, reloj anual
#   D.4 Magicicada (cigarra)    17.0a   err 2.8%  -- semelpara + W_met ninfal
#  VERTEBRADOS (Anexo, NIVEL 2):
#   D.5 Xenopus (rana)          15.0a   err 3.6%  -- K basico + W_met sequia
#   D.6 Alligator (caiman)      50.0a   err 5.2%  -- K basico (NIVEL 2)
#  VERTEBRADOS (nuevos):
#   D.7 Aptenodytes (pinguino)  20.0a   err 1.9%  -- K basico, polar
#   D.8 Ursus maritimus (oso)   25.0a   err 10.7% -- K basico, polar
#   D.9 Andrias (salamandra)    60.0a   err 3.0%  -- K basico (contraste axolote)
#   D.10 Chelonoidis niger      100.0a  err 0.4%  -- K basico, tortuga Galapagos
#
# NOTA: Magicicada -- el cese a 17a es por semelparia (muerte tras emerger y
# reproducirse). La latencia ninfal subterranea (16a) es W_met bajo; se modela
# con t_lat=16 (acumulacion neta inicia al final de la fase ninfal). Andrias es
# pariente del axolote (Ambystoma, caso abierto en 10B) pero SIN su disputa de
# modo: Andrias senesce normalmente -> contraste util.
# ==============================================================================
Dc <- 0.60
fecha_hoy <- format(Sys.time(), "%Y%m%d_%H%M")

simular_K <- function(R, N_sto=400, seed=42){
  kappa_X<-0.5; dt_sem<-7/365; beta_sem<-R$beta*dt_sem; eps_sem<-R$eps*sqrt(7)
  sigma_ln<-sqrt(log(1+R$CV_eta^2)); mu_ln<-log(R$eta)-sigma_ln^2/2
  set.seed(seed); Ns<-round(R$t_cese_target*4*365/7); tc<-rep(NA,N_sto)
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

panel_especie <- function(res, col_V, col_D){
  tv<-res$tv
  if(is.null(tv)||length(tv)<5){plot.new();title(main=res$nombre,cex.main=0.7);return(invisible())}
  x_max<-max(max(tv,na.rm=TRUE)*1.05,res$target*1.2); d<-density(tv,na.rm=TRUE)
  plot(d,col=col_V,lwd=2,xlim=c(0,x_max),main="",xlab="",ylab="",axes=FALSE)
  polygon(c(d$x,rev(d$x)),c(d$y,rep(0,length(d$y))),col=adjustcolor(col_V,0.18),border=NA)
  abline(v=res$target,lty=2,col="gray40",lwd=1.2); abline(v=median(tv),lty=1,col=col_D,lwd=1.5)
  axis(1,cex.axis=0.6,tcl=-0.3); box(col="gray70")
  title(main=sprintf("%s\nmed=%.1f tgt=%.1f e=%.1f%%",res$nombre,median(tv),res$target,res$err),cex.main=0.6) }

# R(especie)
R_D1 <-list(nombre="Apis(reina)",a=0.0988,b=0.0960,mu=3,n=2,t_lat=0.8,Xc=3.0,eta=0.12,beta=54.75,CV_eta=0.30,t_cese_target=2.0,eps=0.0318)
R_D2 <-list(nombre="Eisenia(lombriz)",a=0.0988,b=0.0960,mu=3,n=2,t_lat=2.0,Xc=3.0,eta=0.065,beta=54.75,CV_eta=0.30,t_cese_target=4.25,eps=0.0318)
R_D3 <-list(nombre="Homarus(langosta)",a=0.0988,b=0.0960,mu=3,n=2,t_lat=3.5,Xc=3.0,eta=0.05,beta=54.75,CV_eta=0.30,t_cese_target=6.5,eps=0.0318)
R_D4 <-list(nombre="Magicicada(cigarra)",a=0.0988,b=0.0960,mu=3,n=2,t_lat=16.5,Xc=1.0,eta=5.0,beta=54.75,CV_eta=0.10,t_cese_target=17.0,eps=0.0106)
R_D5 <-list(nombre="Xenopus(rana)",a=0.0988,b=0.0960,mu=3,n=2,t_lat=3.0,Xc=3.0,eta=0.012,beta=54.75,CV_eta=0.40,t_cese_target=15.0,eps=0.0318)
R_D6 <-list(nombre="Alligator(caiman)",a=0.0988,b=0.0960,mu=3,n=2,t_lat=10.0,Xc=3.0,eta=0.0032,beta=54.75,CV_eta=0.40,t_cese_target=50.0,eps=0.0318)
R_D7 <-list(nombre="Aptenodytes(pinguino)",a=0.0988,b=0.0960,mu=3,n=2,t_lat=5.0,Xc=3.0,eta=0.009,beta=54.75,CV_eta=0.40,t_cese_target=20.0,eps=0.0318)
R_D8 <-list(nombre="Ursus(oso polar)",a=0.0988,b=0.0960,mu=3,n=2,t_lat=6.0,Xc=3.0,eta=0.0065,beta=54.75,CV_eta=0.40,t_cese_target=25.0,eps=0.0318)
R_D9 <-list(nombre="Andrias(salamandra)",a=0.0988,b=0.0960,mu=3,n=2,t_lat=10.0,Xc=3.0,eta=0.0025,beta=54.75,CV_eta=0.40,t_cese_target=60.0,eps=0.0318)
R_D10<-list(nombre="Chelonoidis niger",a=0.0988,b=0.0960,mu=3,n=2,t_lat=15.0,Xc=3.0,eta=0.0016,beta=54.75,CV_eta=0.40,t_cese_target=100.0,eps=0.0318)

cat("=== TMCSA [EM-1] -- Canonico 10D (animales diversos) v1 ===\n\n")
Rs<-list(R_D1,R_D2,R_D3,R_D4,R_D5,R_D6,R_D7,R_D8,R_D9,R_D10)
resultados<-lapply(Rs, simular_K)
for(r in resultados)
  cat(sprintf("  %-22s target=%6.1fa  mediana=%6.1f  err=%.1f%%\n",r$nombre,r$target,r$mediana,r$err))

# GRAFICO 1: paneles
ruta_paneles <- paste0("TMCSA_EM1_10D_paneles_", fecha_hoy, ".png")
png(ruta_paneles, width=3000, height=2400, res=190)
par(mfrow=c(3,4), mar=c(3.5,1.5,3,1))
cols_v<-c(rep("#2166AC",4),rep("#762A83",6))  # azul invertebrados, morado vertebrados
cols_d<-c("navy","blue3","dodgerblue4","steelblue","purple4","darkviolet","mediumorchid4","magenta4","orchid4","purple")
for(i in 1:10) panel_especie(resultados[[i]], cols_v[i], cols_d[i])
dev.off()
cat(sprintf("\nPNG paneles: %s\n", ruta_paneles))

# GRAFICO 2: paridad log-log
ruta_par <- paste0("TMCSA_EM1_10D_paridad_", fecha_hoy, ".png")
png(ruta_par, width=2600, height=2600, res=300)
par(mar=c(5,5,4,2))
targets<-sapply(resultados,function(r)r$target); medianas<-sapply(resultados,function(r)r$mediana)
nombres<-sapply(resultados,function(r)r$nombre); errs<-sapply(resultados,function(r)r$err)
cols<-c(rep("#2166AC",4),rep("#762A83",6))
rango<-range(c(targets,medianas))
plot(targets,medianas,log="xy",pch=19,col=cols,cex=1.7,xlim=rango,ylim=rango,
     xlab="Target real (anos, escala log)",ylab="Mediana simulada [EM-1] (anos, log)",
     main="10D Animales (10 especies): paridad [EM-1]\nde 2 anos (abeja) a 100 anos (tortuga Galapagos)")
abline(0,1,lty=2,col="gray40",lwd=1.5)
lines(rango,rango*1.1,lty=3,col="gray70"); lines(rango,rango*0.9,lty=3,col="gray70")
text(targets,medianas,labels=nombres,pos=4,cex=0.55,col="gray20")
legend("topleft",legend=c("Invertebrados","Vertebrados","Prediccion y=x","Banda +-10%"),
       col=c("#2166AC","#762A83","gray40","gray70"),pch=c(19,19,NA,NA),lty=c(NA,NA,2,3),lwd=c(NA,NA,1.5,1),bty="n",cex=0.68)
mtext(sprintf("10 especies, invertebrados y vertebrados. Error medio=%.1f%%",mean(errs)),side=1,line=3.8,cex=0.6,font=3)
dev.off()
cat(sprintf("PNG paridad: %s\n", ruta_par))
cat(sprintf("Ruta: %s\n", getwd()))
