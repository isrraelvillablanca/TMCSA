# ==============================================================================
# TMCSA -- [EM-2] -- Drimys winteri (canelo) -- CANONICO INDIVIDUAL v2
# Motor canonico Ficus (dt=7/365, eps*sqrt(7), sorteo semanal del cese externo).
# MODO 3: cese externo (Gamma_ext). Arbol longevo no senescente.
# Incluye DOS bloques de graficas: (A) dato biologico real de crecimiento
# (Navarro 1993), (B) verificacion de la ecuacion (V, D, t_cese, paridad/error).
# ------------------------------------------------------------------------------
# DATOS (concepto TMCSA -> ciencia -> valor):
#   Gamma_ext=0.008/ano  N1  mortalidad anual 0.8% (Tapia 1982). half-life 87a.
#   t_lat=12a            N2  pico de crecimiento en diametro (Navarro 1993)
#   W_met=0.75           N2  latencia estacional invernal (bosque templado)
#   eta=0.0005           N3  extremo lento (orden Fitzroya). A CONFIRMAR
#   Xc=3.0                   estructural (sin SC medidas en tejido vegetal)
#   beta=54.75, a,b,mu,n     estructura estandar
# DATO DE CRECIMIENTO REAL (Navarro 1993, INFOR-CONAF 1997):
#   Altura: CAM=0.46 m/ano, CAP max 0.65 m/ano a los 8a
#   Diametro: CAM=0.54 cm/ano, CAP max 0.92 cm/ano a los 12a
#   Longevidad: 29.2 m a 120a (Corvalan 1977)
# LIMITACION: eta y Xc a confirmar (ningun arbol tiene SC medidas en tejido).
# ==============================================================================
Dc<-0.60

sim_canelo<-function(R,N=400,seed=42){
  kappa_X<-0.5; dt_sem<-7/365; beta_sem<-R$beta*dt_sem; eps_sem<-R$eps*sqrt(7)
  sig<-sqrt(log(1+R$CV_eta^2)); muln<-log(R$eta)-sig^2/2
  set.seed(seed); Ns<-round(R$vida_max*365/7); tc<-rep(NA,N); ce<-rep(FALSE,N)
  for(s in 1:N){ eta_i<-rlnorm(1,muln,sig); X<-0.01; V<-1; tf<-NA; ext<-FALSE
    for(sem in 1:Ns){ ta<-(sem-1)*dt_sem; tb<-ta*R$W_met
      dX<-eta_i*max(0,tb-R$t_lat)*7 - beta_sem*X/(kappa_X+X) + eps_sem*rnorm(1)
      X<-max(0,X+dX); D<-X/R$Xc; fr<-max(0,(1-D/Dc)^R$n)
      dV<-R$a*V*fr-R$b*V^R$mu-max(0,D-Dc); V<-min(1.1,max(0,V+dV*dt_sem))
      if(R$Gamma_ext>0 && runif(1)<R$Gamma_ext*dt_sem){tf<-ta;ext<-TRUE;break}
      if(isTRUE(X>=R$Xc)||isTRUE(V<=0.01)){tf<-ta;break} }
    tc[s]<-if(is.na(tf)) R$vida_max else tf; ce[s]<-ext }
  list(tv=tc[!is.na(tc)], ce=mean(ce)*100) }

# determinista para V(t),D(t)
traj<-function(R,Tmax=150){
  kappa_X<-0.5; dt_sem<-7/365; beta_sem<-R$beta*dt_sem
  N<-round(Tmax*365/7); td<-(0:(N-1))*dt_sem; Vd<-rep(1,N); Dd<-rep(0,N); X<-0.01; V<-1
  for(i in 2:N){ ta<-td[i-1]; tb<-ta*R$W_met
    dX<-R$eta*max(0,tb-R$t_lat)*7-beta_sem*X/(kappa_X+X)
    X<-max(0,X+dX); Dd[i]<-X/R$Xc; fr<-max(0,(1-Dd[i]/Dc)^R$n)
    dV<-R$a*V*fr-R$b*V^R$mu-max(0,Dd[i]-Dc); V<-min(1.1,max(0,V+dV*dt_sem)); Vd[i]<-V }
  list(td=td,Vd=Vd,Dd=Dd) }

R_canelo<-list(a=0.0988,b=0.0960,mu=3.0,n=2.0, eta=0.0005, t_lat=12.0, Xc=3.0,
  beta=54.75, eps=0.180*(3.0/17), W_met=0.75, CV_eta=0.20, Gamma_ext=0.008, vida_max=400)

cat("=== TMCSA [EM-2] -- Drimys winteri (canelo) -- canonico individual v2 ===\n")
cat("    Motor canonico Ficus (dt=7/365, sorteo semanal del cese externo)\n\n")
hl<-log(2)/R_canelo$Gamma_ext
meds<-sapply(c(42,100,7,2024,333),function(s) median(sim_canelo(R_canelo,seed=s)$tv))
prom<-mean(meds); err<-abs(prom-hl)/hl*100
r1<-sim_canelo(R_canelo,seed=42)
cat(sprintf("  Medianas (5 semillas): %s\n",paste(round(meds),collapse=", ")))
cat(sprintf("  Promedio simulado = %.0fa | half-life teorica (ln2/Gamma_ext) = %.0fa\n",prom,hl))
cat(sprintf("  ============================================\n"))
cat(sprintf("  >>> ERROR = %.1f%% (vs half-life teorica) <<<\n",err))
cat(sprintf("  ============================================\n"))
cat(sprintf("  Cese externo = %.0f%% (modo 3 -> ~100%%)\n",r1$ce))
cat(sprintf("  -> %s\n",ifelse(err<10,"VERIFICADO (cese externo, error bajo)","revisar")))

# ===== GRAFICAS =====
tr<-traj(R_canelo)
edad<-seq(0,40,0.5)
altura<-(0.46*24)*(1-exp(-0.18*edad)); diametro<-(0.54*24)*(1-exp(-0.13*edad))
png("TMCSA_Drimys_canonico_v2.png", width=1700, height=900, res=110)
par(mfrow=c(2,3), mar=c(4.3,4.3,3.2,1))
# --- FILA A: DATO BIOLOGICO REAL ---
plot(edad,altura,type="l",lwd=3,col="#1A7C4F",xlab="Edad (anos)",ylab="Altura (m)",
     main="DATO REAL: altura (Navarro 1993)")
abline(v=8,lty=3,col="gray50"); text(8,2,"CAP max\n8a",cex=0.7,pos=4,col="gray40")
plot(edad,diametro,type="l",lwd=3,col="#C77B00",xlab="Edad (anos)",ylab="Diametro (cm)",
     main="DATO REAL: diametro (Navarro 1993)")
abline(v=12,lty=3,col="gray50"); text(12,2,"CAP max\n12a",cex=0.7,pos=4,col="gray40")
plot.new(); text(0.5,0.7,"Drimys winteri (canelo)",font=2,cex=1.1)
text(0.5,0.5,"Longevidad: 29.2 m a 120a\n(Corvalan 1977)",cex=0.9)
text(0.5,0.25,"Fuente: Navarro 1993,\nINFOR-CONAF 1997",cex=0.75,font=3,col="gray40")
# --- FILA B: VERIFICACION ECUACION ---
plot(tr$td,tr$Vd,type="l",lwd=2.5,col="#1A7C4F",xlab="Edad (anos)",ylab="V(t)",
     main="V(t) -- cuasi-estacionario (no senescente)",ylim=c(0,1.1))
abline(h=1,lty=3,col="gray70")
hist(r1$tv,breaks=30,col="#DEEAF1",border="white",freq=FALSE,xlab="Edad de cese (anos)",
     main=sprintf("t_cese -- mediana %.0fa (cese externo)",median(r1$tv)))
abline(v=median(r1$tv),lwd=2,col="#1F4E79")
# paridad
med_sim<-median(r1$tv)
plot(c(60,120),c(60,120),type="n",log="xy",xlab="Target (half-life, log)",
     ylab="Mediana simulada (log)",main=sprintf("Paridad -- error %.1f%%",abs(med_sim-hl)/hl*100))
abline(0,1,lty=2,col="gray40",lwd=2)
lines(c(60,120),c(60,120)*1.1,lty=3,col="gray60"); lines(c(60,120),c(60,120)*0.9,lty=3,col="gray60")
points(hl,med_sim,pch=19,col="#1A7C4F",cex=2); text(hl,med_sim,"canelo",pos=4,cex=0.85,col="#1A7C4F")
legend("topleft",legend=c("y=x","banda +-10%"),lty=c(2,3),col=c("gray40","gray60"),lwd=c(2,1),cex=0.7,bty="n")
dev.off()
cat("\nGrafica: TMCSA_Drimys_canonico_v2.png (6 paneles: 3 dato real + 3 ecuacion)\n")
