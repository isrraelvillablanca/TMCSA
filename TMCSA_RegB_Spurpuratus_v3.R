# ==============================================================================
# TMCSA -- [EM-2] -- Strongylocentrotus purpuratus (erizo de mar purpura) -- v3
# REGIMEN B (proporcional) + cese externo (Gamma_ext). Motor v2 (orden Gext
# primero, identico al canonico simular_K real). v3 anade graficas completas.
# ------------------------------------------------------------------------------
# COMPLEMENTA, NO REEMPLAZA, a B.EM1-7 del canonico 10B original (target=6.93a,
# Gamma_ext=0.10): esa ficha mide MORTALIDAD POBLACIONAL bajo depredacion
# marina real (mismo patron que Euphausia, "presa base"). Esta ficha mide
# la CAPACIDAD BIOLOGICA DE NO-SENESCENCIA en condiciones sin depredacion
# (laboratorio/acuario, Bodnar & Coffman 2016). Ambas son correctas.
#
# MODO (Paso 0): negligible senescence (regimen B) + cese por Gamma_ext (modo 3).
#   Du et al. 2013 / Bodnar & Coffman 2016 (Mech Ageing Dev): protein
#   carbonyls, 4-HNE, 8-OHdG sin incremento general con la edad en 6 tejidos
#   (musculo, nervio, esofago, gonada, celomocitos, ampulas). Capacidad
#   antioxidante (SOD, total) y proteasoma MANTENIDAS con la edad.
#
# DATOS (concepto -> ciencia -> valor):
#   eta_basal=0.04   produccion basal -- NIVEL 3 (censo, sin cuantificacion
#                    absoluta de carbonilos/tiempo)
#   rho_rep=4.0      reparacion proporcional -- NIVEL 3
#                    -> X_eq=0.0100, MS=180x -> no colapsa por daño
#   Xc=3.0 estructural (molde Myotis brandtii); Dc=0.60
#   Gamma_ext=ln(2)/50=0.01386 -> half-life=50a (NIVEL 2, capacidad biologica
#     en cautiverio, Bodnar 2011-2016)
#   CV_eta=0.30      poblacion silvestre, heterogeneidad intermedia
# CORRECCION v2->v3: orden de Gamma_ext PRIMERO en el bucle (identico al
#   motor simular_K). Verificado por el usuario en R real: error 4.4% (v1,
#   orden anterior). v2/v3 corrigen el orden -- error esperado menor.
# LIMITACION: eta_basal y rho_rep son estimados de censo (NIVEL 3), no
#   medidos directamente. Pendiente: cuantificacion absoluta comparativa.
# ==============================================================================
Dc<-0.60

sim_individuo<-function(R,N=500,seed=42){
  dt_sem<-7/365; eps_sem<-R$eps*sqrt(7); Gext_sem<-R$Gamma_ext*dt_sem
  sig<-sqrt(log(1+R$CV_eta^2)); muln<-log(R$eta_basal)-sig^2/2
  set.seed(seed); Ns<-round(R$vida_max*365/7); tc<-rep(NA,N); ce<-rep(FALSE,N)
  for(s in 1:N){ eta_i<-rlnorm(1,muln,sig); X<-0.01; V<-1; tf<-NA; ext<-FALSE
    for(sem in 1:Ns){ ta<-(sem-1)*dt_sem
      if(Gext_sem>0 && runif(1)<Gext_sem){tf<-ta;ext<-TRUE;break}
      dX<-(eta_i - R$rho_rep*X)*dt_sem + eps_sem*rnorm(1)   # REGIMEN B
      X<-max(0,X+dX); D<-X/R$Xc; fr<-max(0,(1-D/Dc)^R$n)
      dV<-R$a*V*fr-R$b*V^R$mu-max(0,D-Dc); V<-min(1.1,max(0,V+dV*dt_sem))
      if(isTRUE(X>=R$Xc)||isTRUE(V<=0.01)){tf<-ta;break} }
    tc[s]<-if(is.na(tf)) R$vida_max else tf; ce[s]<-ext }
  list(tv=tc[!is.na(tc)], ce=mean(ce)*100) }

traj<-function(R,Tmax=80){
  dt_sem<-7/365; N<-round(Tmax*365/7); td<-(0:(N-1))*dt_sem; Vd<-rep(1,N); Dd<-rep(0,N); X<-0.01; V<-1
  for(i in 2:N){ dX<-(R$eta_basal - R$rho_rep*X)*dt_sem
    X<-max(0,X+dX); Dd[i]<-X/R$Xc; fr<-max(0,(1-Dd[i]/Dc)^R$n)
    dV<-R$a*V*fr-R$b*V^R$mu-max(0,Dd[i]-Dc); V<-min(1.1,max(0,V+dV*dt_sem)); Vd[i]<-V }
  list(td=td,Vd=Vd,Dd=Dd) }

R_purp<-list(nombre="Strongylocentrotus purpuratus (capacidad biologica)",
  a=0.0988,b=0.0960,mu=3.0,n=2.0, eta_basal=0.04, rho_rep=4.0,
  Xc=3.0, eps=0.180*(3.0/17), W_met=1.0, CV_eta=0.30,
  Gamma_ext=log(2)/50, t_cese_target=50, vida_max=300)

cat("=== TMCSA [EM-2] -- S. purpuratus -- v3 (regimen B, con graficas) ===\n\n")
X_eq<-R_purp$eta_basal/R_purp$rho_rep
rho_crit<-R_purp$eta_basal/(R_purp$Xc*Dc)
MS_val<-R_purp$rho_rep/rho_crit
cat(sprintf("  X_eq=%.4f | Dc*Xc=%.2f | MS=%.1f -> %s\n",
    X_eq, Dc*R_purp$Xc, MS_val, ifelse(X_eq<Dc*R_purp$Xc,"no colapsa por daño (correcto)","COLAPSA")))
meds<-sapply(c(42,100,7,2024,333),function(s) median(sim_individuo(R_purp,seed=s)$tv))
prom<-mean(meds); err<-abs(prom-R_purp$t_cese_target)/R_purp$t_cese_target*100
r1<-sim_individuo(R_purp,seed=42)
cat(sprintf("  Medianas (5 semillas): %s\n",paste(round(meds),collapse=", ")))
cat(sprintf("  Promedio=%.1fa | Target=%.0fa (NIVEL 2)\n",prom,R_purp$t_cese_target))
cat(sprintf("  ============================================\n"))
cat(sprintf("  >>> ERROR = %.1f%% <<<\n",err))
cat(sprintf("  ============================================\n"))
cat(sprintf("  Cese externo=%.0f%%\n",r1$ce))
cat(sprintf("  -> %s\n",ifelse(err<10,"VERIFICADO (regimen B, NIVEL 2)","revisar")))

tr<-traj(R_purp)
fecha<-format(Sys.Date(),"%Y%m%d")
png(sprintf("TMCSA_Spurpuratus_canonico_v3_%s.png",fecha), width=1700, height=460, res=110)
par(mfrow=c(1,4), mar=c(4.3,4.3,3.2,1))
plot(tr$td,tr$Vd,type="l",lwd=2.5,col="#6A0DAD",xlab="Edad (años)",ylab="V(t)",
     main="V(t) -- no senescente (régimen B)",ylim=c(0,1.1)); abline(h=1,lty=3,col="gray70")
plot(tr$td,tr$Dd,type="l",lwd=2.5,col="#C77B00",xlab="Edad (años)",ylab="D=X/Xc",
     main="D(t) -- estabiliza en X_eq (no cruza Dc)",ylim=c(0,Dc*1.1))
abline(h=Dc,lty=2,col="#B2182B"); text(60,Dc-0.05,"Dc=0.60",col="#B2182B",cex=0.8)
abline(h=X_eq/R_purp$Xc,lty=3,col="#1A7C4F"); text(60,0.05,"X_eq (equilibrio)",col="#1A7C4F",cex=0.7)
hist(r1$tv,breaks=40,col="#E6D9F2",border="white",freq=FALSE,xlab="Edad de cese (años)",
     main=sprintf("t_cese -- mediana %.0fa (Γ_ext)",median(r1$tv)), xlim=c(0,150))
abline(v=median(r1$tv),lwd=2,col="#4B0082")
plot(c(30,75),c(30,75),type="n",log="xy",xlab="Target (half-life, log)",ylab="Mediana sim (log, promedio 5 semillas)",
     main=sprintf("Paridad -- error %.1f%% (promedio 5 semillas)",err))
abline(0,1,lty=2,col="gray40",lwd=2); lines(c(30,75),c(30,75)*1.1,lty=3,col="gray60"); lines(c(30,75),c(30,75)*0.9,lty=3,col="gray60")
points(R_purp$t_cese_target,prom,pch=19,col="#6A0DAD",cex=2); text(R_purp$t_cese_target,prom,"S.purpuratus",pos=4,cex=0.85,col="#6A0DAD")
legend("topleft",legend=c("y=x","banda +-10%"),lty=c(2,3),col=c("gray40","gray60"),lwd=c(2,1),cex=0.7,bty="n")
dev.off()
cat("\nGrafica: TMCSA_Spurpuratus_canonico_v3.png\n")
