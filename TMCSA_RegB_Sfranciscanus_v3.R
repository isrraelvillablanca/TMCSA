# ==============================================================================
# TMCSA -- [EM-2] -- S./M. franciscanus (erizo de mar rojo) -- v3
# REGIMEN B (proporcional) + cese externo (Gamma_ext). Motor v2 (orden Gext
# primero, identico al canonico simular_K real). v3 anade graficas completas.
# ------------------------------------------------------------------------------
# MODO (Paso 0): negligible senescence robusta (regimen B) + Gamma_ext (modo 3).
#   "Uno de los animales mas longevos de la Tierra... sin incremento de
#   mortalidad relacionado con la edad ni declive de capacidad reproductiva"
#   (Ebert; Bodnar 2015). Evidencia MAS fuerte que S. purpuratus: menor
#   daño oxidativo relativo y mayor estabilidad transcripcional (Bodnar &
#   Coffman 2016).
#
# DATOS (concepto -> ciencia -> valor):
#   eta_basal=0.03   produccion basal -- NIVEL 3 (menor que purpuratus,
#                    coherente con menor daño oxidativo relativo medido)
#   rho_rep=6.0      reparacion proporcional -- NIVEL 3 (mayor que
#                    purpuratus, coherente con mayor longevidad)
#                    -> X_eq=0.0050, MS=360x (2x el margen de purpuratus)
#   Xc=3.0 estructural; Dc=0.60
#   Gamma_ext=ln(2)/100=0.006931 -> half-life=100a (NIVEL 2, Ebert)
#   CV_eta=0.25 (estimado, ligeramente menor que purpuratus)
# CORRECCION v2->v3: orden de Gamma_ext PRIMERO en el bucle. Verificado por
#   el usuario en R real: error 1.0% (v1, orden anterior).
# LIMITACION: eta_basal y rho_rep estimados de censo (NIVEL 3). Pendiente
#   cuantificacion absoluta comparativa con S. purpuratus.
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

traj<-function(R,Tmax=150){
  dt_sem<-7/365; N<-round(Tmax*365/7); td<-(0:(N-1))*dt_sem; Vd<-rep(1,N); Dd<-rep(0,N); X<-0.01; V<-1
  for(i in 2:N){ dX<-(R$eta_basal - R$rho_rep*X)*dt_sem
    X<-max(0,X+dX); Dd[i]<-X/R$Xc; fr<-max(0,(1-Dd[i]/Dc)^R$n)
    dV<-R$a*V*fr-R$b*V^R$mu-max(0,Dd[i]-Dc); V<-min(1.1,max(0,V+dV*dt_sem)); Vd[i]<-V }
  list(td=td,Vd=Vd,Dd=Dd) }

R_fran<-list(nombre="S./M. franciscanus (capacidad biologica)",
  a=0.0988,b=0.0960,mu=3.0,n=2.0, eta_basal=0.03, rho_rep=6.0,
  Xc=3.0, eps=0.180*(3.0/17), W_met=1.0, CV_eta=0.25,
  Gamma_ext=log(2)/100, t_cese_target=100, vida_max=400)

cat("=== TMCSA [EM-2] -- S./M. franciscanus -- v3 (regimen B, con graficas) ===\n\n")
X_eq<-R_fran$eta_basal/R_fran$rho_rep
rho_crit<-R_fran$eta_basal/(R_fran$Xc*Dc)
MS_val<-R_fran$rho_rep/rho_crit
cat(sprintf("  X_eq=%.4f | Dc*Xc=%.2f | MS=%.1f -> %s\n",
    X_eq, Dc*R_fran$Xc, MS_val, ifelse(X_eq<Dc*R_fran$Xc,"no colapsa por daño (correcto)","COLAPSA")))
meds<-sapply(c(42,100,7,2024,333),function(s) median(sim_individuo(R_fran,seed=s)$tv))
prom<-mean(meds); err<-abs(prom-R_fran$t_cese_target)/R_fran$t_cese_target*100
r1<-sim_individuo(R_fran,seed=42)
cat(sprintf("  Medianas (5 semillas): %s\n",paste(round(meds),collapse=", ")))
cat(sprintf("  Promedio=%.1fa | Target=%.0fa (NIVEL 2)\n",prom,R_fran$t_cese_target))
cat(sprintf("  ============================================\n"))
cat(sprintf("  >>> ERROR = %.1f%% <<<\n",err))
cat(sprintf("  ============================================\n"))
cat(sprintf("  Cese externo=%.0f%%\n",r1$ce))
cat(sprintf("  -> %s\n",ifelse(err<10,"VERIFICADO (regimen B, NIVEL 2)","revisar")))

tr<-traj(R_fran)
fecha<-format(Sys.Date(),"%Y%m%d")
png(sprintf("TMCSA_Sfranciscanus_canonico_v3_%s.png",fecha), width=1700, height=460, res=110)
par(mfrow=c(1,4), mar=c(4.3,4.3,3.2,1))
plot(tr$td,tr$Vd,type="l",lwd=2.5,col="#C0392B",xlab="Edad (años)",ylab="V(t)",
     main="V(t) -- no senescente (régimen B)",ylim=c(0,1.1)); abline(h=1,lty=3,col="gray70")
plot(tr$td,tr$Dd,type="l",lwd=2.5,col="#C77B00",xlab="Edad (años)",ylab="D=X/Xc",
     main="D(t) -- estabiliza en X_eq (no cruza Dc)",ylim=c(0,Dc*1.1))
abline(h=Dc,lty=2,col="#B2182B"); text(120,Dc-0.05,"Dc=0.60",col="#B2182B",cex=0.8)
abline(h=X_eq/R_fran$Xc,lty=3,col="#1A7C4F"); text(120,0.03,"X_eq (equilibrio)",col="#1A7C4F",cex=0.7)
hist(r1$tv,breaks=40,col="#FADBD8",border="white",freq=FALSE,xlab="Edad de cese (años)",
     main=sprintf("t_cese -- mediana %.0fa (Γ_ext)",median(r1$tv)), xlim=c(0,300))
abline(v=median(r1$tv),lwd=2,col="#922B21")
plot(c(60,150),c(60,150),type="n",log="xy",xlab="Target (half-life, log)",ylab="Mediana sim (log, promedio 5 semillas)",
     main=sprintf("Paridad -- error %.1f%% (promedio 5 semillas)",err))
abline(0,1,lty=2,col="gray40",lwd=2); lines(c(60,150),c(60,150)*1.1,lty=3,col="gray60"); lines(c(60,150),c(60,150)*0.9,lty=3,col="gray60")
points(R_fran$t_cese_target,prom,pch=19,col="#C0392B",cex=2); text(R_fran$t_cese_target,prom,"S.franciscanus",pos=4,cex=0.85,col="#C0392B")
legend("topleft",legend=c("y=x","banda +-10%"),lty=c(2,3),col=c("gray40","gray60"),lwd=c(2,1),cex=0.7,bty="n")
dev.off()
cat("\nGrafica: TMCSA_Sfranciscanus_canonico_v3.png\n")
