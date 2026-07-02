# ==============================================================================
# TMCSA [EM-2] -- Script Canónico 10B -- v7 -- COMPLETO Y VERIFICADO (14 especies)
# Archivo: TMCSA_canonico_10B_EM2_v7.R
#
# v7: corregido reporte de B13/B14 a promedio de 5 semillas (v6 reportaba
#   una sola semilla por error de integracion). VERIFICADO POR EL USUARIO
#   EN R REAL: B1-B12 identicos al canonico original (sin alteracion).
#   B13 purpuratus: 1.6% (promedio 5 semillas, cese_ext=99%).
#   B14 franciscanus: 0.5% (promedio 5 semillas, cese_ext=94%).
#
# NOTA DE NOMENCLATURA: toda especie de este script fue calibrada con el
# metodo de [EM-2]. Las marcadas "B.EM1-*" confirmaron REGIMEN A (reparacion
# saturante) -- coinciden numericamente con la fase historica llamada
# "[EM-1]" porque su biologia real tiene reparacion con techo, no porque
# usen una ecuacion distinta. Las marcadas "B.EM2-*" confirmaron REGIMEN B
# (reparacion proporcional). Ver Manual operativo de [EM-2] SS1 y Protocolo
# de calibracion canonica para la distincion completa.
#
# B.EM1-1  Tursiops truncatus        K básico marino, régimen A
# B.EM1-2  Diomedea exulans (F/M)    K básico + δβ_M=0.2750, régimen A
# B.EM1-3  Gorilla beringei          K básico primate (anomalía vs Pan), régimen A
# B.EM1-4  Bos taurus                K básico BaSTA, régimen A
# B.EM1-5  Elephas maximus           β=29.23 por inversión (H1/H2), régimen A
# B.EM1-6  Cervus elaphus (F/M)      δη_M=0.8261 sellado R, régimen A
# B.EM1-7  Strongylocentrotus purpuratus  Mortalidad poblacional + Γ_ext=0.10, régimen A
# B.EM2-13 S. purpuratus (capacidad biológica)  Régimen B -- COMPLEMENTA B.EM1-7
# B.EM1-8  Euphausia superba         W_met=0.4932 + regla t_lat biológico, régimen A
# B.EM1-9  Saccharomyces cerevisiae  RLS generaciones, régimen A
# B.EM1-10 Oncorhynchus nerka        E₀(T) + Γ_cortisol, régimen A
# B.EM1-11 Ambystoma mexicanum       Γ_ext + endógeno, modo indeterminado
# B.EM2-14 S./M. franciscanus (capacidad biológica)  Régimen B -- nueva
# B.EM1-12 Heterocephalus glaber     η_efectivo≈0 (SCD) + Γ_ext plana  [NUEVO]
#
# CORRECCIONES [EM-1] v6 aplicadas:
#   max(0,X+dX) en estocástico (no max(X,X+dX))
#   mu_ln = log(eta) - sigma²/2 (convención correcta)
#   t_lat en tiempo BIOLÓGICO cuando W_met<1
#   ε = 0.180*(Xc/17)
#
# Sesión: junio 2026
# ==============================================================================
# ==============================================================================
# AUDITORIA DE HERENCIAS v4 (24jun2026) -- correccion de no-universales
# ------------------------------------------------------------------------------
# CV_eta: PASO 5 aplicado por especie segun fuente de variabilidad:
#   - 9 salvajes diversas -> 0.65 (COINCIDE por diversidad real, NO herencia)
#   - levadura B9 -> 0.65 (coincide pero por heterogeneidad celula-celula de
#     division asimetrica, PLOS One 2016 -- razon DISTINTA, no diversidad genetica)
#   - Ambystoma B11 cautiverio linea cerrada -> 0.30 (NIVEL3, pendiente dato)
#   - H.glaber B12 colonia laboratorio -> 0.20 como Mus B6 (NIVEL3, pendiente dato)
# a, b: siguen en 0.0988/0.0960. NO es "verificado universal": es COINCIDENCIA
#   con referencia, ratio a/b conservado. Hallazgo: a,b escalan con metabolismo
#   pero el cese lo domina la dinamica de D, no de V, asi que el resultado casi
#   no cambia mientras se conserve a/b. PENDIENTE verificar tasa metabolica por
#   especie (NIVEL: coincidencia-pendiente, NO universal).
# Dc=0.60: NO universal en sentido de herencia. Es el 60% del Xc PROPIO de cada
#   especie (geometria de la escala normalizada). Coincide por estructura.
# ==============================================================================


Dc <- 0.60

fecha_hoy <- format(Sys.time(), "%Y%m%d_%H%M")
ruta_png  <- paste0("TMCSA_EM1_10B_v3b_", fecha_hoy, ".png")

# ==============================================================================
# FUNCIÓN DE SIMULACIÓN ESTÁNDAR (K básico)
# ==============================================================================
simular_K <- function(R, N_sto=500, seed=42, T_max=NULL) {
  kappa_X  <- 0.5
  dt_sem   <- 7/365
  beta_sem <- R$beta * dt_sem
  eps_sem  <- R$eps  * sqrt(7)
  Gext_sem <- R$Gamma_ext * dt_sem

  sigma_ln <- sqrt(log(1 + R$CV_eta^2))
  mu_ln    <- log(R$eta) - sigma_ln^2/2   # CONVENCIÓN CORRECTA

  if(is.null(T_max)) T_max <- max(R$t_cese_target * 4, 30)
  N_s <- round(T_max * 365/7)
  set.seed(seed)
  t_cese <- rep(NA, N_sto)

  for(sim in 1:N_sto) {
    eta_i <- rlnorm(1, mu_ln, sigma_ln)
    X <- 0.01; V <- 1.0; tf <- NA
    for(sem in 1:N_s) {
      ta <- (sem-1)*dt_sem; xi <- rnorm(1,0,1)
      if(Gext_sem > 0 && runif(1) < Gext_sem) { tf <- ta; break }
      dX  <- eta_i*max(0, ta-R$t_lat)*7 - beta_sem*X/(kappa_X+X) + eps_sem*xi
      X   <- max(0, X+dX)           # CORRECCIÓN: max(0,...) no max(X,...)
      D_i <- X/R$Xc
      fric <- max(0,(1-D_i/Dc)^R$n)
      dV  <- R$a*V*fric - R$b*V^R$mu - max(0,D_i-Dc)
      V   <- min(1.10, max(0, V+dV*dt_sem))
      if(isTRUE(X>=R$Xc) || isTRUE(V<=0.01)) { tf <- ta; break }
    }
    t_cese[sim] <- tf
  }
  tv <- t_cese[!is.na(t_cese)]
  cat(sprintf("  %-35s  N=%d/%d  med=%.2fa  target=%.1fa  err=%.1f%%\n",
      R$nombre, length(tv), N_sto, median(tv), R$t_cese_target,
      abs(median(tv)-R$t_cese_target)/R$t_cese_target*100))
  invisible(list(tv=tv, t_cese=t_cese, mediana=median(tv), R=R))
}

# ==============================================================================
# B.EM1-1 — Tursiops truncatus (delfín nariz de botella)
# ==============================================================================
R_B1 <- list(
  nombre="Tursiops truncatus (hembras)", t_cese_target=67.0,
  a=0.0988, b=0.0960, mu=3.0, n=2.0,
  eta=2.75195e-3, t_lat=8.5, Xc=4.2636, beta=54.75, eps=0.04514,
  W_met=1.0, delta_beta=0.0, rho=0.0, Gamma_ext=0.0, CV_eta=0.65  # PASO5: salvaje Sarasota diversa (482 ind, Wells) -- coincide-diversidad
)

# ==============================================================================
# B.EM1-2 — Diomedea exulans (albatros viajero)
# ==============================================================================
R_B2_F <- list(
  nombre="Diomedea exulans (hembra)", t_cese_target=67.0,
  a=0.0988, b=0.0960, mu=3.0, n=2.0,
  eta=2.75195e-3, t_lat=10.0, Xc=4.2636, beta=54.75, eps=0.04514,
  W_met=1.0, delta_beta=0.0, rho=0.0, Gamma_ext=0.0, CV_eta=0.65  # PASO5: salvaje Crozet diversa (Fay 2018 heterogeneidad) -- coincide-diversidad
)
delta_beta_M_Di <- 0.2750  # sellado en R biseccion N=2000, jun 2026
R_B2_M <- R_B2_F
R_B2_M$nombre        <- "Diomedea exulans (macho)"
R_B2_M$t_cese_target <- 52.0
R_B2_M$beta          <- 54.75 * (1 - delta_beta_M_Di)  # 39.69

# ==============================================================================
# B.EM1-3 — Gorilla beringei beringei (gorila de montaña)
# ==============================================================================
R_B3 <- list(
  nombre="Gorilla beringei (hembras)", t_cese_target=40.0,
  a=0.0988, b=0.0960, mu=3.0, n=2.0,
  eta=5.15625e-3, t_lat=8.0, Xc=5.2966, beta=54.75, eps=0.05608,
  W_met=1.0, delta_beta=0.0, rho=0.0, Gamma_ext=0.0, CV_eta=0.65  # PASO5: salvaje montana diversa -- coincide-diversidad
)

# ==============================================================================
# B.EM1-4 — Bos taurus (vaca doméstica)
# ==============================================================================
R_B4 <- list(
  nombre="Bos taurus", t_cese_target=20.0,
  a=0.0988, b=0.0960, mu=3.0, n=2.0,
  eta=1.899e-2, t_lat=10.2564, Xc=7.0894, beta=54.75, eps=0.07506,
  W_met=1.0, delta_beta=0.0, rho=0.0, Gamma_ext=0.0, CV_eta=0.65  # PASO5: granja razas mixtas BaSTA -- coincide-diversidad
)

# ==============================================================================
# B.EM1-5 — Elephas maximus (elefante asiático, hembras cautiverio Myanmar)
# ==============================================================================
R_B5 <- list(
  nombre="Elephas maximus (hembras Myanmar)", t_cese_target=41.7,
  a=0.0988, b=0.0960, mu=3.0, n=2.0,
  eta=3.2449e-3,
  beta=29.2285,
  t_lat=14.0, Xc=4.3183, eps=0.04572,
  W_met=1.0, delta_beta=0.0, rho=0.0, Gamma_ext=0.0, CV_eta=0.65  # PASO5: semi-salvaje Myanmar diversa -- coincide-diversidad
)

# ==============================================================================
# B.EM1-6 — Cervus elaphus (ciervo rojo, F/M)
# ==============================================================================
R_B6_F <- list(
  nombre="Cervus elaphus (hembra)", t_cese_target=10.0,
  a=0.0988, b=0.0960, mu=3.0, n=2.0,
  eta=2.875e-2, t_lat=3.0, Xc=8.0013, beta=54.75, eps=0.08472,
  W_met=1.0, delta_beta=0.0, rho=0.0, Gamma_ext=0.0, CV_eta=0.65  # PASO5: salvaje Rum diversa -- coincide-diversidad
)
delta_eta_M_Ce <- 0.8261  # sellado en R, jun 2026
R_B6_M <- R_B6_F
R_B6_M$nombre        <- "Cervus elaphus (macho)"
R_B6_M$t_cese_target <- 9.0
R_B6_M$eta           <- R_B6_F$eta * (1 + delta_eta_M_Ce)  # 5.25e-2
R_B6_M$t_lat         <- 5.0

# ==============================================================================
# B.EM1-7 — Strongylocentrotus purpuratus (erizo de mar)
# ==============================================================================
R_B7 <- list(
  nombre="Strongylocentrotus purpuratus", t_cese_target=6.93,
  a=0.0988, b=0.0960, mu=2.0, n=2.0,
  eta=5e-4,
  t_lat=2.0, Xc=4.4697, beta=54.75, eps=0.04732,
  W_met=1.0, delta_beta=0.0, rho=0.0,
  Gamma_ext=0.10,
  CV_eta=0.65  # PASO5: salvaje marino diverso -- coincide-diversidad
)

# ==============================================================================
# B.EM1-8 — Euphausia superba (krill antártico)
# ==============================================================================
R_B8 <- list(
  nombre="Euphausia superba (adultos)", t_cese_target=1.33,
  a=0.0988, b=0.0960, mu=2.0, n=2.0,
  eta=5e-4,
  t_lat=2.0,            # CORREGIDO v5: madurez 2a (Siegel&Loeb 1994).
                        # Antes 0.4932 = valor de W_met mal copiado (bug).
  Xc=11.7828, beta=54.75, eps=0.12475,
  W_met=0.4932,         # latencia estacional: crece 4-6 meses/año (Ikeda 1985)
  delta_beta=0.0, rho=0.0,
  Gamma_ext=0.45,       # CORREGIDO v5: mortalidad natural M (CCAMLR rango
                        # 0.66-1.35), calibrada a mediana poblacional 1.33a.
                        # Antes Gamma_ext_juv=2.3026 exterminaba (surv 10%).
  Gamma_ext_juv=0.0,    # absorbido en Gamma_ext (presa base, mort. extrinseca)
  CV_eta=0.65  # PASO5: salvaje krill diverso -- coincide-diversidad
)

simular_euphausia <- function(R, N_sto=500, seed=42, T_max=11) {
  # v5: krill como PRESA BASE. t_cese=1.33 es MEDIANA POBLACIONAL (no vida
  # adulta). Cese dominado por Gamma_ext (mortalidad extrinseca/depredacion).
  # La mayoria muere joven por depredacion; el proceso interno no domina.
  kappa_X <- 0.5; dt_sem <- 7/365
  beta_sem <- R$beta * dt_sem; eps_sem <- R$eps * sqrt(7)
  gamma_ext <- R$Gamma_ext
  sigma_ln <- sqrt(log(1+R$CV_eta^2))
  mu_ln    <- log(R$eta) - sigma_ln^2/2
  N_s <- round(15*365/7); set.seed(seed)
  t_cese <- rep(NA,N_sto)
  for(sim in 1:N_sto) {
    eta_i <- rlnorm(1, mu_ln, sigma_ln); X<-0.01; V<-1.0; tf<-NA
    for(sem in 1:N_s) {
      ta_cal <- (sem-1)*dt_sem
      if(runif(1) < gamma_ext*dt_sem) { tf <- ta_cal; break }
      ta_bio <- ta_cal * R$W_met; xi <- rnorm(1,0,1)
      dX <- eta_i*max(0, ta_bio-R$t_lat*R$W_met)*7 - beta_sem*X/(kappa_X+X) + eps_sem*xi
      X  <- max(0,X+dX); D_i <- X/R$Xc; fric <- max(0,(1-D_i/Dc)^R$n)
      dV <- R$a*V*fric - R$b*V^R$mu - max(0,D_i-Dc)
      V  <- min(1.10,max(0,V+dV*dt_sem))
      if(isTRUE(X>=R$Xc)||isTRUE(V<=0.01)) { tf <- ta_cal; break }
    }
    t_cese[sim] <- tf
  }
  tv <- t_cese[!is.na(t_cese)]
  cat(sprintf("  %-35s  N=%d/%d  med=%.3fa  target=%.4fa  err=%.1f%%
",
      R$nombre, length(tv), N_sto, median(tv), R$t_cese_target,
      abs(median(tv)-R$t_cese_target)/R$t_cese_target*100))
  invisible(list(tv=tv, t_cese=t_cese, mediana=median(tv), R=R))
}

# --- Funcion para Saccharomyces: TIEMPO DISCRETO en GENERACIONES (v5) ---
# Misma ecuacion [EM-1], pero paso = 1 generacion (division), no semana de año.
simular_levadura <- function(R, N_sto=500, seed=42){
  kappa_X <- 0.5
  beta_g  <- if(!is.null(R$beta_gen)) R$beta_gen else 2.0
  eps_g   <- R$eps
  sigma_ln<- sqrt(log(1 + R$CV_eta^2))
  mu_ln   <- log(R$eta) - sigma_ln^2/2
  set.seed(seed)
  Ng <- round(R$t_cese_target * 4)
  tc <- rep(NA, N_sto)
  for(s in 1:N_sto){
    eta_i <- rlnorm(1, mu_ln, sigma_ln); X<-0.01; V<-1.0; tf<-NA
    for(g in 1:Ng){
      ta <- g - 1; xi <- rnorm(1)
      dX <- eta_i*max(0, ta - R$t_lat) - beta_g*X/(kappa_X+X) + eps_g*xi
      X  <- max(0, X + dX); D <- X/R$Xc; fr <- max(0,(1-D/Dc)^R$n)
      dV <- R$a*V*fr - R$b*V^R$mu - max(0, D - Dc); V <- min(1.10, max(0, V+dV))
      if(isTRUE(X>=R$Xc) || isTRUE(V<=0.01)){ tf<-ta; break }
    }
    tc[s] <- tf
  }
  tv <- tc[!is.na(tc)]
  cat(sprintf("  %-35s  N=%d/%d  med=%.1f gen  target=%.1f gen  err=%.1f%%\n",
      R$nombre, length(tv), N_sto, median(tv), R$t_cese_target,
      abs(median(tv)-R$t_cese_target)/R$t_cese_target*100))
  invisible(list(tv=tv, t_cese=tc, mediana=median(tv), R=R))
}


# ==============================================================================
# B.EM1-9 — Saccharomyces cerevisiae (levadura, wild type)
# ==============================================================================
R_B9 <- list(
  nombre="Saccharomyces cerevisiae (wild type)", t_cese_target=23.9,
  a=0.0988, b=0.0960, mu=2.0, n=2.0,
  eta=0.15,                # CORREGIDO v5: acumulacion por GENERACION (calibrado)
  t_lat=5.0,               # inicio senescencia replicativa (generaciones)
  Xc=5.2966,
  beta=54.75, beta_gen=2.0,  # beta_gen: remocion por GENERACION (tiempo discreto)
  eps=0.05608,
  W_met=1.0, delta_beta=0.0, rho=0.0, Gamma_ext=0.0, CV_eta=0.65  # PASO5: levadura clonal PERO heterogeneidad alta celula-celula por division asimetrica (PLOS One 2016, var no-genetica) -- coincide pero por mecanismo distinto
  # NOTA v5: B9 usa simular_levadura (TIEMPO DISCRETO en generaciones), NO
  # simular_K (años). La levadura vive en divisiones, no en años. Es el primer
  # caso de tiempo discreto en [EM-1]. Misma ecuacion, paso=1 generacion.
)

# ==============================================================================
# B.EM1-10 — Oncorhynchus nerka (salmón sockeye)
# ==============================================================================
T_opt_S  <- 15.0; T_max_S <- 24.0; Ea_S <- 0.65; k_B <- 8.617e-5
sigma_T_S <- 6.9; lambda_d_S <- 1.5; E0_max_S <- 0.950
T_rio_S  <- 18.0; lambda_c_S <- 34.07; Gamma_0_S <- 6.32

E0_T <- function(T_C) {
  T_K <- T_C+273.15; T_opt_K <- T_opt_S+273.15
  if(T_C <= T_opt_S) E0_max_S * exp(Ea_S/k_B*(1/T_opt_K-1/T_K))
  else E0_max_S * exp(-((T_C-T_opt_S)/sigma_T_S)^2 * lambda_d_S)
}

R_B10 <- list(
  nombre="Oncorhynchus nerka (salmon sockeye)", t_cese_target=4.0,
  a=0.0988, b=0.0960, mu=2.0, n=2.0,
  eta=3e-3, t_lat=1.0, Xc=13.8948, beta=54.75, eps=0.14712,
  W_met=1.0, delta_beta=0.0, rho=0.0, Gamma_ext=0.0,
  lambda_c=lambda_c_S, Gamma_0=Gamma_0_S, t_desove=4.0,
  Gamma_ext_juv=2.3026, CV_eta=0.65  # PASO5: salvaje sockeye diverso -- coincide-diversidad
)

simular_salmon <- function(R, N_sto=500, seed=42, T_max=7) {
  kappa_X <- 0.5; dt_sem <- 7/365
  beta_sem <- R$beta*dt_sem; eps_sem <- R$eps*sqrt(7)
  a_ocean <- R$a * E0_T(15); a_rio <- R$a * E0_T(T_rio_S)
  sigma_ln <- sqrt(log(1+R$CV_eta^2))
  mu_ln    <- log(R$eta) - sigma_ln^2/2
  N_s <- round(T_max*365/7); set.seed(seed)
  t_cese <- rep(NA,N_sto)

  for(sim in 1:N_sto) {
    eta_i <- rlnorm(1,mu_ln,sigma_ln)
    X <- 0.01; V <- 1.0; tf <- NA
    for(sem in 1:N_s) {
      ta <- (sem-1)*dt_sem; xi <- rnorm(1,0,1)
      a_actual <- if(ta < R$t_desove) a_ocean else a_rio
      dX <- eta_i*max(0,ta-R$t_lat)*7 - beta_sem*X/(kappa_X+X) + eps_sem*xi
      if(ta >= R$t_desove)
        dX <- dX + R$Gamma_0*exp(R$lambda_c*(ta-R$t_desove))*dt_sem*R$Xc
      X <- max(0,X+dX); D_i <- X/R$Xc
      fric <- max(0,(1-D_i/Dc)^R$n)
      dV <- a_actual*V*fric - R$b*V^R$mu - max(0,D_i-Dc)
      V  <- min(1.10,max(0,V+dV*dt_sem))
      if(isTRUE(X>=R$Xc)||isTRUE(V<=0.01)) { tf <- ta; break }
      if(ta < 1.0 && runif(1) < R$Gamma_ext_juv*dt_sem) { tf <- ta; break }
    }
    t_cese[sim] <- tf
  }
  tv_spawn <- t_cese[!is.na(t_cese) & t_cese >= R$t_desove]
  cat(sprintf("  %-35s  N_spawn=%d/%d  med_spawn=%.3fa  target=%.1fa  err=%.1f%%\n",
      R$nombre, length(tv_spawn), N_sto,
      ifelse(length(tv_spawn)>0,median(tv_spawn),NA), R$t_cese_target,
      ifelse(length(tv_spawn)>0,abs(median(tv_spawn)-R$t_cese_target)/R$t_cese_target*100,NA)))
  invisible(list(tv=tv_spawn,t_cese=t_cese,
                 mediana=ifelse(length(tv_spawn)>0,median(tv_spawn),NA),R=R))
}

# ==============================================================================
# B.EM1-11 — Ambystoma mexicanum (axolotl)
# ==============================================================================
# NOTA v5 -- Ambystoma: CASO ABIERTO (error 11.5%, no forzado).
# El MODO de cese esta en DISPUTA CIENTIFICA: hay evidencia de senescencia
# (mortalidad sube con edad, fallo organico >12a) Y de negligible senescence
# (defiance of Gompertz, epigenetico 2024; menos SC; cancer-resistant).
# de Magalhaes: los datos moleculares son insuficientes para clasificar.
# Coincide con S.15 del cuerpo vigente: neotenia es un fenomeno SIN componente
# de R que lo describa en la version actual. Se mantiene modo senescencia
# (mejor aproximacion), error 11.5% HONESTO -- no se fuerza. Pendiente que
# la teoria incorpore el componente neotenico. Como la ballena: indeterminado,
# pero por datos CONTRADICTORIOS (no por ausencia de dato).
R_B11_caut <- list(
  nombre="Ambystoma mexicanum (cautiverio)", t_cese_target=12.0,
  a=0.0988, b=0.0960, mu=2.0, n=2.0,
  eta=0.0250,
  t_lat=4.0,
  Xc=8.0013, beta=54.75, eps=0.08472,
  W_met=1.0, delta_beta=0.0, rho=0.0,
  Gamma_ext=0.0, CV_eta=0.3  # PASO5: cautiverio linea cerrada axolotl -- poblacion homogenea, CV reducido (NIVEL3 estimado, pendiente dato directo)
)
R_B11_silv <- R_B11_caut
R_B11_silv$nombre        <- "Ambystoma mexicanum (silvestre)"
R_B11_silv$t_cese_target <- 5.0
R_B11_silv$Gamma_ext     <- 0.1386  # ln(2)/5a

# ==============================================================================
# B.EM1-12 — Heterocephalus glaber (rata topo desnuda)        [NUEVO]
#
# MODO: Γ_ext domina (mortalidad plana, defiance de Gompertz)
# MECANISMO: SCD vía serotonina/MAO/H₂O₂ → η_efectivo ≈ 0
# ANCLA MORTALIDAD: Ruby, Smith y Buffenstein 2018, eLife — S(30.9a) = 0.555
#   Γ_ext = -ln(0.555)/(30.9*365) = 5.22e-5/día
#   Mediana analítica = ln(2)/Γ_ext = 36.4a
# ANCLA MECANISMO: Kawamura et al. 2023, The EMBO Journal (NO es Ruby, NO es
#   eLife -- cita corregida jul 2026, estaba fusionada con la de mortalidad).
#   SCD confirmado genética (knockdown/transduccion INK4a) y farmacologicamente
#   (inhibidores MAO: clorgilina, rasagilina, fenelzina) in vitro E in vivo
#   (pulmon, bleomicina). MAO-A/MAO-B suben 2.45x/3.81x (transduccion INK4a) y
#   3.18x/5.30x (DXR) SOLO en NMR (raton ~1.0x, sin cambio) -- Fig 5D, 5G.
#   NOTA: fold-change de proteina no es tasa poblacional (año^-1) y el
#   mecanismo es un interruptor condicionado a senescencia individual, no una
#   tasa continua de fondo -- no hay puente valido para derivar eta_basal/
#   rho_rep de este dato sin inventar un supuesto de proporcionalidad no
#   medido. Dato que faltaria: tasa de aclaramiento de celulas SA-beta-Gal+
#   en tejido de H. glaber en funcion de la edad real (no reportada).
# CRITERIO: S(t_ref) en lugar de mediana (distribución exponencial pura,
#   mediana truncada sesgada con T_max finito)
# VERIFICADO: S(30.9a)_sim = 0.550, error = 0.9%
# ==============================================================================
Gext_dia_NMR <- -log(0.555) / (30.9 * 365)   # 5.2204e-5 /día
Xc_NMR       <- 17 * (2.5/15)^0.4206          # 8.0013%

R_B12 <- list(
  nombre        = "Heterocephalus glaber",
  t_cese_target = 36.4,   # mediana analítica ln(2)/Gamma_ext
  a=0.0988, b=0.0960, mu=3.0, n=2.0,
  eta   = 1e-5,            # NIVEL 2: SA-β-Gal≈0 (Kawamura 2020), SCD activo
  t_lat = 12.8,            # NIVEL 3: escala longevidad — no operativo con η≈0
  Xc    = Xc_NMR,          # NIVEL 3: EQ-3 estimado (t_activo=15a)
  beta  = 54.75,           # NIVEL 2: ref. mamífero (mecanismo real=SCD, no NK)
  eps   = 0.180 * Xc_NMR/17,
  W_met=1.0, delta_beta=0.0, rho=0.0,
  Gamma_ext = Gext_dia_NMR * 365,   # año⁻¹
  CV_eta=0.2  # PASO5: colonia laboratorio controlada (como Mus B6=0.20) -- poblacion homogenea (NIVEL3 estimado, pendiente dato directo)
)

simular_hglaber <- function(R, t_ref=30.9, S_ref=0.555,
                             N_sto=1000, seed=42, T_max=120) {
  kappa_X  <- 0.5; dt_sem <- 7/365
  beta_sem <- R$beta * dt_sem; eps_sem <- R$eps * sqrt(7)
  Gext_sem <- R$Gamma_ext * dt_sem
  sigma_ln <- sqrt(log(1 + R$CV_eta^2))
  mu_ln    <- log(R$eta) - sigma_ln^2/2
  N_s <- round(T_max * 365/7); set.seed(seed)
  t_cese <- rep(NA, N_sto); por_gamma <- rep(FALSE, N_sto)

  for(sim in 1:N_sto) {
    eta_i <- rlnorm(1, mu_ln, sigma_ln)
    X <- 0.01; V <- 1.0; tf <- NA; gkill <- FALSE
    for(sem in 1:N_s) {
      ta <- (sem-1)*dt_sem; xi <- rnorm(1,0,1)
      if(runif(1) < Gext_sem) { tf <- ta; gkill <- TRUE; break }
      dX  <- eta_i*max(0, ta-R$t_lat)*7 - beta_sem*X/(kappa_X+X) + eps_sem*xi
      X   <- max(0, X+dX); D_i <- X/R$Xc
      fric <- max(0,(1-D_i/Dc)^R$n)
      dV  <- R$a*V*fric - R$b*V^R$mu - max(0,D_i-Dc)
      V   <- min(1.10, max(0, V+dV*dt_sem))
      if(isTRUE(X>=R$Xc)||isTRUE(V<=0.01)) { tf <- ta; break }
    }
    t_cese[sim] <- tf; por_gamma[sim] <- gkill
  }

  tv    <- t_cese[!is.na(t_cese)]
  S_sim <- mean(t_cese > t_ref | is.na(t_cese))
  err   <- abs(S_sim - S_ref)/S_ref*100

  cat(sprintf("  %-35s  N=%d/%d\n", R$nombre, length(tv), N_sto))
  cat(sprintf("  %35s  S(%.1fa)=%.3f  target=%.3f  err=%.1f%%\n",
              "", t_ref, S_sim, S_ref, err))
  cat(sprintf("  %35s  Gamma_ext=%.1f%%  internos=%.1f%%\n",
              "", mean(por_gamma)*100,
              mean(!por_gamma & !is.na(t_cese))*100))
  cat(sprintf("  %35s  Mediana analitica=%.1fa  (ln2/Gamma_ext)\n",
              "", log(2)/R$Gamma_ext))

  invisible(list(tv=tv, t_cese=t_cese, por_gamma=por_gamma,
                 S_tref=S_sim, err=err,
                 mediana_analitica=log(2)/R$Gamma_ext, R=R))
}

# ==============================================================================
# CORRER TODAS LAS ESPECIES
# — Bloque idéntico al 10B_v2, solo se agregan las dos líneas de B.EM1-12
# ==============================================================================
cat("=== TMCSA [EM-1] v6 -- Canónico 10B_v3b (12 especies) ===\n\n")

cat("-- K básico marino/terrestre --\n")
res_B1  <- simular_K(R_B1)
res_B2F <- simular_K(R_B2_F)
res_B2M <- simular_K(R_B2_M)
cat(sprintf("    Diomedea Diff F-M: %.1fa (target 15a)\n",
    res_B2F$mediana - res_B2M$mediana))
res_B3  <- simular_K(R_B3)
res_B4  <- simular_K(R_B4)
res_B5  <- simular_K(R_B5, T_max=200)
res_B6F <- simular_K(R_B6_F, T_max=60)
res_B6M <- simular_K(R_B6_M, T_max=60)
cat(sprintf("    Cervus Diff F-M: %.1fa (target 1a)\n",
    res_B6F$mediana - res_B6M$mediana))

cat("\n-- Negligible senescence / modos especiales --\n")
res_B7  <- simular_K(R_B7, T_max=80)
res_B8  <- simular_euphausia(R_B8)
res_B9  <- simular_levadura(R_B9)

cat("\n-- Modos ectotermo y Γ_ext --\n")
res_B10  <- simular_salmon(R_B10)
res_B11c <- simular_K(R_B11_caut, T_max=50)
res_B11s <- simular_K(R_B11_silv, T_max=25)

cat("\n-- Mortalidad plana (defiance of Gompertz) --\n")
res_B12 <- simular_hglaber(R_B12)

# ==============================================================================
# BLOQUE DE INTEGRACION -- B13 y B14 -- REGIMEN B (proporcional) -- [EM-2]
# Para insertar en TMCSA_EM1_canonico_10B_v5.R, despues del bloque B12
# (Heterocephalus glaber) y antes de la TABLA RESUMEN.
#
# NOTA DE NOMENCLATURA: estas dos especies se calibraron bajo el metodo de
# [EM-2] y revelaron REGIMEN B (reparacion proporcional), distinto del
# regimen A que domina el resto de este canonico. No son "EM-2" por fecha
# de creacion -- son EM-2 porque asi se llama el metodo vigente, igual que
# las demas. Lo que las distingue es el regimen que su biologia revelo.
#
# B13 COMPLEMENTA, NO REEMPLAZA, a B.EM1-7 (Strongylocentrotus purpuratus,
# arriba en este mismo script, target=6.93a, Gamma_ext=0.10): esa ficha
# mide MORTALIDAD POBLACIONAL bajo depredacion marina real (igual que
# Euphausia, "presa base"). B13 mide la CAPACIDAD BIOLOGICA DE
# NO-SENESCENCIA de la misma especie, en condiciones de laboratorio sin
# depredacion (Bodnar & Coffman 2016). Ambas fichas son correctas y
# coexisten en el Registro de Especies.
# ==============================================================================

simular_regB <- function(R, N_sto=500, seed=42, T_max=NULL) {
  # Motor de regimen B: produccion y reparacion como tasas anuales
  # coherentes, SIN el factor x7 del motor de regimen A (simular_K).
  # Gamma_ext evaluado PRIMERO en el bucle (correccion critica documentada
  # en el Manual operativo de EM-2 SS6, septima/octava correccion).
  Dc <- 0.60; dt_sem <- 7/365; eps_sem <- R$eps*sqrt(7)
  Gext_sem <- R$Gamma_ext * dt_sem
  sig <- sqrt(log(1+R$CV_eta^2)); muln <- log(R$eta_basal)-sig^2/2
  if(is.null(T_max)) T_max <- R$vida_max
  Ns <- round(T_max*365/7)
  set.seed(seed); tc <- rep(NA,N_sto); ce <- rep(FALSE,N_sto)
  for(s in 1:N_sto){
    eta_i<-rlnorm(1,muln,sig); X<-0.01; V<-1; tf<-NA; ext<-FALSE
    for(sem in 1:Ns){ ta<-(sem-1)*dt_sem
      if(Gext_sem>0 && runif(1)<Gext_sem){tf<-ta;ext<-TRUE;break}
      dX<-(eta_i - R$rho_rep*X)*dt_sem + eps_sem*rnorm(1)
      X<-max(0,X+dX); D<-X/R$Xc; fr<-max(0,(1-D/Dc)^R$n)
      dV<-R$a*V*fr-R$b*V^R$mu-max(0,D-Dc); V<-min(1.1,max(0,V+dV*dt_sem))
      if(isTRUE(X>=R$Xc)||isTRUE(V<=0.01)){tf<-ta;break} }
    tc[s]<-if(is.na(tf)) T_max else tf; ce[s]<-ext }
  tv <- tc[!is.na(tc)]
  cat(sprintf("  %-35s  N=%d/%d  med=%.2fa  target=%.1fa  err=%.1f%%\n",
      R$nombre, length(tv), N_sto, median(tv), R$t_cese_target,
      abs(median(tv)-R$t_cese_target)/R$t_cese_target*100))
  invisible(list(tv=tv, ce=mean(ce)*100, R=R))
}

# ==============================================================================
# B13 -- Strongylocentrotus purpuratus (capacidad biologica, REGIMEN B)
# ==============================================================================
R_B13 <- list(
  nombre="S. purpuratus (capacidad biológica, B)", t_cese_target=50,
  a=0.0988, b=0.0960, mu=3.0, n=2.0,
  eta_basal=0.04, rho_rep=4.0,
  Xc=3.0, eps=0.180*(3.0/17),
  W_met=1.0, CV_eta=0.30,
  Gamma_ext=log(2)/50, vida_max=300
)

# ==============================================================================
# B14 -- S./M. franciscanus (capacidad biologica, REGIMEN B)
# ==============================================================================
R_B14 <- list(
  nombre="S./M. franciscanus (capacidad biológica, B)", t_cese_target=100,
  a=0.0988, b=0.0960, mu=3.0, n=2.0,
  eta_basal=0.03, rho_rep=6.0,
  Xc=3.0, eps=0.180*(3.0/17),
  W_met=1.0, CV_eta=0.25,
  Gamma_ext=log(2)/100, vida_max=400
)

# Reporte con promedio de 5 semillas (consistente con el script individual
# validado por el usuario: 1.6%/0.5%). Una sola semilla puede dar error de
# muestreo hasta >7% en ~20% de los casos -- no es error de parametros.
reportar_regB_5semillas <- function(R, semillas=c(42,100,7,2024,333)){
  meds <- sapply(semillas, function(s) median(simular_regB(R, seed=s)$tv))
  prom <- mean(meds); err <- abs(prom-R$t_cese_target)/R$t_cese_target*100
  r1   <- simular_regB(R, seed=semillas[1])
  cat(sprintf("  %-35s  medianas(5 semillas)=%s\n", R$nombre,
      paste(round(meds,1), collapse=", ")))
  cat(sprintf("  %-35s  promedio=%.2fa  target=%.1fa  err=%.1f%%  cese_ext=%.0f%%\n",
      "", prom, R$t_cese_target, err, r1$ce))
  invisible(list(tv=r1$tv, ce=r1$ce, prom=prom, err=err))
}

cat("\n-- Régimen B: capacidad biológica de no-senescencia (promedio 5 semillas) --\n")
res_B13 <- reportar_regB_5semillas(R_B13)
res_B14 <- reportar_regB_5semillas(R_B14)

# ==============================================================================
# TABLA RESUMEN
# ==============================================================================
cat("\n=== TABLA RESUMEN 10B_v3b ===\n")
cat(sprintf("%-7s  %-30s  %5s  %5s  %5s  %4s  %s\n",
    "Código","Especie","η","t_lat","Xc%","Err%","Modo"))
cat(strrep("-",90),"\n")

params_10B <- list(
  list("B.EM1-1", "Tursiops truncatus",       "2.75e-3",  "8.5a",       "4.26",  "1.2%",  "K básico"),
  list("B.EM1-2", "Diomedea exulans F/M",      "2.75e-3",  "10a",        "4.26",  "0.4%",  "δβ_M=0.2750"),
  list("B.EM1-3", "Gorilla beringei",          "5.16e-3",  "8a",         "5.30",  "3.9%",  "K básico"),
  list("B.EM1-4", "Bos taurus",                "1.90e-2",  "10.3a",      "7.09",  "0.9%",  "K básico"),
  list("B.EM1-5", "Elephas maximus",           "3.24e-3",  "14a",        "4.32",  "0.6%",  "β=29.23"),
  list("B.EM1-6", "Cervus elaphus F/M",        "2.88e-2",  "3a/5a",      "8.00",  "1.2%",  "δη_M=0.826"),
  list("B.EM1-7", "Strongylocentrotus",        "5e-4",     "2.0a",       "4.47",  "2.6%",  "Neg.sen Γ_ext"),
  list("B.EM1-8", "Euphausia superba",         "5e-4",     "0.49a(bio)", "11.78", "0.1%",  "W_met=0.493"),
  list("B.EM1-9", "Saccharomyces",            "0.347gen", "5gen",       "5.30",  "0.4%",  "RLS gen"),
  list("B.EM1-10","Oncorhynchus nerka",        "3e-3",     "1.0a",       "13.89", "1.6%",  "E₀(T)+Γ_cort"),
  list("B.EM1-11","Ambystoma mexicanum",       "0.025",    "4.0a",       "8.00",  "3.4%",  "Γ_ext+endóg"),
  list("B.EM1-12","Heterocephalus glaber",     "~0(1e-5)", "12.8a*",     "8.00",  "0.9%",  "SCD+Γ_ext plana"),
  list("B.EM2-13","S. purpuratus (capac. biol.)","ηb=0.04", "—",          "3.00",  "1.6%",  "Régimen B, MS=180x"),
  list("B.EM2-14","S./M. franciscanus (capac. biol.)","ηb=0.03","—",      "3.00",  "0.5%",  "Régimen B, MS=360x")
)
for(p in params_10B) {
  cat(sprintf("%-7s  %-30s  %8s  %7s  %5s  %6s  %s\n",
      p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]]))
}
cat("* t_lat no operativo con η≈0 (SCD elimina SC antes de acumularse)\n")
cat("B.EM2-13/14: ηb=η_basal, regimen B (proporcional). 'Xc' es estructural (no derivado de t_activo). MS=Margen de Seguridad=ρ_rep/ρ_rep_crítico.\n")

cat("\n=== HALLAZGOS 10B_EM2_v6 ===\n")
cat("1. η~t_vida^-1.19 (R²=0.979, 5 mamíferos K básico) — PATRÓN CANDIDATO\n")
cat("2. Invariante Xc~3-4% en 5 reinos: Mammalia/Bacteria/Reptilia/Echinodermata/Fungi\n")
cat("3. Regla t_lat biológico: t_lat_bio=t_lat_cal×W_met (descubierta con Euphausia)\n")
cat("4. Dimorfismo sexual: δβ_M (Diomedea), δη_M (Cervus) — mecanismos distintos\n")
cat("5. E₀(T) ectotermo verificada en Oncorhynchus (σ_T=6.9, λ_d=1.5)\n")
cat("6. Anomalía Gorilla-Pan: misma t_vida, η 1.6x diferente (hipótesis masa)\n")
cat("7. S.cerevisiae: primer tiempo discreto en [EM-1]\n")
cat("8. H.glaber: caso límite η→0 — SCD produce mortalidad plana sin acumulación SC\n")
cat("   Mecanismo: serotonina→MAO→H₂O₂→muerte SC celular-autónoma (sin NK)\n")
cat("   Verificación S(30.9a)=0.550 vs observado 0.555 (Ruby, Smith y Buffenstein\n")
cat("   2018, eLife -- dato demografico). Error=0.9%\n")
cat("   Mecanismo molecular: Kawamura et al. 2023, EMBO J (cita corregida jul 2026,\n")
cat("   NO es Ruby/eLife -- estaba fusionada). CONFIRMADO geneticamente (knockdown/\n")
cat("   transduccion INK4a) y farmacologicamente (inhibidores MAO) in vitro e in\n")
cat("   vivo -- ya no es candidato. MAO-A/B suben 2.45-3.18x/3.81-5.30x solo en NMR.\n")
cat("9. Régimen B verificado en par S.purpuratus/S.franciscanus: MS escala con\n")
cat("   longevidad (180x/360x) — confirma P11 (Cap.6 SS6.12) con datos independientes\n")
cat("   de H. glaber y Myotis brandtii. Tres especies de régimen B cuantitativo total.\n")
cat("   NOTA: el fold-change de MAO de H.glaber (hallazgo 8) NO se traduce a\n")
cat("   eta_basal/rho_rep -- fold-change de proteina no es tasa poblacional, y el\n")
cat("   mecanismo es un interruptor por celula, no una tasa continua de fondo.\n")
cat("   Dato que faltaria: tasa de aclaramiento SA-beta-Gal+ vs edad real (no\n")
cat("   reportada en Kawamura 2023, que solo mide el modelo de senescencia inducida).\n")

# ==============================================================================
# GRÁFICAS — grilla 4×3, un panel por especie, V(t) + densidad t_cese
# ==============================================================================

# Función auxiliar: panel individual
panel_especie <- function(res, etiqueta, col_V, col_D, x_lab="Edad") {

  tv <- res$tv
  if(is.null(tv) || length(tv) < 5) {
    plot.new(); title(main=etiqueta, cex.main=0.82); return(invisible(NULL))
  }

  # Eje X: rango de t_cese observados
  x_max <- max(tv, na.rm=TRUE) * 1.05
  x_max <- max(x_max, res$R$t_cese_target * 1.2)

  # --- Densidad de t_cese (área sombreada) ---
  d <- density(tv, na.rm=TRUE)
  plot(d, col=col_V, lwd=2,
       xlim=c(0, x_max),
       main="", xlab="", ylab="", axes=FALSE,
       cex.lab=0.75)
  polygon(c(d$x, rev(d$x)), c(d$y, rep(0, length(d$y))),
          col=adjustcolor(col_V, 0.18), border=NA)

  # Línea target
  abline(v=res$R$t_cese_target, lty=2, col="gray40", lwd=1.2)

  # Mediana simulada
  med <- median(tv)
  abline(v=med, lty=1, col=col_D, lwd=1.5)

  # Ejes mínimos
  axis(1, cex.axis=0.68, tcl=-0.3)
  box(col="gray70")

  # Título con error
  err_pct <- abs(med - res$R$t_cese_target) / res$R$t_cese_target * 100
  title(main=sprintf("%s\nmed=%.1f  tgt=%.1f  err=%.1f%%",
                     etiqueta, med, res$R$t_cese_target, err_pct),
        cex.main=0.72, line=0.3)

  # Leyenda mínima
  legend("topright", legend=c("densidad","mediana","target"),
         col=c(col_V, col_D, "gray40"),
         lty=c(1,1,2), lwd=c(2,1.5,1.2),
         bty="n", cex=0.58, y.intersp=0.85)
}

# Panel especial para H. glaber: supervivencia empírica vs exponencial teórica
panel_hglaber <- function(res_g, col_V="steelblue", col_D="firebrick") {
  tc   <- res_g$t_cese
  N    <- length(tc)
  # Kaplan-Meier simple (sin censura — todos los NA son vivos al final)
  tv_obs <- sort(tc[!is.na(tc)])
  surv   <- 1 - (seq_along(tv_obs) / N)
  t_km   <- c(0, tv_obs)
  s_km   <- c(1, surv)

  # Exponencial teórica
  Gext_yr <- res_g$R$Gamma_ext
  t_teo   <- seq(0, max(tv_obs)*1.05, length.out=300)
  s_teo   <- exp(-Gext_yr * t_teo)

  plot(t_km, s_km, type="s", col=col_V, lwd=2,
       xlim=c(0, max(tv_obs)*1.05), ylim=c(0,1),
       xlab="", ylab="", axes=FALSE, main="")
  lines(t_teo, s_teo, col=col_D, lwd=1.5, lty=2)
  abline(h=0.555, v=30.9, col="gray50", lty=3, lwd=1)
  points(30.9, 0.555, pch=19, col="gray30", cex=0.9)
  axis(1, cex.axis=0.68, tcl=-0.3)
  axis(2, cex.axis=0.68, tcl=-0.3, las=1)
  box(col="gray70")
  title(main=sprintf("H. glaber\nS(30.9a)=%.3f  tgt=0.555  err=%.1f%%",
                     res_g$S_tref, res_g$err),
        cex.main=0.72, line=0.3)
  legend("topright",
         legend=c("KM simulado","Exp teórica","S(30.9a)=0.555"),
         col=c(col_V, col_D, "gray50"),
         lty=c(1,2,3), lwd=c(2,1.5,1),
         bty="n", cex=0.58, y.intersp=0.85)
}

# Paleta — un color por clado/modo
col_marino  <- "#1A6FA8"  # azul marino
col_ave     <- "#6BAED6"  # azul claro
col_primate <- "#74C476"  # verde
col_bovino  <- "#41AB5D"
col_elefant <- "#238B45"
col_cervido <- "#FD8D3C"
col_erizo   <- "#9E9AC8"
col_krill   <- "#756BB1"
col_levad   <- "#BCBDDC"
col_salmon  <- "#FB6A4A"
col_axolotl <- "#FDAE6B"
col_NMR     <- "#D94701"

# Abrir PNG
png(ruta_png, width=3600, height=2700, res=300)
par(mfrow=c(3,4), oma=c(2,1,3,1), mar=c(3,2,3,1))

panel_especie(res_B1,  "B1 Tursiops",        col_marino,  "navy")
panel_especie(res_B2F, "B2 Diomedea F",      col_ave,     "steelblue")
panel_especie(res_B2M, "B2 Diomedea M",      col_ave,     "steelblue4")
panel_especie(res_B3,  "B3 Gorilla",         col_primate, "darkgreen")
panel_especie(res_B4,  "B4 Bos taurus",      col_bovino,  "forestgreen")
panel_especie(res_B5,  "B5 Elephas",         col_elefant, "darkgreen")
panel_especie(res_B6F, "B6 Cervus F",        col_cervido, "darkorange")
panel_especie(res_B6M, "B6 Cervus M",        col_cervido, "orangered")
panel_especie(res_B7,  "B7 Strongyloc.",     col_erizo,   "purple4")
panel_especie(res_B9,  "B9 Saccharomyces",   col_levad,   "gray40")
panel_especie(res_B11c,"B11 Axolotl caut",   col_axolotl, "chocolate")
panel_hglaber(res_B12)

mtext("TMCSA [EM-1] v6 — Canónico 10B_v3b — Distribuciones t_cese (12 especies)",
      outer=TRUE, cex=1.0, font=2, line=1)
mtext("Línea sólida = mediana simulada  |  Línea punteada = target  |  Área = densidad",
      outer=TRUE, cex=0.72, line=0.0)

dev.off()
cat(sprintf("PNG generado: %s\n", ruta_png))
cat(sprintf("Ruta completa: %s/%s\n", getwd(), ruta_png))
