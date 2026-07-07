# ==============================================================================
# TMCSA -- SENSIBILIDAD GLOBAL (SOBOL) -- Familia 2 (Gamma_ext dominante)
# Especie representativa: Heterocephalus glaber (R_B12, script canonico 10B v8)
# Archivo: TMCSA_SensibilidadGlobal_Sobol_Hglaber_v2.R
#
# CAMBIO v2 (jul 2026), tras hallazgo de la v1:
#   La v1 uso t_cese y T_activo como unicos observables, con N_ind=20. Los
#   ST resultantes quedaron pegados a ~1.0 para los 10 parametros -- no es
#   un hallazgo biologico plausible, es la firma de un estimador de Sobol
#   inestable cuando la señal real es pequeña frente al ruido de muestreo
#   (Gamma_ext fijo absorbe casi toda la varianza de t_cese en esta especie).
#
#   Una exploracion rapida en Python con un tercer observable -- la FRACCION
#   de individuos cuyo cese fue por causa externa (Gamma_ext) vs. interna --
#   mostro que esa fraccion NO es constante: media~0.935 pero rango real
#   0.80-1.00 al mover los parametros +-30%. Es decir, en los extremos del
#   rango explorado, una fraccion no trivial de individuos SI puede colapsar
#   por causa interna antes de que Gamma_ext los alcance -- señal real, no
#   solo ruido, aunque los indices S1/ST de esa prueba rapida (N_ind=60,
#   base_N=32) seguian siendo demasiado ruidosos para confiar en ellos
#   (S1 negativos, ST hasta 1.90).
#
#   v2 agrega frac_ext como TERCER observable (ademas de t_cese y T_activo,
#   que se mantienen por comparabilidad con Familia 1), con N_ind=50 (mas
#   que los 20 originales, sin llegar a 100 que ya se probo demasiado lento
#   en Python puro para esta ventana temporal de ~146 años).
#
# DECISION DOCUMENTADA -- rango de eta:
#   eta=1e-5 es un parametro NIVEL 2 en el limite del espacio biologicamente
#   plausible (SCD: produccion de daño casi nula). Un rango +-30% simetrico
#   alrededor de un ancla tan cercana a cero puede no representar bien la
#   incertidumbre real de este parametro especifico. Se usa +-20% para eta
#   (en vez de +-30%), y se aplica ademas un piso explicito en 0 (pmax con
#   0) como salvaguarda adicional, aunque con +-20% sobre un ancla positiva
#   el rango ya es positivo por construccion. Los otros 9 parametros
#   mantienen +-30%, igual que en Familia 1, para maxima comparabilidad.
#
# REGLA B.12 -- parsear antes de correr:
#   Rscript -e 'invisible(parse(file="TMCSA_SensibilidadGlobal_Sobol_Hglaber_v2.R"))'
# ==============================================================================

if (!requireNamespace("sensitivity", quietly=TRUE)) install.packages("sensitivity")
library(sensitivity)

N_base <- 1000
N_ind  <- 50   # v2: aumentado de 20 a 50, sin llegar a 100 (demasiado lento)

set.seed(123)
Dc <- 0.60

Xc_NMR <- 17 * (2.5/15)^0.4206
BASE <- list(a=0.0988, b=0.0960, mu=3.0, n=2.0, eta=1e-5, t_lat=12.8,
             Xc=Xc_NMR, beta=54.75, eps=0.180*Xc_NMR/17, CV_eta=0.20)
nombres <- names(BASE)
p <- length(nombres)

GAMMA_EXT_FIJO <- (-log(0.555)/(30.9*365)) * 365
TARGET_HGLABER <- 36.4

# --- Rangos: +-30% para todos, EXCEPTO eta (+-20%, con piso en 0) -----------
rango_bajo  <- sapply(BASE, function(v) v * 0.7)
rango_alto  <- sapply(BASE, function(v) v * 1.3)
rango_bajo["eta"]  <- max(0, BASE$eta * 0.8)   # +-20% para eta, piso en 0
rango_alto["eta"]  <- BASE$eta * 1.2

cat("Rangos usados por parametro:\n")
for (nm in nombres) cat(sprintf("  %-8s [%.6g, %.6g]\n", nm, rango_bajo[nm], rango_alto[nm]))

generar_muestra <- function(N) {
  m <- matrix(NA, nrow=N, ncol=p)
  colnames(m) <- nombres
  for (j in 1:p) m[,j] <- runif(N, rango_bajo[j], rango_alto[j])
  as.data.frame(m)
}

X1 <- generar_muestra(N_base)
X2 <- generar_muestra(N_base)

# ------------------------------------------------------------------------------
# MOTOR -- igual estructura que v1, ahora registrando tambien la fraccion
# de cese por causa externa (frac_ext) en la misma corrida.
# ------------------------------------------------------------------------------
simular_par <- function(params, N_ind=50, seed=1) {
  a <- params[["a"]]; b <- params[["b"]]; mu <- params[["mu"]]; n <- params[["n"]]
  eta <- params[["eta"]]; t_lat <- params[["t_lat"]]; Xc <- params[["Xc"]]
  beta <- params[["beta"]]; eps <- params[["eps"]]; CV_eta <- params[["CV_eta"]]

  kappa_X  <- 0.5
  dt_sem   <- 7/365
  beta_sem <- beta * dt_sem
  eps_sem  <- eps * sqrt(7)
  Gext_sem <- GAMMA_EXT_FIJO * dt_sem
  sig      <- sqrt(log(1 + CV_eta^2))
  muln     <- log(eta) - sig^2/2
  T_ventana <- TARGET_HGLABER * 4
  T_sem    <- round(T_ventana * 365/7)

  set.seed(seed)
  t_ceses   <- numeric(N_ind)
  T_activos <- numeric(N_ind)
  causas_ext <- numeric(N_ind)

  for (i in 1:N_ind) {
    eta_i <- rlnorm(1, muln, sig)
    X <- 0.0; V <- 1.0; T_act <- 0.0; tf <- NA; fue_ext <- 0
    for (sem in 1:T_sem) {
      ta <- (sem - 1) * dt_sem

      if (runif(1) < Gext_sem) { tf <- ta; fue_ext <- 1; break }

      xi <- rnorm(1)
      dX <- eta_i * max(0, ta - t_lat) * 7 -
             beta_sem * X/(kappa_X + X) + eps_sem * xi
      X  <- max(0, X + dX)
      D_i <- X / Xc
      fric <- max(0, (1 - D_i/Dc)^n)
      dV <- a * V * fric - b * V^mu - max(0, D_i - Dc)
      V  <- min(1.10, max(0, V + dV * dt_sem))
      T_act <- T_act + V * dt_sem
      if (X >= Xc || V <= 0.01) { tf <- ta; fue_ext <- 0; break }
    }
    t_ceses[i]    <- if (!is.na(tf)) tf else T_ventana
    T_activos[i]  <- T_act
    causas_ext[i] <- fue_ext
  }
  list(t_cese = median(t_ceses), T_activo = median(T_activos),
       frac_ext = mean(causas_ext))
}

# ------------------------------------------------------------------------------
# DISENO SOBOL -- tres objetos (t_cese, T_activo, frac_ext), mismo diseno
# verificado, corridos una sola vez.
# ------------------------------------------------------------------------------
cat(sprintf("\nFamilia 2 (H. glaber) v2 -- p=%d, N_base=%d, N_ind=%d -> evaluaciones = %d\n",
            p, N_base, N_ind, (p+2)*N_base))

sob_tcese    <- soboljansen(model=NULL, X1=X1, X2=X2, nboot=200)
sob_Tactivo  <- soboljansen(model=NULL, X1=X1, X2=X2, nboot=200)
sob_fracext  <- soboljansen(model=NULL, X1=X1, X2=X2, nboot=200)
stopifnot(isTRUE(all.equal(as.matrix(sob_tcese$X), as.matrix(sob_Tactivo$X))))
stopifnot(isTRUE(all.equal(as.matrix(sob_tcese$X), as.matrix(sob_fracext$X))))
cat("Verificado: los tres disenos son identicos fila por fila.\n")

X_diseno <- sob_tcese$X
n_eval <- nrow(X_diseno)
Y_tcese   <- numeric(n_eval)
Y_Tactivo <- numeric(n_eval)
Y_fracext <- numeric(n_eval)

t_inicio <- Sys.time()
for (i in 1:n_eval) {
  res <- simular_par(X_diseno[i, ], N_ind=N_ind, seed=i)
  Y_tcese[i]   <- res$t_cese
  Y_Tactivo[i] <- res$T_activo
  Y_fracext[i] <- res$frac_ext
  if (i %% 100 == 0 || i == n_eval) {
    transcurrido <- as.numeric(difftime(Sys.time(), t_inicio, units="secs"))
    cat(sprintf("  %d/%d  (%.1fs transcurridos, %.3fs/eval, frac_ext_media_hasta_ahora=%.4f)\n",
                i, n_eval, transcurrido, transcurrido/i, mean(Y_fracext[1:i])))
  }
  if (i %% 500 == 0) {
    saveRDS(list(Y_tcese=Y_tcese[1:i], Y_Tactivo=Y_Tactivo[1:i],
                 Y_fracext=Y_fracext[1:i], ultimo_i=i),
            "TMCSA_Sobol_Hglaber_v2_progreso_parcial.rds")
  }
}
cat(sprintf("\nTiempo total: %.1f min\n", as.numeric(difftime(Sys.time(), t_inicio, units="mins"))))
cat(sprintf("Fraccion externa observada: media=%.4f sd=%.4f min=%.4f max=%.4f\n",
            mean(Y_fracext), sd(Y_fracext), min(Y_fracext), max(Y_fracext)))

sob_tcese   <- tell(sob_tcese,   Y_tcese)
sob_Tactivo <- tell(sob_Tactivo, Y_Tactivo)
sob_fracext <- tell(sob_fracext, Y_fracext)

cat("\n", strrep("=", 78), "\nRESULTADOS -- SOBOL t_cese\n", strrep("=", 78), "\n", sep="")
print(sob_tcese)
cat("\n", strrep("=", 78), "\nRESULTADOS -- SOBOL T_activo\n", strrep("=", 78), "\n", sep="")
print(sob_Tactivo)
cat("\n", strrep("=", 78), "\nRESULTADOS -- SOBOL fraccion de cese externo\n", strrep("=", 78), "\n", sep="")
print(sob_fracext)

hacer_tabla <- function(sob, sufijo) {
  df <- data.frame(
    parametro = nombres,
    S1 = sob$S$original, S1_min = sob$S$`min. c.i.`, S1_max = sob$S$`max. c.i.`,
    ST = sob$T$original, ST_min = sob$T$`min. c.i.`, ST_max = sob$T$`max. c.i.`
  )
  setNames(df, c("parametro", paste0(c("S1","S1_min","S1_max","ST","ST_min","ST_max"), "_", sufijo)))
}
tabla_final <- Reduce(function(x,y) merge(x,y,by="parametro"),
  list(hacer_tabla(sob_tcese,"tcese"), hacer_tabla(sob_Tactivo,"Tactivo"), hacer_tabla(sob_fracext,"fracext")))
write.csv(tabla_final, "TMCSA_SensibilidadGlobal_Hglaber_v2_resultados.csv", row.names=FALSE)
cat("\nGuardado en TMCSA_SensibilidadGlobal_Hglaber_v2_resultados.csv\n")
