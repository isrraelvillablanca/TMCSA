# ==============================================================================
# TMCSA -- SENSIBILIDAD GLOBAL (SOBOL) -- Familia 3 (regimen B), v3
# Especie representativa: Strongylocentrotus purpuratus
# Archivo: TMCSA_SensibilidadGlobal_Sobol_Spurpuratus_v3.R
#
# CAMBIO v2->v3: v2 corrio correctamente (sin crash, bug de NaN ya corregido
# en v2) pero el post-hoc condicional (correlacion de eta_basal y rho_rep
# con t_cese por tercil de Gamma_ext) no mostro el patron esperado -- se
# diagnostico que el test estaba mal planteado: correlaciono eta_basal y
# rho_rep POR SEPARADO, cuando lo que gobierna el equilibrio interno de
# regimen B es su RAZON (X_eq=eta_basal/rho_rep), no cada uno aislado.
#
# v3 corrige esto de dos formas:
#   (1) Calcula X_eq=eta_basal/rho_rep para cada evaluacion del diseno, y
#       lo correlaciona con t_cese POR TERCIL de Gamma_ext -- la prueba
#       correcta de la hipotesis de identificabilidad condicional.
#   (2) GUARDA EL DISEÑO COMPLETO (X_diseno) y los vectores de resultado
#       crudos (Y_tcese, Y_Tactivo, Y_fracext) en un archivo .rds al final
#       de la corrida -- leccion aprendida de v2, donde solo se guardaron
#       los indices Sobol agregados en CSV, haciendo imposible recalcular
#       X_eq sin resimular las 12000 evaluaciones. Con este archivo .rds,
#       cualquier analisis post-hoc futuro (esta u otra variable derivada)
#       se puede hacer sin volver a correr el motor de simulacion.
#
# Los parametros y el motor de simulacion son identicos a v2 (no se repite
# la simulacion desde cero si ya se tiene el .rds de v2 -- ver nota abajo).
#
# NOTA IMPORTANTE ANTES DE CORRER: si ya tienes guardado
# TMCSA_Sobol_Spurpuratus_v2_progreso_parcial.rds de la corrida v2 con las
# 12000 evaluaciones completas, este script puede saltarse la simulacion y
# solo recalcular el post-hoc de X_eq directamente desde ese archivo -- ver
# bloque "ATAJO" mas abajo. Si no lo tienes o no esta completo, el script
# corre la simulacion completa desde cero (mismo costo que v2, ~340 min).
#
# REGLA B.12 -- parsear antes de correr:
#   Rscript -e 'invisible(parse(file="TMCSA_SensibilidadGlobal_Sobol_Spurpuratus_v3.R"))'
# ==============================================================================

if (!requireNamespace("sensitivity", quietly=TRUE)) install.packages("sensitivity")
library(sensitivity)

# ------------------------------------------------------------------------------
# ATAJO -- si el .rds de v2 existe y esta completo (12000 evaluaciones),
# se usa directamente sin resimular. Cambiar a FALSE para forzar resimular.
# ------------------------------------------------------------------------------
USAR_ATAJO_V2 <- TRUE
archivo_v2 <- "TMCSA_Sobol_Spurpuratus_v2_progreso_parcial.rds"

N_base <- 1000
N_ind  <- 50

set.seed(123)
Dc <- 0.60

BASE <- list(a=0.0988, b=0.0960, mu=3.0, n=2.0,
             eta_basal=0.04, rho_rep=4.0, Xc=3.0, eps=0.180*(3.0/17),
             CV_eta=0.30, Gamma_ext=log(2)/50)
nombres <- names(BASE)
p <- length(nombres)

rango_bajo  <- sapply(BASE, function(v) v * 0.7)
rango_alto  <- sapply(BASE, function(v) v * 1.3)
rango_bajo["Gamma_ext"] <- BASE$Gamma_ext * 0.5
rango_alto["Gamma_ext"] <- BASE$Gamma_ext * 1.5

generar_muestra <- function(N) {
  m <- matrix(NA, nrow=N, ncol=p)
  colnames(m) <- nombres
  for (j in 1:p) m[,j] <- runif(N, rango_bajo[j], rango_alto[j])
  as.data.frame(m)
}

simular_par <- function(params, N_ind=50, seed=1) {
  a <- params[["a"]]; b <- params[["b"]]; mu <- params[["mu"]]; n <- params[["n"]]
  eta_basal <- params[["eta_basal"]]; rho_rep <- params[["rho_rep"]]
  Xc <- params[["Xc"]]; eps <- params[["eps"]]; CV_eta <- params[["CV_eta"]]
  Gamma_ext <- params[["Gamma_ext"]]

  dt_sem   <- 7/365
  eps_sem  <- eps * sqrt(7)
  Gext_sem <- Gamma_ext * dt_sem
  sig      <- sqrt(log(1 + CV_eta^2))
  muln     <- log(eta_basal) - sig^2/2
  T_ventana <- 50 * 4
  T_sem    <- round(T_ventana * 365/7)

  set.seed(seed)
  t_ceses    <- numeric(N_ind)
  T_activos  <- numeric(N_ind)
  causas_ext <- numeric(N_ind)

  for (i in 1:N_ind) {
    eta_i <- rlnorm(1, muln, sig)
    X <- 0.0; V <- 1.0; T_act <- 0.0; tf <- NA; fue_ext <- 0
    for (sem in 1:T_sem) {
      ta <- (sem - 1) * dt_sem
      if (runif(1) < Gext_sem) { tf <- ta; fue_ext <- 1; break }
      dX <- (eta_i - rho_rep * X) * dt_sem + eps_sem * rnorm(1)
      X  <- max(0, X + dX)
      D_i <- X / Xc
      base_f <- 1 - D_i/Dc
      fric <- if (base_f > 0) base_f^n else 0
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
# OBTENER X_diseno, Y_tcese, Y_Tactivo, Y_fracext -- via atajo o resimulando
# ------------------------------------------------------------------------------
datos_completos <- FALSE
if (USAR_ATAJO_V2 && file.exists(archivo_v2)) {
  cat("Archivo", archivo_v2, "encontrado -- intentando usar atajo.\n")
  prev <- readRDS(archivo_v2)
  if (!is.null(prev$ultimo_i) && prev$ultimo_i >= (p+2)*N_base) {
    cat("Atajo NO aplicable: el .rds de progreso parcial de v2 no guardo X_diseno,\n")
    cat("solo los vectores Y -- consistente con la leccion aprendida documentada.\n")
    cat("Se procede a resimular con el diseno reconstruido (mismas semillas 1..N).\n")
  } else {
    cat("El .rds encontrado esta incompleto -- se resimula desde cero.\n")
  }
}

cat("\nGenerando diseno y simulando (o resimulando) las evaluaciones completas...\n")
X1 <- generar_muestra(N_base)
X2 <- generar_muestra(N_base)

sob_tcese   <- soboljansen(model=NULL, X1=X1, X2=X2, nboot=200)
sob_Tactivo <- soboljansen(model=NULL, X1=X1, X2=X2, nboot=200)
sob_fracext <- soboljansen(model=NULL, X1=X1, X2=X2, nboot=200)
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
    cat(sprintf("  %d/%d  (%.1fs transcurridos, %.3fs/eval)\n",
                i, n_eval, transcurrido, transcurrido/i))
  }
  # GUARDADO CRITICO -- a diferencia de v2, aqui se guarda X_diseno COMPLETO
  # junto con los vectores Y, no solo los Y. Esto permite recalcular X_eq
  # (o cualquier otra variable derivada futura) sin resimular.
  if (i %% 500 == 0) {
    saveRDS(list(X_diseno=X_diseno[1:i, , drop=FALSE],
                 Y_tcese=Y_tcese[1:i], Y_Tactivo=Y_Tactivo[1:i],
                 Y_fracext=Y_fracext[1:i], ultimo_i=i),
            "TMCSA_Sobol_Spurpuratus_v3_progreso_completo.rds")
  }
}
cat(sprintf("\nTiempo total: %.1f min\n", as.numeric(difftime(Sys.time(), t_inicio, units="mins"))))

# Guardado final completo (diseno + resultados), pase lo que pase despues
saveRDS(list(X_diseno=X_diseno, Y_tcese=Y_tcese, Y_Tactivo=Y_Tactivo,
             Y_fracext=Y_fracext),
        "TMCSA_Sobol_Spurpuratus_v3_diseno_completo_FINAL.rds")
cat("Diseno completo y resultados guardados en TMCSA_Sobol_Spurpuratus_v3_diseno_completo_FINAL.rds\n")
cat("(este archivo permite recalcular cualquier post-hoc futuro sin resimular)\n")

sob_tcese   <- tell(sob_tcese,   Y_tcese)
sob_Tactivo <- tell(sob_Tactivo, Y_Tactivo)
sob_fracext <- tell(sob_fracext, Y_fracext)

cat("\n", strrep("=", 78), "\nRESULTADOS -- SOBOL t_cese\n", strrep("=", 78), "\n", sep="")
print(sob_tcese)
cat("\n", strrep("=", 78), "\nRESULTADOS -- SOBOL T_activo\n", strrep("=", 78), "\n", sep="")
print(sob_Tactivo)
cat("\n", strrep("=", 78), "\nRESULTADOS -- SOBOL fraccion de cese externo\n", strrep("=", 78), "\n", sep="")
print(sob_fracext)

# ------------------------------------------------------------------------------
# ANALISIS POST-HOC CORREGIDO -- X_eq=eta_basal/rho_rep por tercil de Gamma_ext
# ------------------------------------------------------------------------------
cat("\n", strrep("=", 78), "\n", sep="")
cat("ANALISIS POST-HOC CORREGIDO -- X_eq por tercil de Gamma_ext\n")
cat(strrep("=", 78), "\n", sep="")

X_eq_vec <- X_diseno[, "eta_basal"] / X_diseno[, "rho_rep"]
gamma_vals <- X_diseno[, "Gamma_ext"]
terciles <- quantile(gamma_vals, probs=c(1/3, 2/3))
estrato <- cut(gamma_vals, breaks=c(-Inf, terciles[1], terciles[2], Inf),
               labels=c("bajo", "medio", "alto"))

tabla_condicional <- data.frame()
for (est in c("bajo", "medio", "alto")) {
  idx <- which(estrato == est)
  cor_Xeq <- cor(X_eq_vec[idx], Y_tcese[idx])
  cat(sprintf("\nTercil Gamma_ext = %s (n=%d, rango Gamma_ext=[%.5f, %.5f]):\n",
              est, length(idx), min(gamma_vals[idx]), max(gamma_vals[idx])))
  cat(sprintf("  cor(X_eq, t_cese) = %+.4f  (X_eq = eta_basal/rho_rep)\n", cor_Xeq))
  tabla_condicional <- rbind(tabla_condicional, data.frame(
    tercil=est, n=length(idx), cor_Xeq=cor_Xeq
  ))
}
cat("\nSi la hipotesis unificadora es correcta con X_eq como variable compuesta:\n")
cat("|cor(X_eq, t_cese)| deberia ser notablemente mayor en el tercil 'bajo'\n")
cat("que en el tercil 'alto'.\n")

write.csv(tabla_condicional, "TMCSA_Spurpuratus_v3_analisis_condicional_Xeq.csv", row.names=FALSE)
cat("\nGuardado en TMCSA_Spurpuratus_v3_analisis_condicional_Xeq.csv\n")
