# Fitting and Predicting VaR basd on an ARMA-GARCH process
# Marius Hofert

library(rugarch)
library(qrmtools)

###############################################################################
## Especificación del modelo
###############################################################################
nu <- 3 # grados de libertad de la distribución estandarizada de Z_t
fixed.p <- list(mu = 0,       # media (intercepto)
                ar1 = 0.5,    # phi_1 (AR(1) parametro de mu_t)
                ma1 = 0.3,    # theta_1 (MA(1) parametro de mu_t)
                omega = 4,    # alpha_0 (intercepto)
                alpha1 = 0.4, # alpha_1 (GARCH(1) parametro de sigma_t^2))
                beta1 = 0.2,  # beta_1 (GARCH(1) parametro de sigma_t^2)
                shape = nu)   # grados de livertad para las inovaciones estandarizadas de t_nu
armaOrder<- c(1,1)            # especificación del modelo ARMA
garchOrder <- c(1,1)          # especificación del modelo GARCH
varModel <- list(model = "sGARCH", garchOrder = garchOrder)
spec <- ugarchspec(varModel, mean.model = list(armaOrder = armaOrder),
                   fixed.pars = fixed.p, distribution.model = "std") # residuales t estandarizados

###############################################################################
# Simulación
###############################################################################
n <- 1000                    # tamaño de la muestra
x <- ugarchpath(spec, n.sim = n, m.sim = 1, rseed = 271) # n.sim largo de la trayectoria  simulada, m.sim numero de caminos

# Extraer las series resultantes
X <- fitted(x) # Proceso simulado de X_t = mu_t + epsilon_t para epsilon_t = sigma_t * Z_t
sig <- sigma(x) # volatilidad de sigma_t (desviaciones estandar condicionales)
eps <- x@path$residSim # residuos no estandarizados de epsilon_t = sigma_t * Z_t

# Checar que existan los valores
stopifnot(all.equal(X, x@path$seriesSim, check.attributes = FALSE),
          all.equal(sig, x@path$sigmaSim, check.attributes = FALSE))

plot(X, type = "l", xlab = "t", ylab = expression(X[t]))

plot(sig, type = "h", xlab = "t", ylab = expression(sigma[t]))

plot(eps, type = "h", xlab = "t", ylab = expression(epsilon[t]))

###############################################################################
# Estimar un modelo ARMA-GARCH a la serie simulada
###############################################################################
spec <- ugarchspec(varModel, mean.model = list(armaOrder = armaOrder),
                   distribution.model = "std") # sin parametros fijos

# Para utilizar los datos del IPC, quitar comentario a la siguiente linea:
#X <- ipc.ret
fit <- ugarchfit(spec, data = X)

# Extraer las series resultantes
mu. <- fitted(fit) # estimación de hat(mu)_t ( = hat(X)_t)
sig. <- sigma(fit) # estimación de hat(sigma)_t

# Checar que los valores sean válidos
stopifnot(all.equal(as.numeric(mu.), fit@fit$fitted.values),
          all.equal(as.numeric(sig.), fit@fit$sigma))

plot(X, type = "l", xlab = "t",
     ylab = expression("Data"~X[t]~"and fitted values"~hat(mu)[t]))
lines(as.numeric(mu.), col = adjustcolor("blue", alpha.f = 0.5))
legend("bottomright", bty = "n", lty = c(1,1),
       col = c("black", adjustcolor("blue", alpha.f = 0.5)),
       legend = c(expression(X[t]), expression(hat(mu)[t])))

resi <- as.numeric(residuals(fit))
stopifnot(all.equal(fit@fit$residuals, resi))
plot(resi, type = "l", xlab = "t", ylab = expression(epsilon[t]))

# Gráfica Q-Q de los residuales estandarizados de Z_t contra la t especificada
# (t_nu con varianza de 1)
Z <- fit@fit$z
stopifnot(all.equal(Z, as.numeric(resi/sig.)))
qq_plot(Z, FUN = function(p) sqrt((nu-2)/nu) * qt(p, df = nu),
  main = substitute("Q-Q plot of ("*Z[t]*") against a standardized"~t[nu],
                                      list(nu. = nu)))

###############################################################################
# Calcular la serie de tiempo del VaR
###############################################################################
alpha <- 0.99

# Extraer los valores estimados de VaR_alpha
VaR. <- as.numeric(quantile(fit, probs = alpha))

# Construir manualmente y comparar con dos
nu. <- fit@fit$coef["shape"] # extrae los grados de libertad estimados de nu
VaR.. <- as.numeric(mu. + sig. * sqrt((nu.-2)/nu.) * qt(alpha, df = nu.)) # VaR manual
stopifnot(all.equal(VaR.., VaR.))

btest <- VaRTest(1-alpha, actual = -X,
                 VaR = quantile(ugarchfit(spec, data = -X), probs = 1-alpha))
btest$expected.exceed # Number of expected excedances = (1-alpha) * n
btest$actual.exceed # Number of actual exceedances
# Prueba incondicional
btest$uc.H0 # hipótesis nula
btest$uc.Decision # decisión de la prueba
# Prueba condicional
btest$cc.H0 # hipótesis nula
btest$cc.Decision # decisión de la prueba

###############################################################################
# Precedir VaR basado en el modelo estimado
###############################################################################

fspec <- getspec(fit) # especificaciones del proceso estimado
setfixed(fspec) <- as.list(coef(fit)) # usar los parametros de la estimación
m <- ceiling(n/10) # numero de pasos para el pronóstico (roll/iterate m-1 times forward with frequency 1)
pred <- ugarchforecast(fspec, data = X, n.ahead = 1, n.roll = m-1, out.sample = m) # predecir del proceso estimado

# Extraer los resultados de las series
mu.predict <- fitted(pred) # extraer los X_t predecidos
sig.predict <- sigma(pred) # extraer la sigma_t predecida
VaR.predict <- as.numeric(quantile(pred, probs = alpha)) # valores VaR_alpha correspondientes con las predicciones

# Checar valores
stopifnot(all.equal(mu.predict, pred@forecast$seriesFor, check.attributes = FALSE),
          all.equal(sig.predict, pred@forecast$sigmaFor, check.attributes = FALSE))

VaR.predict. <- as.numeric(mu.predict + sig.predict * sqrt((nu.-2)/nu.)*qt(alpha, df = nu.))
stopifnot(all.equal(VaR.predict., VaR.predict))

###############################################################################
# Simular trayectorias futuras de X_t y calcular el VaR correspondiente
###############################################################################
B <- 1000
set.seed(271)
X.sim.obj <- ugarchpath(fspec, n.sim = m, m.sim = B) #simlar trayectorias futuras

# Calcular el VaR_alpha simulado y sus intervalos de confianza correspondinetes
X.sim <- fitted(X.sim.obj)
sig.sim <- sigma(X.sim.obj)
eps.sim <- X.sim.obj@path$residSim
VaR.sim <- (X.sim - eps.sim) + sig.sim * sqrt((nu.-2)/nu.) * qt(alpha, df= nu.)
VaR.CI <- apply(VaR.sim, 1, function(x) quantile(x, probs = c(0.05, 0.975)))

###############################################################################
# Graficar
###############################################################################

# Set up
yran <- range(X, # trayectoria simulada
              mu., VaR., # valores estimados de la media condicional y VaR_alpha
              mu.predict, VaR.predict, VaR.CI) # valores predecidos
myran <- max(abs(yran))
yran <- c(-myran, myran) # rango Y
xran <- c(1, length(X) + m) # rango X

# Datos simulados de X, valores condicionales estimados de la media, mu_t, y VaR_alpha
plot(X, type = "l", xlim = xran, ylim = yran, xlab = "Time t", ylab = "",
     main = "Simulated ARMA-GARCH, fit, VaR, VaR predictions and CI")
lines(as.numeric(mu.), col = adjustcolor("darkblue", alpha.f = 0.5))
mtext(paste0("Expected exceed.:", btest$expected.exceed, "    ",
             "Actual exceed.:", btest$actual.exceed, "    ",
             "Test:", btest$cc.Decision),
      side = 4, adj = 0, line = 0.5, cex = 0.9)

# Predictions
t. <- length(X) + seq_len(m) # puntos en el tiempo futuro
lines(t., mu.predict, col = "blue")
lines(t., VaR.predict, col = "red")
lines(t., VaR.CI[1,], col = "orange")
lines(t., VaR.CI[2,], col = "orange")
legend("bottomright", bty = "n", lty = rep(1,6), lwd = 1.6,
       col = c("black", adjustcolor("darkblue", alpha.f = 0.5), "blue",
                                    "darkred", "red", "orange"),
       legend = c(expression(X[t]), expression(hat(mu)[t]),
                  expression("Predicted"~mu[t]~"or("~X[t]*")"),
                  substitute(widehat(VaR)[a], list(a = alpha)),
                  substitute("Predicted"~VaR[a], list(a = alpha)),
                  substitute("95%-CI for"~VaR[a], list(a = alpha))))

