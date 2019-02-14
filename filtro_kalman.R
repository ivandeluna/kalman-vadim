ipc <- read.csv("ipc.csv")

X <- ipc.ret

# Función para estimar la raíz del error mediano cuadrático
RMSE <- function(e,r){    # e - valores estimados, r - valores reales
  sqrt(mean((e-r)^2))
}

mu <- mean(X)      # Media de la st
sig <- var(X)    # Sigma de la st
sd <- sqrt(sig)    # D.E. de la st
size <- length(X)  # Tama;o de la st
x_1 <- X[1]        # primer valor de X
z <- rnorm(size, mu, sd) # distribucin Z de X


# Preparacion datos a estimar
xhat <- vector(mode = "numeric", length = size)               # Vector de estimación de x a posteriori
P <- vector(mode = "numeric", length = size)                  # Estimación de e a posteriori
xhatminus <- vector(mode = "numeric", length = size)          # Estimación de x a priori
Pminus  <- vector(mode = "numeric", length = size)            # Estimación de e a priori
K <- vector(mode = "numeric", length = size)                  # Ganancia Kalman o blending factor
xreal <- vector(mode = "numeric", length = size)              # un valor real de x  
z_est <- vector(mode = "numeric", length = size)              # una observacion de x

R <- 0.1^2         # estimación de la varianza

# Estimaciones iniciales
xhat[1] = 0.0     # valor inicial del vector xhat
P[1] <- 0.03434   # valor inicial del vector P
xreal[1] <- x_1[1]     # valor inicial del vector xreal
a <- 1.1           # valor inicial de a

# Filtro Kalman Clásico

for (t in 2:size) {
  # actualización del tiempo
  xhatminus[t] <- xhat[t-1]
  Pminus[t] <- P[t-1]+ sig
  
  # actualización de la medición
  K[t] <- Pminus[t]/(Pminus[t]+R)
  xhat[t] <- xhatminus[t] + K[t] * (z[t] - xhatminus[t])
  P[t] <- (1 - K[t]) * Pminus[t]
}

plot(z)
lines(xhat)

# Filtro Kalman Vadim
Pmax <- max(P)
xhat[1] <- 0.0
P[1] <- 0.3434
xreal[1] <- x_1[1]
a = 1.1

for (t in 2:size) {
  # Sistema real
  xreal[t] <- a * xreal[t-1]
  z_est[t] <- xreal[t] + z[t]
  
  # Actualización del tiempo
  xhatminus[t] <- a * xhat[t-1]
  Pminus[t] <- (a^2) * P[t-1] + sig
  
  # Actualización de la medición
  K[t] <- Pminus[t]/(Pminus[t] + R)
  xhat[t] <- xhatminus[t] + K[t] * (z[t] - xhatminus[t])
  P[t] <- Pmax
} 

plot(z)
lines(xhat)



RMSE(xhat, z)
