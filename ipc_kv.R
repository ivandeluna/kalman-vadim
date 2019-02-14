# Aplicación del Filtro Kalman-Vadim
# Para la estimación de la volatilidad en
# Series de Tiempo del IPC

# Iván de Luna


# Lista de paquetes requeridos y checa si se encuentran instalados
# si no se encuentran, se instalan.

lop <- c("MASS", "stats", "tseries", "DAAG",
         "xts","forecast", "urca")
np <- lop[!(lop %in% installed.packages()[,"Package"])]
if (length(np)) install.packages(np)

# Se cargan los paquetes requeridos
lapply(lop, require, character.only = TRUE)

# Se carga el archivo ipc.csv y amxl.csv con los datos del IPC
# y de América Movil de los últimos 5 años

ipc <- read.csv("ipc.csv")
amx <- read.csv("amxl.csv")

# Checar los datos que contiene la variable ipc
head(ipc)
head(amx)

# Darle formato a las fechas tipo POSIX

ipc.dates <- as.POSIXct(ipc$Date, format = "%Y-%m-%d")
ipc.dates <- as.Date(ipc$Date)

amx.dates <- as.POSIXct(amx$Date, format = "%Y-%m-%d")
amx.dates <- as.Date(amx$Date)

# Checar que las fechas estén en el formato deseado

head(ipc.dates)
head(amx.dates)

# creamos un vector del cierre del ipc y amx

ipc.close <- as.numeric(ipc$Close)
amx.close <- as.numeric(amx$Close)

# Graficamos las series
plot(ipc.close)
plot(amx.close)

plot.ts(ipc.close)
plot.ts(amx.close)

# Pero las series no tienen tiempo, por eso hay que agregarles el
# vector del tiempo

ipc.s <- data.frame(ipc.dates, ipc.close)
amx.s <- data.frame(amx.dates, amx.close)

plot(ipc.s)
plot(amx.s)


plot.ts(ipc.s$ipc.close)
plot.ts(amx.s$amx.close)

# Los rendimientos se calculan como la diferencia logarítmica
# entre los precios del punto t y t-1, esto se puede hacer
# declarando una función o utilizando alguna ya disponible

rend <- function(x) diff(log(x))

ipc.ret <- rend(ipc.close)
amx.ret <- rend(amx.close)

plot.ts(amx.ret)

# Para conocer la estacionalidad de la serie de tiempo se pueden utilizar
# diversas funciones como la de autocorrelación

acf(amx.ret)
acf(ipc.ret)

pacf(amx.ret)
pacf(ipc.ret)

# también puede utilizarse el test de Dick-Fuller
adf.test(ipc.ret, alternative = "stationary", k=0)
adf.test(amx.ret, alternative = "stationary", k=0)

# Podemos entender que la serie es estacional con la primera diferencia
# vamos a modelar un ARIMA del tipo (0,0,0), (1,0,0), (0,0,1) y (1,0,1)
# para ver la diferencia y cuál queda mejor

# ARIMA (0,0,0)
a1 <- arima(ipc.ret, c(0,0,0), seasonal = list(order = c(0,0,0),
                                               period = 12))

b1 <- Arima()
a1
p1 <- predict(a1, n.ahead = 12)
p1

# ARIMA (1,0,0)
a2 <- arima(ipc.ret, c(1,0,0), seasonal = list(order = c(0,0,0),
                                               period = 12))
a2

p2 <- predict(a2, n.ahead = 12)
p2

# ARIMA (0,0,1)
a3 <- arima(ipc.ret, c(0,0,1), seasonal = list(order = c(0,0,1),
                                               period = 12))
a3

p3 <- predict(a3, n.ahead = 12)
p3

# ARIMA (1,0,1)
a4 <- arima(ipc.ret, c(1,0,1), seasonal = list(order = c(1,0,1),
                                               period = 12))
a4

p4 <- predict(a4, n.ahead = 12)
p4

# ARIMA (0,1,1)
a5 <- arima(ipc.ret, c(0,1,1), seasonal = list(order = c(0,1,1),
                                               period = 12))
a5

p5 <- predict(a5, n.ahead = 12)
p5

# ARIMA (1,1,1)
