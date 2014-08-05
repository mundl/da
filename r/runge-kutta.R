# differential equation y'=sin(x)^2 * y
dy <- function(x, y, parms, ...) {
  return(sin(x)^2 * y)
}

# exact solution 
exact <- function(x) 2 * exp(0.5 * (x - sin(x) * cos(x)))

# euler's method
euler <- function(x, y, h, fun) {
  y1 <- y + h * fun(x, y)
  return(c(x + h, y1))
}

# heun's method
heun <- function(x, y, h, fun) {
  yp <- y + h * fun(x, y)
  y1 <- y + 0.5 * h * (fun(x, y) + fun(x+h, yp))
  return(c(x + h, y1))
}

# classical Runge–Kutta method
runge <- function(x, y, h, fun) {
  k <- numeric(4)
  k[1] <- h * fun(x, y)
  k[2] <- h * fun(x + h / 2, y + k[1] / 2)
  k[3] <- h * fun(x + h / 2, y + k[2] / 2)
  k[4] <- h * fun(x + h, y + k[3]) 
  
  return(c(x + h, y + 1 / 6 * sum(k * c(1, 2, 2, 1))))
}

# step size = 0.5, last value = 5
h <- 0.5
niter <- 5/h
run <- eul2 <- eul <- heu <- data.frame(x=0, y=exact(0))

for(i in seq_len(niter)+1) {
  eul[i, ] <- euler(x=eul$x[i-1], y=eul$y[i-1], h=h, fun=dy)
  heu[i, ] <- heun (x=heu$x[i-1], y=heu$y[i-1], h=h, fun=dy)
  run[i, ] <- runge(x=run$x[i-1], y=run$y[i-1], h=h, fun=dy)
}

# euler's method with reduced step size
h <- 0.5/2
niter <- 5/h
for(i in seq_len(niter)+1) {
  eul2[i, ] <- euler(x=eul2$x[i-1], y=eul2$y[i-1], h=h, fun=dy)
}

# evaluating exact solution at 
x <- seq(0, 5, 0.1)

# concatenating the methods into a data.frame
odesolve <- rbind(data.frame(x=x, y=exact(x), method="Exact Solution"),
                  data.frame(run,  method="Runge-Kutta method"),
                  # data.frame(heu,  method="Heun's method"),
                  data.frame(eul2, method="Euler's method (reduced step size)"),
                  data.frame(eul,  method="Euler's method"))

# translating into german
odesolve$method <- factor(odesolve$method, 
                          levels=c("Exact Solution", "Runge-Kutta method", 
                                   "Heun's method", 
                                   "Euler's method (reduced step size)", 
                                   "Euler's method"),
                          labels=c("Exakte Lösung", "Runge-Kutta", 
                                   "Heun", "Euler (halbe Schrittweite)", 
                                   "Euler"))

p <- ggplot(odesolve, aes(x=x, y=y, col=method)) +   geom_line() + 
  geom_point(data=subset(odesolve, as.numeric(method)!=1)) +
  scale_color_discrete("") + 
  mytheme + 
  theme(plot.title=element_blank())

ggsave("figure/runge-kutta.pdf", width=width, height=width/2, plot=p, unit="mm")
