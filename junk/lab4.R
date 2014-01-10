#+ setup, include = FALSE
opts_chunk$set(fig.path = 'figure/lab4-')
options(show.signif.stars = F)
par.old <- par(no.readonly = TRUE)

#' # Lab 4: Kriging 
library(geoR)
library(ggplot2)
library(reshape2)

#' ### Generating a gaussian random field
#' `Read about the command grf() wich can be used to generate Gaussian random fields. Choose some variogram function and generate a realisation of a grf with` $n$`=105. Keep five points in a seperate dataset, you will use them for model comparison.`

#my.data <- grf(n=105, cov.pars=c(1, 0.3))
load("grf2.RData")

mask <- rep(c(T, F), times=c(100, 5))
x.train <- subset(my.data, mask)
x.test  <- subset(my.data, !mask)

#' ### Variogram estimation
#' `Now try to estimate the variogram. Plot the observed data. Plot the empirical variogram. Try different functional forms an starting parameters. Try both, variofit() and likfit() for each functional form. Describe your results.`
 
#+ data, fig.width=8, fig.height=6
ggplot(as.data.frame(my.data), aes(x, y, col=data)) + geom_point(size=5)

#+ emp-vario, fig.width=8, fig.height=4, echo=2
par(mfrow=c(1, 2), mar=c(2, 4, 2, 2) + 0.1)
plot(my.data, plot.location=T, ask=F, main="empirical variogram")
par(par.old)

#' ### Testing various functional forms
vario.emp <- variog(my.data)

forms <- c("matern", "exponential", "gaussian", "spherical", "circular", 
           "cubic", "powered.exponential", "cauchy", 
           "gneiting", "pure.nugget")

#' The blue line the estimation done by `variofit()` whereas the black one is the one done by `likfit()`. 
#+ est-vario, fig.width=10, fig.height=7, echo=-c(1, 4), tidy=F
par(mfrow=c(3, 4), mar=c(2, 2, 2, 2) + 0.1)
krige.obj <- list()
for (i in forms) {
  plot(vario.emp, main=i, xlab="", ylab="")
  krige.obj[["variofit"]][[i]] <- 
    variofit(variog(x.train, messages=F), cov.model=i, messages=F)
  lines(krige.obj[["variofit"]][[i]], col=4)
  
  krige.obj[["likfit"]][[i]] <- 
    likfit(x.train, cov.model=i, 
           ini.cov.pars=krige.obj[["variofit"]][[i]]$cov.pars+0.01, 
           messages=F)
  lines(krige.obj[["likfit"]][[i]])
}
par(par.old)

#' ### Testing different initial values
extract_estimates <- function(x) {
  round(c(nugget=x$nugget, range=x$cov.pars[1], sill=x$cov.pars[2]), 5)
}

inits <- expand.grid(seq(0, 2, length.out=5), seq(0.1, 2, length.out=5))
estimates <- list()
for (i in seq_len(nrow(inits))) {
  estimates[[i]] <- 
    extract_estimates(suppressWarnings(variofit(variog(x.train, messages=F), 
                                              cov.model="matern", 
                                              ini.cov.pars=inits[i, ],
                                              messages=F)))
}
estimates <- do.call(rbind, estimates)
estimates

#' Only fifth digit after the comma changes, so the starting values are practically without effect. 

#' ### Do the actual kriging
#' `Choose about five different variogram functions and perform conventional kriging (krige.conv()) with them. Plot the resulting surfaces.`
forms <- c("matern", "exponential", "spherical", "powered.exponential", 
           "cauchy", "pure.nugget")
grid <- expand.grid(seq(0, 1, length.out=100), seq(0, 1, length.out=100))

#+ krige, fig.width=10, fig.height=7
par(mfrow=c(2, 3), mar=c(2, 2, 2, 2) + 0.1)
for (i in forms) {
  kc <- suppressWarnings(
    krige.conv(x.train, locations=grid, 
               output=list(messages=F),
               krige=krige.control(obj.model=krige.obj[["likfit"]][[i]])))
  image(kc, main=i, xlab="", ylab="")
  points(my.data, add=T)
}
par(par.old)

#' ### Estimates the five points and their residuals
#' `For each model, estimate the predicted and observed values for the five points put aside. What can you say about these models?`
resid <- list()
for (i in forms) {
  kc <- suppressWarnings(
    krige.conv(x.train, locations=x.test$coords, 
               output=list(messages=F),
               krige=krige.control(obj.model=krige.obj[["likfit"]][[i]])))
  resid[[i]] <- x.test$data - kc$predict
}
resid <- do.call(rbind, resid)

#+ residuals, fig.width=8, fig.height=6
ggplot(melt(resid), aes(Var2, value, fill=Var1)) + geom_bar(stat="identity", position="dodge") + 
  scale_fill_discrete("model") + labs(x="point number", y="residual")

#' Except for the pure nugget variogramm, all the models have similiar residuals and the generated surfaces are more or less the same. 

#' ### Cross Validation
#' `Perform one-out cross validation (xvalid()). What model would you choose as your best model? How close is it to the true one?`
xv <- list()
for (i in forms) {
  xv[[i]] <- xvalid(geodata=x.train, model=krige.obj[["likfit"]][[i]], 
                    messages=F)
}

me <- function(error, ...) mean(error)
mse <- function(error, ...) mean(error^2)
msdr <- function(error, krige.var, ...) mean((error)^2/krige.var)
ssr <- function(error, ...) sum((error)^2)
funs <- list(me=me, mse=mse, msdr=msdr, ssr=ssr)

diagnostics <- sapply(xv, function(x) sapply(funs, function(y) do.call(y, x)))
round(diagnostics, 3)

#' For the plot, we are going to leave out the "pure nugget" model, to have sensible limits.
tmp <- melt(diagnostics[, colnames(diagnostics)!="pure.nugget"])
#+ diagnostics, fig.width=8, fig.height=5
ggplot(tmp, aes(value, Var2)) + 
  geom_point() + 
  facet_wrap(~Var1, scales = "free_x") + 
  labs(x = "", y = "")

#' We'd choose the cauchy model as the best one, because it has the lowest MSE.