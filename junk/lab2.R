#+ setup, include = FALSE
opts_chunk$set(fig.path = 'figure/lab2-')
options(show.signif.stars = F)
par.old <- par(no.readonly = TRUE)

#' # Lab 2: Analysis of Point Processes in R
#' First, load the packages,
library(ggplot2)
library(reshape2)
library(spatstat)

#' do some global assignments,
n <- 1000    # number of iterations used for Monte Carlo and bootstrapping


#' and define general useful functions.
format_p <- function(p, n) ifelse(p < (1/n),  paste("<", 1/n, sep=""), p)
p_2sided <- function(sim, obs) mean(abs(sim-mean(sim)) >= abs(obs - mean(sim)))


#' ## Question 1: Examining the Bramble Cane Dataset.
#' ### 1a) Exploratory analysis
#' `Read the brief description of the bramble can data. Plot the unmarked data and describe what you see. Do you think there are patterns, clusters, features of interest?`

#' To avoid typing, shorten the name ;-)
brambles <- bramblecanes

#+ plot-brambles, fig.width=8, fig.height=8, tidy=F
ggplot(as.data.frame(brambles), aes(x, y, shape = marks, col = marks)) + 
  geom_point() + coord_fixed() + 
  scale_color_discrete("age") + scale_shape_discrete("age") +
  labs(x = "", y = "") 

#' There seems to be a pattern associated with the old plants (horizontal roots, stolon), an envelope of young plants around the old ones. There are gaps in the middle of the sample plot, consistent with a peripatetic growth pattern.

summary(brambles)

#' Only 10% of the plants are two years old.


#' ### 1b) Quadrat test
#' `Perform the quadrat.test and interpret the output. Try different resolutions.`
q.test <- list()
for (i in 2:12) {
  q.test[[as.character(i)]] <- quadrat.test(brambles , nx=i)
}
q.test[[1]]
p.val <- data.frame(nx=as.numeric(names(q.test)), 
                    p=sapply(q.test, function(x) x$p.value))

#+ 	quadrat-res, fig.width=8, fig.height=3, tidy=F				  
ggplot(p.val, aes(as.factor(nx), p, group=1)) +  
  geom_line() + geom_point() + 
  labs(x="quadrat counts in x resp. y", y="p value")

#' The `quadrat.test()` function  is testing the null hypothesis that the data pattern is spatially random. The null hypothesis is rejected because the p-value is well below a significance level of $\alpha=0.05$. Even for the smallest meaningful resolution of $nx=ny=2$. 

#' ### 1c) G-function
#' `Write a program to evaluate the G-function for your data. Use Gest()-function for comparison.`  
#' 
#' G-function basically is the empirical cumulative density function of the random distance from the typical point to its nearest neighbour. First we have to write a function `dnearest()` which for every point returns the distance to its nearest neighbour. Objects of class ppp get coerced. 

dnearest <- function(x, ...) {
  if (class(x)=="ppp") x <- as.data.frame(x)[, c("x", "y")]
  d <- as.matrix(dist(x, ...))
  diag(d) <- Inf 
  
  r <- apply(d, 1, min)
  return(r)
}

#' `my_G()` returns the empirical densities for the distances to the nearest neighbours. 

my_G <- function(x, ...) {
  r <- dnearest(x, ...)
  y <- data.frame(r, G=ecdf(r)(r))
  y <- y[order(y$r), ]
  rownames(y) <- NULL
  return(y)
}

#' `Use plot(envelope(brambles, Gest)) to plot the estimated G-function. How would you interpret the plot? Is there evidence against complete spatial randomness?`  
 
#' In the following figure the G-function `Gest()` provided by package `spatstat` is plotted together with `my_G()` (red dots) and envelopes of the theoretical G-function (grey area). 

#+ g-fun, fig.width=8, fig.height=5
#invisible(plot(spatstat::envelope(Y=brambles, fun=Gest, verbose=F), main=""))
invisible(plot(envelope(Y=brambles, fun=Gest, verbose=F), main=""))
points(G ~ r, data=my_G(x=brambles), col=2, lwd=2, cex=0.5) 

#' As the observed G-function is everywhere outside of the simulated envelope, it is reasonable to assume that our process **is not** spatially independent, based on that evidence. 


#' ### 1d) Monte Carlo test of complete spatial randomness
#' `Perform a MC-test of complete spatial randomness of the bramble cane data, using the average distance to the nearest neighbour as your statistic. Carefully describe each step of your test. What are your results and conclusions?`  

#' * generate $n$ samples with $x$ and $y$ following a uniform distribution,
#' * for each sample: compute the average distance to the nearest neighbour
#' * compute p-value for the null hypothesis of CSR (mean of the number of samples where $d_{obs} >= d_{sim}$)

d.obs <- mean(dnearest(brambles))
d.sim <- replicate(n, mean(dnearest(data.frame(x=runif(brambles$n), 
                                               y=runif(brambles$n)))))

p.val <- mean(d.sim <= d.obs)
format_p(p.val, n)
#+ dmean-dens, fig.width=8, fig.height=3, tidy=F
qplot(d.sim, geom = "density", xlim=c(0, max(d.sim))) + 
  geom_vline(xintercept = d.obs, col = 2, linetype = 2)

#' The average distance to nearest neighbour is a maximum for a Poisson process. As there is not a single sample in the Monte Carlo simulation where the average distance to nearest neighbour is smaller than in the observed data, the p-value is smaller than $1/n$. Therefore the sample is not spatially random.  
# '  

#' `Repeat the previous steps by applying a maximum statistic test to the G-function. What are your results and conclusions?`

dist2G <- function(x, agg=max) {
  Gsim <- Gest(x)
  match.fun(agg)(abs(Gsim$rs - Gsim$theo))
  match.fun(agg)(Gsim$rs - Gsim$theo)
}

d.obs <- dist2G(brambles, agg=max)
d.sim <- replicate(n, dist2G(ppp(runif(brambles$n), runif(brambles$n)), agg=max))

p.val <- mean(d.sim >= d.obs)
format_p(p.val, n)

#+ dmax-dens, fig.width=8, fig.height=3, tidy=F
qplot(d.sim, geom = "density") + 
  geom_vline(xintercept = d.obs, col = 2, linetype = 2)

#' For a maximum statistic the null hypothesis states that there is complete spatial randomness. A small maximum distance of the observed G-function to the theoretical G-function of a Poisson process, indicates spatial randomness. Because the p-value is less than $1/n$, this is not the case. 

#' ## Question 2: Examining a Poisson Process
#' `Simulate a realisation of SCR process with `$n=100$` using the ` ppp() ` function. Perform the density analysis, estimate G-, J-, K-, and L-functions and the maximum statistic of the G-function for your realisation.`

set.seed(pi)
csr.sim <- ppp(runif(100), runif(100))

#+ fig.width=8, fig.height=6
plot(density(csr.sim))

#+ gjkl-funs, fig.width=8, fig.height=6, echo=2, tidy=F
par(mfrow = c(2, 2), mar = c(2, 4, 0, 2) + 0.1)
invisible(lapply(list(Gest, Jest, Kest, Lest), 
                 function(x) plot(x(csr.sim), main = "", 
                                  legendargs = list(cex = 0.7))))
par(par.old)

#+ maxstat
g.sim <- Gest(csr.sim)
max(abs(g.sim$theo - g.sim$km))

#' `What would you conclude with regard to the hypothesis of complete spatial randomness? Why do your estimated functions differ from the theoretical one?` 
#'    

#' The theoretical {G,J,K,L}-functions are similar to the observed ones, any variation is due to randomness. The maximum statistic also shows a rather small value. Therefore it is reasonable to accept the null hypothesis. 


#' ## Question 3: Examining a marked point process
#' ### 3a) relationship of different ages
#' `Use the Gcross() function to investigate relationship between bramble cane plants of different ages. What conclusions can you draw?`
pairs <- combn(levels(brambles$marks), 2)

#+ gcross, echo=2
par(mfrow=c(3, 1), mar=c(2, 4, 0, 2) + 0.1)
invisible(apply(pairs, 2, function(x) plot(Gcross(X=brambles, i=x[1], j=x[2]), 
                                           main="", xlim=c(0, 0.09))))
par(par.old)



#' ### 3b) Monte Carlo Test
#' `Construct and implement an MC-test to check whether younger plants tend to cluster more than bramble cane plants on average.`  

#' If there is clustering than at least for some points in our sample the distances to their nearest neighbours have to be smaller, compared to unclustered observations. This means, if a G-function is shifted to the left (smaller distances occur more frequently), there is clustering.  
#' 

#' Plants of the age 0 are considered to be young and are compared to the average (all plants). A maximum statistic is used to test the null hypothesis that younger plants are not more clustered than all the plants.  
#' 
#' $H_0: G_{young} - G_{all} \le 0$  
#' $H_1: G_{young} - G_{all}  > 0$
#' 
#' For the Monte Carlo sampling, we radomly choose 359 points and assume that they are of age 0. To assure, that the G-functions always get evaluated at the same distances $r$, we setup this vector.
r <- seq(0, 0.102, length.out=513)
g.obs <- data.frame(r, 
                    all=Gest(brambles, r=r)$km, 
                    old=Gest(brambles[brambles$marks!=0, ], r=r)$km,
                    young=Gest(brambles[brambles$marks==0, ], r=r)$km)
max.obs <- max(g.obs$young - g.obs$all)

sample_young <- function(x) x[sample(x=x$n, size=sum(x$marks==0)), ]
g.young.sim  <- replicate(n, Gest(sample_young(brambles), r=r)$km)

#+ mc-gfun, fig.height = 5, fig.width = 8, tidy = F
ggplot(melt(g.obs, id.vars = "r"), aes(r, value, col = variable)) +
  geom_line(data = melt(data.frame(r, g.young.sim), id.vars = "r"), 
            mapping = aes(r, value, group = variable), col = "grey")+ 
  geom_line() +
  labs(x = "distance", y = "cumulative probability", title = "G-function") + 
  scale_color_discrete("")

max.sim <- apply(g.young.sim, 2, function(x) max(x - g.obs$all))
#+ mc-dens, fig.height = 3, fig.width = 8, tidy = F
qplot(max.sim, geom = "density") + 
  geom_vline(xintercept = max.obs, col = 2, linetype = 2)

p.val <- mean(max.sim >= max.obs)
format_p(p.val, n)

#' As we can also see in the plot of the G-functions, all plants are always more clustered than any subset (young plants, old plants). Because the p-value is larger than $\alpha=0.05$ we accept the null hypothesis. 


#' ### 3c) Bootstrap test of independent marking
#' `Use bootstrap to test the hypothesis of independent marking. Carefully describe your assumptions and procedure. What does the hypothesis mean in non-statistical terms? Interpret the results of the test.`  
#' 

#' For the null hypothesis we assume that marks are independent, that is G-functions are the same.  
#' 
#' $H_0: G_{i,j} = G_{j,i} = G_{j,k} = G_{i,k}$ ... 

#' 
#' to obtain the bootstrap, do $n$ times:   
#' * For all nine combinations of levels the $Gcross()$ function will be computed.  
#' * Calculate the spread/range of the envelope for every $r$, just keep the value. 
#'    

max_diff_env <- function(x){
  pairs <- expand.grid(levels(x$marks), levels(x$marks))
  tmp <- apply(pairs, 1, function(y) Gcross(X=x, i=y[1], j=y[2], r=seq(0, 0.22, length.out=1000)))
  g.sim <- sapply(tmp, function(x) x$km)
  max(diff(apply(g.sim, 1, range)))
}

random_marks <- function(x)  {
  x$marks <- sample(x$marks)
  x
}

bootstrap <- replicate(n, max_diff_env(random_marks(brambles)))
diff.obs <- max_diff_env(brambles)

#+ boot-dens, fig.height = 3, fig.width = 8, tidy = F
qplot(bootstrap, geom = "density") + 
  geom_vline(xintercept = diff.obs, col = 2, linetype = 2)

p.val <- p_2sided(bootstrap, diff.obs) 
format_p(p.val, n)

#' As the p-value is larger than $\alpha=0.05$ we cannot reject the null hypthesis.

