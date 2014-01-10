#+ setup, include = FALSE
opts_chunk$set(fig.path = 'figure/lab3-')
options(show.signif.stars = F)

#' # Lab 3: Bayesian statistics and WinBUGS

library(ggplot2)

#' ## Question 1: Diagnosing a patient
#' ### 1a) exact solution
#' `A patient comes to test for a rare disease (only c.1 in 100 000 people in the general population has it.) The test has sensitivity (probalility of true positiv of 0.99) and specificity (probability of true negative of 0.95) and gives a positive result. Given the positive result of the test, what is the probability that the patient has the disease? What if the disease is common (c.1 in 10 people in the general population has the disease)?`  
#' 

cond_prob <- function(prevalence, sens, spec) {
  prevalence * sens / ( prevalence * sens +  (1-prevalence) * (1-spec))
}

tmp <- cond_prob(prevalence=c(1/100000, 1/100, 1/10), sens=0.99, spec=0.95)
tmp

prev <- seq(1e-3, 1, 1e-3)
prob <- cond_prob(prevalence=prev, sens=0.99, spec=0.95)
#+ prev, fig.width=8, fig.height=3, tidy=F
qplot(prev, prob, geom = "line") + 
  labs(x = "population prevalence", 
       y = "probability of a patient\nhaving the disease")


#' ### 1b) Results using Winbugs
#' `Repeat the analysis using WinBUGS.`  
#' The trace generated with WinBUGS is read and the probability of the parameter $\theta$ (patient having the disease when diagnosed positively) calculated.
prob <- read.table("winbugs/trace", sep="\t", header=F)
names(prob) <- c("index", "theta")
prob$theta <- as.logical(prob$theta)
mean(prob$theta)
#' The probability, using the whole trace, is similar to the one obtained in 1a).  

#' Now we plot the change of probability with length of simulation. The horizontal dashed lines shows the exact solution. We see, that it took about 5000 interations for the function to converge.
prob$prob <- sapply(seq_len(nrow(prob)), function(x) mean(prob$theta[1:x]))
#+ trace, fig.width=8, fig.height=3, tidy=F
ggplot(prob, aes(index, prob)) + geom_line() + 
  coord_cartesian(ylim = c(0.1, 0.25)) +
  geom_hline(yintercept = tmp[2], col = 2, linetype = 2) + 
  labs(x = "number of iterations", 
       y = "probability of a patient\nhaving the disease")


