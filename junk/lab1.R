#+ setup, include = FALSE
opts_chunk$set(fig.path = 'figure/lab1-')
options(show.signif.stars = F)

#' # Lab 1: Introduction to R and review of some useful statistical techniques

#' First, load the packages,
library(ggplot2)
library(reshape2)

#' do some global assignments,
n <- 1000    # number of iterations used for Monte Carlo and bootstrapping

#' and define general useful functions.
format_p <- function(p, n) ifelse(p < (1/n),  paste("<", 1/n, sep=""), p)
p_2sided <- function(sim, obs) mean(abs(sim-mean(sim)) >= abs(obs - mean(sim)))

#' ## Question 1: A binomial proportion test vs. Monte Carlo test
#' `Suppose you have tossed a coin 20 times and observed 'heads' 15 times. You want to do a statistical test to see whether your coin is fair.`
#' ### 1a) hypothesises and assumptions
#' `Write down explicitly the null and the alternative hypothesis. Are there any assumptions to make?`
#' 

#' The null hypothesis ($H_0$) states that the coin is fair, giving on average 
#' 50% of the tosses head respectively tails (probability $p  =  0.5$). We are 
#' going to accept the null hyptothesis if the p-value of the test is less than 
#' our chosen significance level of $\alpha = 0.05$.  
#' 
#' $H_0: p  =  0.5$  
#' $H_1: p \ne  0.5$  
#' 
#' The following is assumed:
#' * data follow binomial distribution
#' * all the tosses are independent
#'    

#' ### 1b) exact test 
#' `Use a binomial proportions test to test the hypothesis in 1a). How would you interpret the result?`
#+ tidy = F  
heads.obs <- 15
tosses <- 20
binom.test(x = heads.obs, n = tosses, p = 0.5, 
           alternative = "two.sided", conf.level = 0.95)

#' As the p-value is less than 0.05, so we reject the null hypothesis. Therefore 
#' we can say that for the choosen significance level of $0.05$ the coin is not 
#' fair. 
#' The test output furthermore states, that the empirical probability of success
#' is $\frac{15}{20}  =  0.75$. We also get a confidence interval stating that in at
#' least 95% of the cases $p$ is within the interval $[0.51, 0.94]$
#' (which excludes $p  =  0.50$).
#'

#' ### 1c) Monte Carlo
#' `Use a Monte Carlo test to test the hypothesis in 1a). Write down your algorithm before writing the program. Report and interpret the result.`

#' * get a new sample $n$ times from a binomial distribution, 
#' * for each sample: count the heads and determine how 'extreme' this result 
#' is, 
#' * compute p-value (mean of the number of extreme samples)

set.seed(2718)
heads <- rbinom(n = n, size = tosses, prob = 0.5)
#+ tidy = F
sim <- data.frame(toss.no = 1:n,  
                  heads = factor(heads, levels = seq_len(tosses)),  
                  extreme = (heads >=  heads.obs | heads <=  tosses-heads.obs))

#+ echo = 2
options(width = 85)
table(sim$heads)
options(width = 75)

#' In the following plot frequencies of the number of heads for each sample
#' (20 tosses) is shown. The bins of extreme events (heads observed more than 14
#' times or less than 6 times) are colored. 
#+ hist_mc, fig.height = 4, fig.width = 8, tidy = F
ggplot(sim, aes(heads, fill = extreme)) +  
  geom_bar(binwidth = 1, col = 1, alpha = 0.4, show_guide = F) + 
  scale_fill_manual(name = "", values = c(NA, 2)) +
  scale_x_discrete(drop = F) +
  labs(x = "number of heads in one sample", y = "count of samples")

#' The p-value is the 'read areas' share of the 'total area'.
p.val <- mean(sim$extreme)
format_p(p.val, n)


#' To see how the p-value changes with the number of samples drawn in the Monte
#' Carlo simulation, we plot the development of the p-value against the sample 
#' size. The dashed horizontal line indicates our choosen significance level of 
#' $\alpha=0.05$.
sim$p.vals <- sapply(sim$toss.no, function(x) mean(sim$extreme[1:x]))
#+ mc-p_value, fig.width = 8, fig.height = 4, tidy = F
ggplot(sim, aes(toss.no, p.vals)) + 
  geom_line() + 
  geom_hline(yintercept = 0.05, col = 2, linetype = 2) + 
  labs(x = "number of samples in Monte Carlo simulation", y = "p value") + 
  if(max(sim$p.vals > 0.25)) coord_cartesian(ylim = c(0, 0.25))

#' Aditionally compute quantiles which cover 95% of the data and check if 
#' observed probability of 0.75 is within this interval.
quants <- quantile(x = heads, probs = c(0.025, 0.975))
quants
as.vector(quants[1] < heads.obs & heads.obs < quants[2])
#' Again the p-value is less than 0.05 as also derived in 1b). As expected, the
#' observed probability **does not** lie inside the interval of the quantiles. 


#' ### 1d) bootstrap 
#' `Use bootstrap to test the hypothesis in 1a). Write down your algorithm before writing the programm. Report and interpret the result.`
#' 
#' * generate the bootstrap (resample our sample `n` times), 
#' * for each realisation: count the heads, 
#' * compute the probability of heads


#+ tidy = F
sample.obs <- c(rep(0, 5), rep(1, 15))
heads <- replicate(n, sum(sample(sample.obs, size = 20, replace = T)))
sim <- data.frame(toss.no = 1:n,  
                  heads = factor(heads, levels = seq_len(tosses)),  
                  extreme = (heads <=  tosses/2))

#+ echo = 2
options(width = 85)
table(sim$heads)
options(width = 75)

#' In the following plot frequencies of the number of heads for each resample
#' (20 tosses) is shown. The bins of extreme events (heads observed less than 11 
#' times) are colored. 
#+ hist_boot, fig.height = 4, fig.width = 8, tidy = F
ggplot(sim, aes(heads, fill = extreme)) +  
  geom_bar(binwidth = 1, col = 1, alpha = 0.4, show_guide = F) + 
  scale_fill_manual(name = "", values = c(NA, 2)) +
  scale_x_discrete(drop = F) +
  labs(x = "number of heads in one sample", y = "count of samples")

#' The p-value is the 'read areas' share of the 'total area'.
p.val <- mean(sim$extreme)
format_p(p.val, n)


#' To see how the p-value changes with the number of samples drawn in the 
#' boostrap simulation, we plot the development of the p-value against the sample 
#' size. The dashed horizontal line indicates our choosen significance level of 
#' $\alpha=0.05$.
sim$p.vals <- sapply(sim$toss.no, function(x) mean(sim$extreme[1:x]))
#+ boot-p_value, fig.width = 8, fig.height = 4, tidy = F
ggplot(sim, aes(toss.no, p.vals)) + 
  geom_line() + 
  geom_hline(yintercept = 0.05, col = 2, linetype = 2) + 
  labs(x = "number of samples in boostrap simulation", y = "p value") + 
  if(max(sim$p.vals > 0.25)) coord_cartesian(ylim = c(0, 0.25))



#' ### 1e) Compare/Contrast results
#' `Compare and contrast the results in 1b) - 1d).`  
#' All the tests (exact test, MC-test, bootstrap test) indicate that the coin is not fair regarding a significance level of $\alpha=0.05$. 



#' ## Question 2: A oneway ANOVA vs. a boostrap ANOVA
#' ### 2a) Examine the dataset
#' `Take a look at the iris dataset and type`
help(iris)

#' `to read the description of the dataset. It contains measurements of petals 
#' and sepals of 150 iris flowers belonging to 3 different species. We want to
#' know, whether sepal width of flowers is the same for all three species. Take
#' a loock at the data, do you think there is a difference?`

#+ boxplot_species, fig.width = 8, fig.height = 3
ggplot(iris, aes(Species, Sepal.Length)) + geom_boxplot() + coord_flip()

#' I'd assume there is a difference, because the medians are quite different and the boxes are quite narrow. 
#' 

#' ### 2b) parametric ANOVA
#' `In order to be valid, ANOVA requires normality and homoscedasticity of
#' residuals. Sometimes, these assumptions do not hold. In fact, take a look:`
lm.obs <- lm(Sepal.Length ~ Species, data = iris)
summary(lm.obs)
anova(lm.obs)

#+ anova-diag, echo=2
par(mfrow = c(2, 2))
plot(lm.obs)
par(mfrow = c(1, 1))

#' ### 2c) bootstrap
#' `Looks like there might be some heteroscedasticity. Let's use non-parametric bootstrap instead:`
#' * `sample 100 times with replacement from the original dataset`
#' * `so ANOVA F-test and record the F-statistic`
#' * `see how extreme your observed value is compared to the possible ones.`
#'    
 
#' hypothesis for F-value (one-sided, there is a difference between the species):  
#' $H_0: F_{obs}  \le  F_{sim}$  
#' $H_1: F_{obs}  >  F_{sim}$   


coeff <- as.list(round(coef(lm.obs), 3))
names(coeff) <-levels(iris$Species)
obs <- data.frame(f.value = anova(lm.obs)$F[1], coeff)
                  
lm.sim <- replicate(n, lm(Sepal.Length ~ sample(Species), data = iris), 
                    simplify=F)
coeff <- data.frame(round(do.call(rbind, lapply(lm.sim[1:20], coef)), 3))
names(coeff) <- levels(iris$Species)
sim <- data.frame(f.value = sapply(lm.sim, function(x) anova(x)$F[1]), 
                  coeff)

p.val <- mean(sim$f.value >= obs$f.value)
format_p(p.val, n)

#' The observed F-value in 2b) of 119 is about 15 times higher than the largest 
#' F-value in the bootstrap. As a consequence the p-value is smaller than $1/n$. The corresponding figure of the density estimation of the F-value can be found in the following section. 


#' ### 2d) compare parametric and bootstrap ANOVA
#' `Compare the species specific estimates of mean sepal length for parametric an for the bootstrap ANOVA. Write down the null and alternative hypotheses. Compare the results of the parametric test an the bootstrap test.`
#' 


#' hypothesis for the coefficients $c$ (two-sided, coefficients are different 
#' from zero):  
#' $H_0: c_{obs}  =  c_{sim}$  
#' $H_1: c_{obs}  \ne  c_{sim}$  


#+ density_anova, fig.height = 6, fig.width = 8, tidy = F
ggplot(melt(sim, id.vars = NULL, ), aes(value, xintercept = value)) + 
  geom_density() + 
  geom_vline(data = melt(obs, id.vars = NULL), 
             mapping = aes(xintercept = value),
             col = 2, linetype = 2) +
  facet_wrap(~variable, scales = "free") + 
  labs(x = "") 

#+ tidy=F
p.val <- sapply(levels(iris$Species), 
                function(x) p_2sided(sim=sim[, x], obs=obs[, x]))
format_p(p.val, n)

#' Because the estimates of the parameteric test (red, dashed line) are very different from the ones in the bootstrap (p-value less than $1/n$) we reject the null hypothesis. This means, there **are** differences for every species. This also corresponds to the parameteric ANOVA where we get very small p-values for every species. 


#' ### 2e) Monte Carlo
#' `Would it be possible to use Monte Carlo test. Why?`  
#'   

#' No. To perform a Monte Carlo test the distribution of the population has to
#' be known.


#' 
#' ## Question 3: Cross-validation of a regression model
#' `Suppose, we want to know, whether species is a better predictor of the sepal length, tha the sepal width is.`  

#' ### 3a) development of models
#' `Write down the two models we are comparing.`


#+ tidy=F
formulas <- list("Species"     = formula(Sepal.Length ~ Species),
                 "Sepal.Width" = formula(Sepal.Length ~ Sepal.Width))
formulas

#' ### 3b) comparison of models
#' `How can you compare their fit? Which is the better model?`  
#' 
#' A models goodness of fit can be specified for example by:  
#' * the absolute values of the residuals
#' * the sum of squared residuals  
#' * the coefficient of determination (i.e. $R^2$) 
#' * the F-statistic 
#' * an information criterion ($AIC$, $BIC$)  
#' 
models <- lapply(formulas, function(x) lm(x, data = iris))
lapply(models, summary)
lapply(models, AIC)

#' For the model `sepal.length ~ Species` we get a smaller $AIC$ and a higher
#' F-value as well as a higher coefficient of determination, stating that 
#' *Species* is a better predictor in our case. 

#' ### 3c) cross-validation
#' `Now use the cross-validation to compare the two models. Discuss the results. Would you generally expect a categorial variable to be a better predictor tha a continuous one? Are there any practical considerations to be taken into account?`

#+ tidy=F
resid.cv <- data.frame()
for (i in seq_along(iris$Sepal.Length)) {
  tmp <- sapply(formulas, 
                function(x) predict(object = lm(x, data = iris, subset = -i),
                                    newdata = iris[i, ]))
  tmp <- data.frame(index    = i, 
                    term     = names(formulas),
                    residual = iris$Sepal.Length[i] - tmp)
  resid.cv <- rbind(resid.cv, tmp)
}


#' Comparing cross-validation residuals with the residuals of the complete 
#' linear model

#+ tidy=F
resids <- function(x) {
  data.frame(index    = seq_along(fitted(x)),
             term     = attr(terms(x), "term.labels"),
             residual = resid(x))
}

comp <- rbind(data.frame(method = "cv", resid.cv), 
              data.frame(method = "lm",
                         do.call(rbind, lapply(models, resids))))

comp <- merge(comp, data.frame(index    = seq_len(nrow(iris)), 
                               observed = iris$Sepal.Length))

tapply(comp$residual, comp[, c("term", "method")], function(x) sum(x^2))

#' Looking at the sum of squared residuals, the model Sepal.Length ~ Species performs better. 
                 

#+ cv_resid, tidy=F, fig.width=8, fig.height=4
ggplot(comp, aes(observed, residual, col=method)) + 
  facet_wrap(~term) + geom_point() + labs(x="observed value")
#' Compared to the residuals from the linear models the absolute values of the cross-validation residuals are slightly larger. That is, because the predicted value itself is not used to fit the model.  
#' Generally one cannot say, that a contiuous variable is a better predictor than a categorial one. In our case, the model using the categorial variable *species* yields smaller residuals, because the variance within a group is relatively small compared to the overall variance. 
