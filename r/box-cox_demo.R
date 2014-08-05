# for box cox transformations
library(car)


# argument must be strictly positive for transformation
# therefore even zero is bad... 
redox <- within(redox, value[abbrev=="O2" & value==0] <- NA)

o2 <- subset(redox, abbrev=="O2", select="value", drop=T)
o2t <- na.omit(bcPower(o2, coef(powerTransform(o2))))

dens <- list(o2 = density(o2, na.rm=T, from=0),
             o2t= density(o2t, na.rm=T))

old.par <- par(no.readonly = T)

x <- seq(mean(o2t) - 3*sd(o2t), mean(o2t) + 3*sd(o2t),length=1000)

dens.df <- lapply(dens, function(x) data.frame(data=x$data.name, x=x$x, y=x$y))
dens.df <- do.call(rbind, dens.df)
dens.df <- rbind(dens.df, 
                 data.frame(data="norm", x=x, y=dnorm(x=x, mean=mean(o2t), sd=sd(o2t))))

mylabs <- list("Messwerte in mg/l",
               bquote(Box-Cox:~lambda==.(round(coef(powerTransform(o2)), 2))),
               bquote(NV:~mu==.(paste(round(mean(o2t), 2), ",", sep=""))~sigma==.(round(sd(o2t), 2))))


p <- ggplot(dens.df, aes(x, y, linetype=data, col=data)) + geom_line() +
  scale_color_discrete("", labels = mylabs) + 
  scale_linetype_discrete("", labels = mylabs) +
  labs(x="", y="Dichte") + 
  mytheme + 
  theme(axis.text.y=element_blank(),
        axis.title.x=element_blank(), plot.title=element_blank(),
        plot.margin=unit(x=c(0, -1, -0.5, 0), unit="line"))

ggsave("figure/box-cox_dichte.pdf", plot=p, 
       width=width, height=30, unit="mm")


tmp <- rbind(data.frame(x=o2,  data="o2"), 
             data.frame(x=o2t, data="o2t"))

mylabs <- list("Messwerte",
               bquote(Box-Cox:~lambda==.(round(coef(powerTransform(o2)), 2))))
grid <- data.frame(a=-10:20, b=1)

p <- ggplot(tmp, aes(sample = x, col=data, shape=data, group=data)) + 
  geom_abline(data=grid, mapping=aes(intercept=a, slope=b), 
              col="lightgrey", size=0.1) + 
  stat_qq() + 
  scale_color_discrete("", labels = mylabs) + 
  scale_shape_discrete("", labels = mylabs) +
  labs(x="theoretisches Quantil", y="Stichprobenquantil") + 
  mytheme +
  theme(plot.title=element_blank(),
        plot.margin=unit(x=c(0, -1, 0, 0), unit="line"))

ggsave("figure/box-cox_qq.pdf", plot=p, 
       width=width, height=60, unit="mm")




### Transformation ----
# estimate the exponents of box-cox transformation
power <- apply(wide[, setdiff(names(wide), c("time", "lab", "year"))], 2, powerTransform)
scale.par <- data.frame(power = sapply(power, function(x) as.numeric(x$lambda)))

# test if power transformation is necessary
# if 0 or 1 are within the confidence interval (est +/- se *1.96)
# see car:::summary.powerTransform
sapply(power, function(x) testTransform(x, lambda=1)$pval) > 0.05
sapply(power, function(x) testTransform(x, lambda=0)$pval) > 0.05

# transform and scale the data
trans <- bcPower(wide[, setdiff(names(wide), c("time", "lab", "year"))], scale.par$power)
colnames(trans) <- sub("\\^.*$", "", colnames(trans))
scale.par$mean <- colMeans(trans, na.rm=T)
scale.par$sd   <- sapply(trans, sd, na.rm=T)

trans <- apply(trans, 2, scale)
