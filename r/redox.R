# todo
# identify and handle (exclude) outlier (single linkage)
# add parameter orp
# quantify quality of classificattion (T- F-values)

library(knitr)
library(plyr)
#opts_chunk$set(fig.path='figure/redox-')
library(reshape2)
library(lubridate)
library(ggdendro)
library(ggplot2)
library(OpenStreetMap)
library(sig)

# for box cox transformations
library(car)

# for parcoord() plot
library(MASS)

load("monitoring.RData")

load("d:/da/r/station.RData")
new <- structure(list(lab = c("90", "91", "2001", "O2", "O7", "1A", 
                              "XXIXu", "ND4.5", "73", "P4", "RF1", "RF2", "O2Damm", "1931,690", 
                              "Uferstr"), 
                      x = c(17348L, 17492L, 22339L, 27059L, 25914L, 
                            14673L, 15292L, 11775L, 17404L, 33307L, 28364L, 24822L, 27315L, 
                            4193L, 27059L), 
                      y = c(335106L, 335163L, 333117L, 332509L, 332580L, 
                            336076L, 338192L, 337628L, 334450L, 330672L, 331771L, 332622L, 
                            332959L, 345242L, 332509L)), 
                 .Names = c("lab", "x", "y"), 
                 row.names = c(1L, 2L, 3L, 4L, 5L, 6L, 8L, 9L, 10L, 
                               11L, 12L, 13L, 14L, 15L, 16L), 
                 class = "data.frame")

new <- data.frame(new, "internal_id"=NA, "z"=NA, "file"=NA, "comment"=NA, "type"=NA)

station <- rbind(station, new)

lab.gw <- c("LSA27", "LSA23", "LSA66",
            "LRD29E", "LSD30A", "LRD41", "LSA73", 
            "LSG32", "LSG9", "LSG11")

lab.ofl <- c("1", "12", "1A", "1B", "34", "42", "48", "51", "60", "73", "87", "90", "91", "EW-M")

p <- showStation(x=c(lab.gw, lab.ofl, grep("^HFB", station$lab, value=T))) +
  scale_x_continuous(limits=c(14250, 18000), minor_breaks=NULL, breaks=NULL) + 
  scale_y_continuous(limits=c(333500, 336700), minor_breaks=NULL, breaks=NULL) + 
  labs(x="", y="")

ggsave("figure/map_monitoring.pdf", height=21, width=29.7, unit="cm", plot=p)


showStation <- function (x, map = NULL, tbl = station) {
  if (is.null(map)) {
    load("osm-map.RData")
    map <- autoplot(map)
  }
  ind <- match(x, tbl$lab)
  map + geom_point(data = tbl[ind, ]) + geom_text(aes(label = lab), 
                                                  data = tbl[ind, ], size = 4, vjust = -1)
}

showStation(c(lab.ofl, lab.gw)) + xlim()


monitoring$lab <- factor(monitoring$lab, levels=lab.gw)
monitoring.raw$lab <- factor(monitoring.raw$lab, levels=lab.gw)

monitoring$is.bdl[is.na(monitoring$is.bdl)] <- ""
monitoring.raw$is.bdl[is.na(monitoring.raw$is.bdl)] <- ""

library(sp)
library(maptools)
library(rgeos)
point <- SpatialPoints(na.omit(station[match(lab.gw, station$lab), c("x", "y")]))
alt <- readShapeLines("N:/H811/WV-IWS.423/1_MA45_Untere_Lobau/gauster/gis/altarm.shp")
donau <- readShapeLines("N:/H811/WV-IWS.423/1_MA45_Untere_Lobau/gauster/gis/donau.shp")


dists <- data.frame(alt=numeric(), donau=numeric())
for (i in seq_along(point)) {
  dists[i, "alt"]   <- gDistance(point[i,], alt)
  dists[i, "donau"] <- gDistance(point[i,], donau)
}

dists$ofg <- apply(dists, 1, min)
dists <- data.frame(lab.gw, round(dists, 0))


# plot(donau)
# lines(alt)
# points(point, pch=3)
# axis(1, line=-2)
# axis(2, line=-2)
# 
# for (i in seq_along(point)) {
#   draw.circle(x=point$x[i], point$y[i], radius=as.numeric(dists[i, ]), border=2:3)
# }
# 


load("osm-map.RData")
redox <- droplevels(subset(monitoring, abbrev %in% c("O2", "NO3", "Mn", "Fe", "SO4", "NH4", "Temp") & is.bdl == ""))

# percentage of samples below detection limit 
ddply(redox, .(abbrev), function(x) round(sum(x$is.bdl=="<")/nrow(x), 3))

# argument must be strictly positive for transformation
# therefore even zero is bad... 
redox <- within(redox, value[abbrev=="O2" & value==0] <- NA)

wide <- dcast(data=redox, lab + time ~ abbrev, value="value")
#wide <- na.omit(wide)
wide <- merge(wide, dists, all.x=T)
wide$year <- year(wide$time)
write.csv(wide, file="redox.csv")

### Transformation ----
# estimate the exponents of box-cox transformation
bc <- apply(wide[, setdiff(names(wide), c("time", "lab", "year"))], 2, powerTransform)
power <- sapply(bc, function(x) as.numeric(x$lambda))

# test if power transformation is necessary
# if 0 or 1 are within the confidence interval (est +/- se *1.96)
# see car:::summary.powerTransform
sapply(bc, function(x) testTransform(x, lambda=1)$pval)
sapply(bc, function(x) testTransform(x, lambda=0)$pval)
            
# transform and scale the data
trans <- bcPower(wide[, setdiff(names(wide), c("time", "lab", "year"))], power)
colnames(trans) <- sub("\\^.*$", "", colnames(trans))
trans <- apply(trans, 2, scale)

# reappend labels
trans <- data.frame(wide[, c("lab", "time")], trans)




### Hierachical Clustering ----
tmp <- subset(trans, year(time) == 2011)

for (method in c("complete", "ward")) {
  hc <- hclust(dist(tmp[, setdiff(names(tmp), c("time", "lab"))]), method=method)
  p <- ggdendrogram(hc, theme_dendro = TRUE, rotate = TRUE) + 
    scale_x_continuous(labels=NULL) + 
    labs(x="", y="", title=paste("linkage:", method)) 
  
  plot(p)
  
  group <- as.factor(cutree(hc, 4))
  ct <- data.frame(month=month(tmp$time, label=T), lab=tmp$lab, group)
  ct <- merge(ct, station[, c("lab", "x", "y")])
  
  
  message(method)
  dcast(ct, month ~ lab, value.var="group")
  p <- ggplot(ct, aes(x=lab, y=month, fill=group, label=group)) + geom_raster() +
    geom_text() + 
    labs(x="", y="", title=paste("linkage:", method)) + 
    theme(legend.position="none")
  plot(p)
  
  p <- ggplot(ct, aes(x=lab, y=month, col=group, group=lab)) + 
    geom_line(size=10) + coord_flip() + 
    labs(x="", y="", title=paste("linkage:", method))
  plot(p)
  
  x <- melt(ct, id.vars=c("time", "lab", "group"), variable.name="abbrev")
  ggplot(x, aes(x=as.factor(group), y=value)) + geom_boxplot() + 
    facet_wrap(~abbrev, scales="free_y") +
    labs(x="", y="")
  
  
  
  #parcoord(wide[, -(1:2)], col=group, main="raw data", var.label=T)
  #parcoord(trans[, -(1:2)], col=group, main="transformed data")
}

# parallel coordinates
# centers muss man nur bei hierachical clustering extra berechnen
library(plyr)
centers <- round(ddply(ct[, setdiff(names(ct), c("time", "lab"))], .(group), colMeans), 3)
centers <- centers[, setdiff(names(centers), "group")]
ord <- apply(centers, 2, order)
ord <- do.call(order, as.data.frame(t(ord)))
centers <- centers[, ord]

pdf("tmp/clust_alpha.pdf")
parcoord(wide[ind, setdiff(names(wide), c("time", "lab"))][, ord], adjustcolor(col=km$cluster, alpha.f=0.2), lwd=1, var.label=T)
parcoord(centers, col=seq_len(nrow(centers)), lwd=2, add=T)
dev.off()




# zeitliche Entwicklung
saveHTML(
for (i in levels(droplevels(ct$month))) {
point <- subset(ct, month==i)
p <- autoplot(map) + 
  geom_point(data=point, mapping=aes(x=x, y=y, col=group), size=5) + 
  scale_color_discrete(drop=FALSE) + 
  geom_text(data=point, mapping=aes(x=x, y=y, label=lab)) + 
  labs(x="", y="", title=i)
plot(p)
})


### Partitioning Clustering ----
ind <- year(trans$time) == 2011
km <- kmeans(trans[ind, setdiff(names(trans), c("time", "lab"))], centers=4)
ct <- data.frame(wide[ind, ], group=km$cluster)

# boxplot
x <- melt(ct, id.vars=c("time", "lab", "group"), variable.name="abbrev")
ggplot(x, aes(x=as.factor(group), y=value)) + geom_boxplot() + 
  facet_wrap(~abbrev, scales="free_y") +
  labs(x="", y="")



### Fuzzy Clustering ----
ind <- year(trans$time) == 2011
library(e1071)
# set starting values
# ufcl depends very much on starting values, really?
cm  <- cmeans(trans[ind, setdiff(names(trans), c("time", "lab"))], centers=4, method="cmean", m=1.75)

tmp <- wide[ind, ]
tmp <- data.frame(month=month(tmp$time, label=T), 
                         lab=tmp$lab, group=cm$membership)
tmp <- melt(tmp, id.vars=c("month", "lab"), variable.name="group")
levels(tmp$group) <- 1:4


# simplified map
b <- 2
a1 <- b / 10

x <- matrix(ncol=2, nrow=3)
x <- (row(pos) + col(pos))*(b+a1) + row(pos)*(b+2*a1)
x <- x - x[1, 1]

y <- seq(from=1, by=b, length.out=4)
y <- y + y*a1
y <- y - y[1]

pos <- data.frame(lab = c("LSA27", "LSA23", "LSA66",
                          "LRD29E", "LSD30A", "LRD41", "LSA73", 
                          "LSG32", "LSG9", "LSG11"),
                  x1 = c(1, 1, 1, 2, 2, 2, 2, 3, 3, 3),
                  x2 = c(1, 1, 2, 1, 1, 2, 2, 1, 1, 1),
                  y  = c(2, 3, 2, 2, 3, 3, 4, 1, 2, 4), 
                  stringsAsFactors=F)

coord <- data.frame(lab = pos$lab,
                    x   = x[cbind(pos$x1, pos$x2)],
                    y   = y[pos$y])

hfb <- data.frame(x=c(mean(coord$x[c(1, 3)]),
                      mean(coord$x[c(4, 6)]),
                      coord$x[10]),
                  y=y[c(1, 3, 3)])


rect <- with(coord, data.frame(xmax=x+b/2, xmin=x-b/2, ymax=y+b/2, ymin=y-b/2))
offset <- expand.grid(dx=c(-1, 1)*b/4, dy=c(1, -1)*b/4)

membership <- cbind(tmp, coord[match(tmp$lab, coord$lab), c("x", "y")])
membership$x1 <- membership$x + offset$dx[membership$group]
membership$y1 <- membership$y + offset$dy[membership$group]


p <- ggplot(membership[membership$value >= 0.2, ], 
            aes(x=x1, y=y1, col=group, size=value)) + 
  geom_point() + scale_size_area("membership", max_size=6, limits=c(0.2, 1)) + 
  annotate(geom="rect", col="gray", fill=NA, 
           xmin=rect$xmin, xmax=rect$xmax, ymin=rect$ymin, ymax=rect$ymax) + 
  annotate(geom="point", x=hfb$x, y=hfb$y, col="darkgrey", shape=3, size=2) + 
  labs(x="", y="") + 
  scale_x_continuous(breaks=hfb$x, 
                     labels=c("Kreuzgrund", "Rohrwörth", "Gänshaufen")) + 
  scale_y_continuous(breaks=c(y[1], tail(y, 1)),
                     labels=c("donaunahe", "altarmnahe")) +
  coord_fixed() + 
  facet_wrap(~month, nrow=3, drop=F) + 
  theme_bw() +
  theme(axis.ticks=element_line(colour=NA), panel.grid.major=element_line(colour=NA))
p
ggsave("tmp/fuzzy_cluster_simplified.pdf", width=29.7, height=21, unit="cm", 
       scale=1.25, plot=p) 

## map, drawn to scale
coord <- subset(station, lab %in% unique(membership$lab), 
                select=c("lab", "x", "y"))

b <- 200
rect <- with(coord, data.frame(xmax=x+b/2, xmin=x-b/2, ymax=y+b/2, ymin=y-b/2))
offset <- expand.grid(dx=c(-1, 1)*b/4, dy=c(1, -1)*b/4)

membership <- cbind(tmp, coord[match(tmp$lab, coord$lab), c("x", "y")])
membership$x1 <- membership$x + offset$dx[membership$group]
membership$y1 <- membership$y + offset$dy[membership$group]

p <- autoplot(map) + 
  annotate(geom="rect", col="gray", fill="white", 
           xmin=rect$xmin, xmax=rect$xmax, ymin=rect$ymin, ymax=rect$ymax) + 
  geom_point(data=membership[membership$value >= 0.2, ], 
         mapping=aes(x=x1, y=y1, col=group, size=value)) + 
   scale_size_area("membership", max_size=6, limits=c(0.2, 1)) + 
  
  scale_x_continuous(limits=c(14500, 18000), minor_breaks=NULL, breaks=NULL) + 
  scale_y_continuous(limits=c(333500, 336200), minor_breaks=NULL, breaks=NULL) +
  facet_wrap(~month, nrow=3, drop=F)

ggsave("tmp/fuzzy_cluster.pdf", width=50, height=29.7, unit="cm", plot=p) 




### PCA ----
pc <- prcomp(subset(trans, year(time) == 2011)[, setdiff(names(trans), c("time", "lab"))], scale=T)
biplot(pc, choice=1:2)
summary(pc)

#plot(hclust(dist(t(predict(pc))), method="single"))

# cluster the PC
ind <- year(trans$time) == 2011
km <- kmeans(predict(pc)[, 1:2], centers=4)
ct <- data.frame(wide[ind, ], group=km$cluster)

# boxplot
x <- melt(ct, id.vars=c("time", "lab", "group"), variable.name="abbrev")
ggplot(x, aes(x=as.factor(group), y=value)) + geom_boxplot() + 
  facet_wrap(~abbrev, scales="free_y") +
  labs(x="", y="")

# cluster the variables
plot(hclust(d=dist(t(trans[, setdiff(names(trans), c("time", "lab"))]))))

# Regressionskoeffizienten
round(cor(wide[, -c(1:2)]), 2)




# Einhüllende Mn/SO4
tmp <- wide[wide$lab != "LSA73" & wide$Mn > 10, ]
tmp1 <- tmp[tmp$Mn > 22, ]


p1 <- tmp1[which.min(tmp1$SO4), c("SO4", "Mn")]
p <- tmp1[-which.min(tmp1$SO4), c("SO4", "Mn")]
k <- max((p$Mn - p1$Mn)/(p$SO4 - p1$SO4))
d <- p1$Mn - k*p1$SO4


ggplot(tmp, aes(SO4, Mn, color=factor(year))) + geom_point() +
  geom_abline(intercept=d, slope=k)


tmp <- subset(monitoring.raw, is.bdl  != "<")
tmp1 <- dcast(tmp, time  + lab ~ abbrev, value.var="value")
ggplot(tmp1, aes(NO3, NO2)) + geom_point()
ggplot(tmp1, aes(O2, NH4)) + geom_point()


tmp1 <- subset(tmp, abbrev %in% c("NO2", "NO3", "NH4", "O2"))
#tmp1 <- subset(tmp, grepl("N-", abbrev))

#tmp1 <- ddply(tmp1, .(abbrev), function(x)  data.frame(x, value.s=scale(x$value)))
#tmp1$abbrev <- as.character(tmp1$abbrev)


library(sig)

conv.par <- read.csv2("N:/H811/WV-IWS.423/2 MA31 Donauinsel Nord/gauster/R/conv-parameter.def", 
                      stringsAsFactors=F)
conv.unit <- read.csv2("N:/H811/WV-IWS.423/2 MA31 Donauinsel Nord/gauster/R/conv-unit.def", 
                       stringsAsFactor=F)

table_rev <- data.frame(from=conv.par$to, to=conv.par$from, factor=1/conv.par$factor, stringsAsFactors=F)

tmp2 <- fixParameter(x=tmp1, to=table_rev$to, table=table_rev)
tmp2 <- fixUnit(x=tmp2, table=conv.unit, good=tvo[, c("abbrev", "unit")])

ggplot(subset(tmp2, abbrev!="O2"),
              aes(value, abbrev, group=time)) + facet_wrap(~lab) + 
  geom_line(col="darkgrey")  + geom_jitter(mapping=aes(col=month(time)))
