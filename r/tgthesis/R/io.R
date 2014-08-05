# the whole monitoring data set ----
load("../../data//monitoring.RData")
load("../../data//osm-map.RData")

lab.gw <- c("LSA27", "LSA23", "LSA66",
            "LRD29E", "LSD30A", "LRD41", "LSA73", 
            "LSG32", "LSG9", "LSG11")

lab.ofl <- c("1", "12", "1A", "1B", "34", "42", "48", "51", "60", "73", "87", "90", "91", "EW-M")

monitoring$lab <- factor(monitoring$lab, levels=lab.gw)
monitoring.raw$lab <- factor(monitoring.raw$lab, levels=lab.gw)

monitoring$is.bdl[is.na(monitoring$is.bdl)] <- ""
monitoring.raw$is.bdl[is.na(monitoring.raw$is.bdl)] <- ""



# just a subset of parameters influencing orp ----
redox <- droplevels(subset(monitoring, abbrev %in% c("O2", "NO3", "Mn", "Fe", "SO4", "NH4", "Temp")))


# definition of monitors 
load("../../data//station.RData")
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


#require(sig)
# p <- showStation(x=c(lab.gw, lab.ofl, grep("^HFB", station$lab, value=T)), 
#                  map=autoplot(map)) +
#   scale_x_continuous(limits=c(14250, 18000), minor_breaks=NULL, breaks=NULL) + 
#   scale_y_continuous(limits=c(333500, 336700), minor_breaks=NULL, breaks=NULL) + 
#   labs(x="", y="")
# 
# ggsave("../../thesis/figure/map_monitoring.pdf", height=21, width=29.7, unit="cm", plot=p)



require(sp)
require(maptools)
require(rgeos)
point <- SpatialPoints(na.omit(station[match(lab.gw, station$lab), c("x", "y")]))
alt <- readShapeLines("../../data//altarm.shp")
donau <- readShapeLines("../../data//donau.shp")


dists <- data.frame(alt=numeric(), donau=numeric())
for (i in seq_along(point)) {
  dists[i, "alt"]   <- gDistance(point[i,], alt)
  dists[i, "donau"] <- gDistance(point[i,], donau)
}

dists$ofg <- apply(dists, 1, min)
dists <- data.frame(lab=lab.gw, round(dists, 0))

station <- merge(station, dists, all.x=T)
rm(new, lab.gw, lab.ofl, point, alt, donau, dists)