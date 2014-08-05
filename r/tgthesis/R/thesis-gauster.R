require(ggplot2)
require(grid)
require(plyr)
require(reshape2)
require(lubridate)
require(OpenStreetMap)

width <- 408.297 * 0.352777778

mytheme <- theme_bw(base_size = 10) + 
  theme(legend.background = element_blank(), legend.key = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(), 
        strip.background = element_blank(), plot.background = element_blank(), 
        axis.line = element_blank(), panel.grid = element_blank(),
        axis.ticks = element_blank())


