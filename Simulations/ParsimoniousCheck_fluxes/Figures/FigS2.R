library(ggplot2)
library(ggpubr)
library(readr)
library(reshape2)
library(tibble)
library(svglite)
library(stringr) 
library(viridis)
library(dplyr)

# set theme
my_theme <- theme_light(base_size = 20) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title = element_text(size = 20)) + 
  theme(rect = element_rect(fill = "transparent")) +
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14))

my_colours <- viridis_pal()(7)
my_colours <- my_colours[c(1, 2, 4, 5, 6, 7)] # exclude one colour for growth+NGAM

# create directories for the plots
dir.create("S2")

# read data
path <- "../fluxes.txt"
data <- read_delim(path, "\t", col_names = TRUE, comment = "#")
data <- subset(data, (`objective 2` != "NGAM (MAX)") | is.na(`objective 2`))

# add extra column with summarised objective
data$objective <- paste(data$`objective 1`, " + ", data$`objective 2`)

# use absolute value of fluxes, directionality doesn't matter here
# and set fluxes that are 0 too small value
data$`average flux` <- abs(data$`average flux`)
data$`average flux`[data$`average flux` <= 1e-9] <- 1e-9

# order data 
data$objective <- factor(data$objective, levels = unique(data$objective))

# divide data set according to parsimonious
dataFalse <- data[data$parsimonious == FALSE, ]
dataTrue <- data[data$parsimonious == TRUE, ]

# check if the order is the same in both subsets
# identical(dataTrue[c(1,2,4,5,6,9,10)], dataFalse[c(1,2,4,5,6,9,10)])

# calculate absolute and relative differences
diff <- dataTrue[ , c(1:10, 13)]
diff$absDiff <- dataTrue$`average flux` - dataFalse$`average flux`

# plot absolute flux difference for phase I
fig <- ggplot(data = subset(diff, phase == "phase I"),
              aes(x = pathway, y = absDiff,
                                 fill = objective, color = objective)) +
  my_theme + 
  geom_hline(yintercept = 0, color = "black") +
  geom_boxplot(fill = "black", color = "black", alpha = 0.2,
               outlier.shape = NA) +
  geom_jitter(alpha = 0.5, size = 0.5, width = 0.3) +
  facet_grid(objective ~ .) +
  xlab("pathway") +
  ylab(expression("absolute flux difference [mmol (gDW h)"^"-1"~", h"^"-1"~"]")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 12)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
  scale_color_manual(values = my_colours) +
  scale_fill_manual(values = my_colours) +
  theme(legend.position =  "none")
ggsave(file = "S2/absDiff_I.svg", plot = fig)

# plot absolute flux difference for phase II
fig <- ggplot(data = subset(diff, phase == "phase II"),
              aes(x = pathway, y = absDiff, fill = objective, color = objective)) +
  my_theme + 
  geom_hline(yintercept = 0, color = "black") +
  geom_boxplot(fill = "black", color = "black", alpha = 0.2,
               outlier.shape = NA) +
  geom_jitter(alpha = 0.5, size = 0.5, width = 0.3) +
  facet_grid(objective ~ .) +
  xlab("pathway") +
  ylab(expression("absolute flux difference [mmol (gDW h)"^"-1"~", h"^"-1"~"]")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 12)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
  scale_color_manual(values = my_colours) +
  scale_fill_manual(values = my_colours) +
  theme(legend.position =  "none")
ggsave(file = "S2/absDiff_II.svg", plot = fig)