library(ggplot2)
library(ggpubr)
library(readr)
library(reshape2)
library(tibble)
library(svglite)
library(stringr) 
library(viridis)
library(dplyr)
library(forcats)

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
my_colours <- my_colours[c(3, 1, 1)] # take only colour for growth and growth+NGAM
my_alpha <- c(0.8, 0.8, 0.4)

# create directories for the plots
dir.create("3A")

# read data
path <- "../fluxes.txt"
data <- read_delim(path, "\t", col_names = TRUE, comment = "#")
data <- subset(data, ((`objective 1` == "growth (MAX)") &
                 (is.na(`objective 2`) | (`objective 2` == "NGAM (MAX)")) &
                 (parsimonious == FALSE)) | 
                 ((`objective 1` == "growth (MAX)") & is.na(`objective 2`) &
                  (parsimonious == TRUE)))

# add extra column with summarised objective
data$objective <- paste(data$`objective 1`, " + ", data$`objective 2`, " (",
                          data$parsimonious, ")")

# use absolute value of fluxes, directionality doesn't matter here
# and set fluxes that are 0 too small value
data$`average flux` <- abs(data$`average flux`)
data$`average flux`[data$`average flux` <= 1e-9] <- 1e-9

# order data 
data$objective <- factor(data$objective, levels = unique(data$objective)[c(2, 3, 1)])

# divide data set into purely maximal growth and +parsimonious or +NGAM
dataBase <- data[(data$objective == "growth (MAX)  +  NA  ( FALSE )"), ]
dataBase <- rbind(dataBase, dataBase)
dataMod <- data[!(data$objective == "growth (MAX)  +  NA  ( FALSE )"), ]
diff <- dataMod[ , c(1:10, 13)]

# calculate absolute and relative differences
diff$absDiff <- dataMod$`average flux` - dataBase$`average flux`
diff$relDiff <- diff$absDiff / dataBase$`average flux`

diff$relDiff[diff$relDiff > 2] <- 2
diff$relDiff[diff$relDiff < -1] <- -1

# plot relative difference in phase I
fig <- ggplot(subset(diff, phase == "phase I"), 
              aes(x = pathway, y = relDiff, 
                  fill = objective)) +
  my_theme + 
  geom_hline(yintercept = 0, color = "black") +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_point(aes(color = objective), position = position_jitterdodge(),
             size = 0.3) +
  xlab("pathway") +
  ylab(expression("relative changes of fluxes")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 12)) +
  scale_x_discrete(breaks = levels(as.factor(data$pathway)),
                   labels = function(x) str_wrap(x, width = 30)) +
  theme(aspect.ratio = 0.4) +
  scale_color_manual(values = my_colours[c(1, 2)]) +
  scale_fill_manual(values = my_colours[c(1, 2)]) +
  guides(color = guide_legend(nrow = 3))
ggsave(file = "3A/relNGAM_I.svg", plot = fig)

# plot relative difference in phase II
fig <- ggplot(subset(diff, phase == "phase II"), 
              aes(x = pathway, y = relDiff, 
                  fill = objective)) +
  my_theme + 
  geom_hline(yintercept = 0, color = "black") +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_point(aes(color = objective), position = position_jitterdodge(),
             size = 0.3) +
  xlab("pathway") +
  ylab(expression("relative changes of fluxes")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 12)) +
  scale_x_discrete(breaks = levels(as.factor(data$pathway)),
                   labels = function(x) str_wrap(x, width = 30)) +
  theme(aspect.ratio = 0.4) +
  scale_color_manual(values = my_colours[c(1, 2)]) +
  scale_fill_manual(values = my_colours[c(1, 2)]) +
  guides(color = guide_legend(nrow = 3))
ggsave(file = "3A/relNGAM_II.svg", plot = fig)
