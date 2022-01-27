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
dir.create("3BD")

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

# plot flux count
fig <- ggplot(data  %>%
                group_by(phase, objective, `flexibility 1`, `flexibility 2`) %>%
                summarize_at("average flux", sum, na.rm = T),
              aes(x = phase, y = `average flux`,
                  fill = objective, alpha = objective)) +
  my_theme + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = objective), position = position_jitterdodge()) +
  xlab("phase") +
  ylab(expression("sum of average fluxes [mmol (gDW h)"^"-1"~", h" ^"-1"~"]")) +
  scale_color_manual(values = my_colours) +
  scale_fill_manual(values = my_colours) +
  scale_alpha_manual(values = my_alpha) +
  scale_x_discrete(limits = c("phase II", "phase I")) +
  theme(aspect.ratio = 1) +
  guides(color = guide_legend(nrow = 3)) +
  coord_flip()
ggsave(file = "3BD/totalFlux.svg", plot = fig)

# plot average generation time (to avoid double counting, we use only phase I here)
fig <- ggplot(subset(data, phase == "phase I") %>%
                group_by(objective, `flexibility 1`, `flexibility 2`) %>%
                summarize_at("average generation time", mean, na.rm = T),
              aes(x = objective, y = `average generation time`,
                  fill = objective, alpha = objective)) +
  my_theme + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = objective), width = 0.2) +
  xlab("objective") +
  ylab("average generation time [h]") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
  scale_color_manual(values = my_colours) +
  scale_fill_manual(values = my_colours) +
  scale_alpha_manual(values = my_alpha) +
  theme(aspect.ratio = 1) +
  guides(color = guide_legend(nrow = 3)) +
  coord_flip()
ggsave(file = "3BD/genTime.svg", plot = fig)
