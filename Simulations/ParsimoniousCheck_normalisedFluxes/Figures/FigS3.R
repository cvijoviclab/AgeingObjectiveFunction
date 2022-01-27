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
dir.create("S3")

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
diff$relDiff <- diff$absDiff / dataFalse$`average flux`

diff$bigChange <- diff$relDiff >= 1 - 0.01 | diff$relDiff <= -1 + 0.01

# get counts of decreasing, increasing and unchanged fluxes
countData_I <- subset(diff, phase == "phase I") %>% 
  group_by(pathway, bigChange) %>%  
  summarise(count = n()) %>% 
  mutate(percentage = count/sum(count))

# plot counts for fluxes that change more than 50% in phase I
fig <- ggplot(data = subset(countData_I, bigChange == 1), 
                            aes(x = pathway, y = percentage * 100)) +
  my_theme + 
  geom_bar(stat = "identity", color = "black", alpha = 0.8, width = 0.7) +
  xlab("pathway") +
  ylab("percentage of fluxes") +
  ylim(c(0, 40)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 12)) +
  scale_x_discrete(breaks = levels(as.factor(data$pathway)),
                   labels = function(x) str_wrap(x, width = 30)) +
  theme(aspect.ratio = 0.4)
ggsave(file = "S3/count_I.svg", plot = fig)


# get counts of decreasing, increasing and unchanged fluxes
countData_II <- subset(diff, phase == "phase II") %>% 
  group_by(pathway, bigChange) %>% 
  summarise(count = n()) %>% 
  mutate(percentage = count/sum(count))

# plot counts for fluxes that change more than 50% in phase I
fig <- ggplot(data = subset(countData_II, bigChange == 1), 
              aes(x = pathway, y = percentage * 100)) +
  my_theme + 
  geom_bar(stat = "identity", color = "black", alpha = 0.8, width = 0.7) +
  xlab("pathway") +
  ylab("percentage of fluxes") +
  ylim(c(0, 40)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 12)) +
  scale_x_discrete(breaks = levels(as.factor(data$pathway)),
                   labels = function(x) str_wrap(x, width = 30)) +
  theme(aspect.ratio = 0.4)
ggsave(file = "S3/count_II.svg", plot = fig)

