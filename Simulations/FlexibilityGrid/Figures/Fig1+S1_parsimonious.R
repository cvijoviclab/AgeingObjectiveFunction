library(ggplot2)
library(ggpubr)
library(readr)
library(reshape2)
library(tibble)
library(svglite)
library(stringr) 
library(cividis)

# set theme
my_theme <- theme_light(base_size = 20) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title = element_text(size = 20)) + 
  theme(rect = element_rect(fill = "transparent")) +
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14))

# create directories for the plots
dir.create("1")
dir.create("S1")

# go through all result folders and make complete set 
folders <- list.files(path = "../")
completeSet <- tibble()

for (p in folders) {
    
  file <- paste("../", p, "/flexibilityGrid.txt", sep = "")

  if (file.exists(file)) {
    
    # read data
    data <- read_delim(file, "\t", col_names = TRUE, comment = "#")
    
    # use only the data with regulation
    data <- data[data$`regulation factor` == 0.04, ]
    
    # prepare generation times to make it better to plot
    data$`average generation time`[data$`average generation time` > 5] <- 5
    
    # add the data without parsimonious solution 
    # to complete data set including all objectives
    if (grepl("_p", p)) {
      completeSet <- rbind(completeSet, data)
    }
  }
}

# separate the objectives
completeSet$`objective 1` <- completeSet$objective
completeSet$`objective 2` <- completeSet$objective

for (i in c(1:nrow(completeSet))) {
  splitObjective <- strsplit(completeSet$`objective 1`[i], split = " [+] ")
  completeSet$`objective 1`[i] <- splitObjective[[1]][1]
  completeSet$`objective 2`[i] <- splitObjective[[1]][2]
  
  if (is.na(completeSet$`objective 2`[i])) {
    completeSet$`objective 2`[i] <- splitObjective[[1]][1]
  }
}

# create order forthe objectives for the plots
order1 <- c("growth (MAX)", "NGAM (MAX)", "ATPProduction (MAX)", "NADHProduction (MIN)")
order2 <- c("growth (MAX)", "NGAM (MAX)", "glucoseUptake (MIN)", "NADHProduction (MIN)", 
            "ATPProduction (MAX)", "ATPProduction (MIN)")

# plot replicative lifespan
fig <- ggplot(data = completeSet, aes(x = `flexibility 1`, y = `flexibility 2`,
                                      color = rls, fill = rls)) +
  my_theme +
  geom_tile() +
  facet_grid(factor(`objective 2`, levels = order2) ~
               factor(`objective 1`, levels = order1)) +
  scale_x_continuous(name = expression("flexibility" ~ epsilon[1]),
                     limits = c(0.025, 0.525),
                     breaks = c(0.1, 0.3, 0.5)) +
  scale_y_continuous(name = expression("flexibility" ~ epsilon[2]),
                     limits = c(0.025, 0.225),
                     breaks = c(0, 0.1, 0.2)) +
  theme(axis.text = element_text(size = 14)) +
  scale_fill_gradientn(colours = c("white", "darkgrey", "#02818a"), 
                       breaks = c(0, 15, 30), limits = c(0, 30), 
                       na.value = "white") +
  scale_color_gradientn(colours = c("white", "darkgrey", "#02818a"), 
                        breaks = c(0, 15, 30), limits = c(0, 30), 
                        na.value = "white") +
  theme(aspect.ratio = 0.5)
ggsave("1/rls_parsimonious.svg", plot = fig)

# plot generation time
fig <- ggplot(data = completeSet, aes(x = `flexibility 1`, y = `flexibility 2`,
                                      color = `average generation time`, 
                                      fill = `average generation time`)) +
  my_theme +
  geom_tile() +
  facet_grid(factor(`objective 2`, levels = order2) ~
               factor(`objective 1`, levels = order1)) +
  scale_x_continuous(name = expression("flexibility" ~ epsilon[1]),
                     limits = c(0.025, 0.525),
                     breaks = c(0.1, 0.3, 0.5)) +
  scale_y_continuous(name = expression("flexibility" ~ epsilon[2]),
                     limits = c(0.025, 0.225),
                     breaks = c(0, 0.1, 0.2)) +
  theme(axis.text = element_text(size = 14)) +
  scale_fill_gradientn(colours = c("white", "#02818a", "darkgrey", "darkgrey"), 
                       breaks = c(0, 1.5, 2.3, 5), na.value = "white") +
  scale_color_gradientn(colours = c("white", "#02818a", "darkgrey", "darkgrey"), 
                        breaks = c(0, 1.5, 2.3, 5), na.value = "white") +
  theme(aspect.ratio = 0.5)
ggsave("1/genTime_parsimonious.svg", plot = fig)

# find the wildtypes 
completeSet$wildtype <- (completeSet$rls > 20) & (completeSet$rls < 30) &
                         (completeSet$`average generation time` > 1.5) &
                         (completeSet$`average generation time` < 2.3)

# plot the wildtypes only
fig <- ggplot(data = completeSet, aes(x = `flexibility 1`, 
                                      y = `flexibility 2`,
                                      color = wildtype, fill = wildtype)) +
  my_theme +
  geom_tile() +
  facet_grid(factor(`objective 2`, levels = order2) ~
               factor(`objective 1`, levels = order1)) +
  scale_x_continuous(name = expression("flexibility" ~ epsilon[1]),
                     limits = c(0.025, 0.525),
                     breaks = c(0.1, 0.3, 0.5)) +
  scale_y_continuous(name = expression("flexibility" ~ epsilon[2]),
                     limits = c(0.025, 0.225),
                     breaks = c(0, 0.1, 0.2)) +
  theme(axis.text = element_text(size = 14)) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("white", "#252525")) +
  scale_fill_manual(values = c("white", "#252525")) +
  theme(aspect.ratio = 0.5)
ggsave("S1/wildtype_parsimonious.svg", plot = fig)

# plot wildtpe count
fig <- ggplot(data = completeSet, aes(x = `objective 2`, y = as.integer(wildtype))) +
  my_theme +
  geom_bar(stat = "identity") +
  facet_grid(. ~ factor(`objective 1`, levels = order1)) +
  xlab("second objective") +
  ylab("wildtype parameter count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 10)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) +
  ylim(c(0, 20)) +
  theme(aspect.ratio = 1)
ggsave("1/wildtypeCount_parsimonious.svg", plot = fig)
