library(ggplot2)
library(ggpubr)
library(readr)
library(reshape2)
library(tibble)
library(svglite)
library(scales)
library(cividis)
library(stringr) 
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

# create directories for the plots
dir.create("4")

# go through all result folders and make complete set 
paths <- list.files(path = "../")
completeSet <- tibble()

for (p in paths) {
    
  file <- paste("../", p, "/flexibilityGrid.txt", sep = "")

  if (file.exists(file)) {
    
    # read data
    data <- read_delim(file, "\t", col_names = TRUE, comment = "#")
    data <- data[data$`regulation factor` == 0.04, ]
    
    # add the data with parsimonious solution 
    # to complete data set including all objectives
    if (grepl("_p", p)) {
      completeSet <- rbind(completeSet, data)
    }
  }
}

# find the wildtypes 
completeSet$wildtype <- (completeSet$rls > 20) & (completeSet$rls < 30) &
                        (completeSet$`average generation time` > 1.5) &
                        (completeSet$`average generation time` < 2.3)
wildtypes <- subset(completeSet, wildtype == TRUE)

# plot (damage)
fig <- ggplot(data = wildtypes) +
  my_theme +
  geom_boxplot(aes(x = fct_rev(`objective`), y = `damage phase 1`), 
                   fill = "darkgrey", outlier.shape = NA, alpha = 0.5) +
  geom_jitter(aes(x = fct_rev(`objective`), y = `damage phase 1`), 
              color = "darkgrey", size = 1.5) +
  geom_boxplot(aes(x = fct_rev(`objective`), y = `damage phase 2`), 
               fill = "black", outlier.shape = NA, alpha = 0.5) +
  geom_jitter(aes(x = fct_rev(`objective`), y = `damage phase 2`), 
              color = "black", size = 1.5) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) + 
  xlab("objective") +
  ylab(expression("damage at end of phase")) +
  coord_flip()
ggsave("4/damage.svg", plot = fig)

# plot (rls)
fig <- ggplot(data = wildtypes) +
  my_theme +
  geom_boxplot(aes(x = fct_rev(`objective`), y = `rls phase 1`), 
               fill = "darkgrey", outlier.shape = NA, alpha = 0.5) +
  geom_jitter(aes(x = fct_rev(`objective`), y = `rls phase 1`), 
              color = "darkgrey", size = 1.5) +
  geom_boxplot(aes(x = fct_rev(`objective`), y = `rls phase 2`), 
               fill = "black", outlier.shape = NA, alpha = 0.5) +
  geom_jitter(aes(x = fct_rev(`objective`), y = `rls phase 2`), 
              color = "black", size = 1.5) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) + 
  xlab("objective") +
  ylab("number of divisions") +
  coord_flip()
ggsave("4/rls.svg", plot = fig)

# plot (average generation time)
fig <- ggplot(data = wildtypes) +
  my_theme +
  geom_boxplot(aes(x = fct_rev(`objective`), y = `time phase 1`), 
               fill = "darkgrey", outlier.shape = NA, alpha = 0.5) +
  geom_jitter(aes(x = fct_rev(`objective`), y = `time phase 1`), 
              color = "darkgrey", size = 1.5) +
  geom_boxplot(aes(x = fct_rev(`objective`), y = `time phase 2`), 
               fill = "black", outlier.shape = NA, alpha = 0.5) +
  geom_jitter(aes(x = fct_rev(`objective`), y = `time phase 2`), 
              color = "black", size = 1.5) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) + 
  xlab("objective") +
  ylab("time") +
  coord_flip()
ggsave("4/time.svg", plot = fig)