### required packages
library(tidyverse)
library(ggplot2)
library(rstatix)


kw_data_file <- '522_kw_test_data_V3.csv'


kw_data <- read.csv(kw_data_file)

# perform kruskal-wallis test on x4.pct values by group
k <- kruskal.test(x4.pct~group, data = kw_data)
dunn_test(kw_data, x4.pct~group, p.adjust.method = "bonferroni")


### create plot with ranks as blue dots and align the letters next to them
kw_data %>% 
  ggplot(aes(x = group, y = x4.pct), ) +
  geom_count() +
  theme(axis.text = element_text(size = 18))+
  # add median as horizontal line
  stat_summary(fun.y = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, col = "blue")
  #labs(title = '521')

ggsave('522_kw_test_V2.png')
