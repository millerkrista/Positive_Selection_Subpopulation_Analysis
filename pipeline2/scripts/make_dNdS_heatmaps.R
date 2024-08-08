#!/usr/bin/env Rscript

#####Script to create heatmaps from files containing dN/dS, P(dS<dN) values
#####for each participant ID mapped to their equivalent HXB2 codon position

# load packages
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(glue)


main <- function() {
  
  # get filenames as list
  dNdS_files <- list.files(path='results/tables/', pattern='*_dNdS.csv')
  
  # for every dNdS file,
  for (file in 1:length(dNdS_files)) {
    
    # get name of dNdS file
    dNdS_file <- dNdS_files[file]
    
    # get annotation name
    annotation <- str_replace(dNdS_file, '_dNdS.csv', '')
    
    # get name of corresponding probability file
    prob_file <- str_replace(dNdS_file, 'dNdS.csv', 'probability.csv')
      
    # call function to make heatmap for dNdS and probability
    dNdS_heatmap <- get_heatmaps(glue('results/tables/{dNdS_file}'), 'dN/dS')
    prob_heatmap <- get_heatmaps(glue('results/tables/{prob_file}'), 'P(dS<dN)')
    
    # make pdf containing both heatmaps
    plot <- ggarrange(plotlist = list(dNdS_heatmap, prob_heatmap), ncol = 2, nrow = 1,
              labels = list('A', 'B'))
    annotate_figure(plot, top = text_grob(annotation))
    
    # save heatmaps image
    ggsave(filename = glue('results/plots/{annotation}_heatmap.png'), height = 25, width = 70, limitsize = FALSE)
  }
}


get_heatmaps <- function(filename, val) {
  
  # read in file
  dat <- read.csv(filename, header = TRUE, row.names = 1, check.names = FALSE)
  
  # get data in format for plotting
  dat2 <- dat %>%
    rownames_to_column() %>%
    gather(colname, value, -rowname)
  
  # plot heatmap
  heat <- ggplot(dat2, aes(x = colname, y = rowname)) +
    labs(y = 'Codon Position', x = 'Participant') +
    geom_tile(data = subset(dat2, !is.na(value)), aes(fill = value)) +
    geom_tile(data = subset(dat2,  is.na(value)), aes(colour = 'NA'),
              linetype = 0, fill = "grey") +
    scale_fill_continuous(name = val, type = 'viridis') +
    guides(colour=guide_legend("", override.aes=list(colour="grey")))

  return(heat)
}

main()


##### Use this to get heatmap for entire gene
h1 <- get_heatmaps('results/tables/pos_selection_sites_dNdS.csv', 'dN/dS')
h2 <- get_heatmaps('results/tables/pos_selection_sites_probability.csv', 'P(dS<dN)')

plot <- ggarrange(plotlist = list(h1, h2), ncol = 2, nrow = 1,
                  labels = list('A', 'B'))
annotate_figure(plot, top = text_grob('Env Gene'))

# save heatmaps image
ggsave(filename = 'results/plots/env_gene_heatmap.png', height = 120, width = 30, limitsize = FALSE)


