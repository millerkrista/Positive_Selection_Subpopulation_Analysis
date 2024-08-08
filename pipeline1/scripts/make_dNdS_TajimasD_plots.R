#!/usr/bin/env Rscript

####### Script to create plots of dN/dS on y-axis1, TajimasD on y-axis2, 
### and codon position on x-axis with gene annotation overylay

# load packages
library(stringr)
library(plotrix)
library(glue)
library(readxl)
library(wesanderson)


main <- function() {
  
  # get lists of dataframes from excel workbooks for 01 and 607
  tajimasD_01 <- read_excel_workbook('data/VRC601_TajimasD_1.xlsx')
  tajimasD_607 <- read_excel_workbook('data/TajimasD.xlsx')
  tajimasD_mexico <- read_excel_workbook('data/mexico_TajimasD.xlsx')
  
  # concatenate the lists of tajimasD data
  tajimasD_data_frame <- c(tajimasD_01, tajimasD_607, tajimasD_mexico)
  
  # load in fubar file names as list
  fubar_files <- list.files(path = 'data/', pattern = '*_simple.txt')
  
  # for every fubar data file
  for (i in 1:length(fubar_files)) {
    
    # call function to:
    dataframes <- process_dataframes(glue('data/{fubar_files[i]}'), tajimasD_data_frame)
    
    # 1. get fubar data as dataframe
    fubar <- dataframes[[1]]
    
    # 2. get participant ID
    individual_id <- dataframes[[2]]
    
    # 3. get tajimasD dataframe in codon position
    tajimasD <- dataframes[[3]]
    
    # read in file of annotations
    consensus_annotation <- read.csv(glue('data/{individual_id}_consensus_annotations.csv'))
    
    # make and save plot
    get_plot(individual_id, fubar, tajimasD, consensus_annotation)
    
  }
}


read_excel_workbook <- function(tajimasD_path) {
  
  # get names of all excel sheets in TajimasD workbook
  sheets <- readxl::excel_sheets(tajimasD_path) 
  
  # read every sheet of TajimasD workbook into list
  df_list <- lapply(sheets, function(x) readxl::read_excel(tajimasD_path, sheet = x)) 
  
  # get each sheet from list as a dataframe
  tajimasD_data_frame <- lapply(df_list, as.data.frame) 
  
  # assign names to data frames 
  names(tajimasD_data_frame) <- sheets
  
  return(tajimasD_data_frame)
}


process_dataframes <- function(fubar_file, tajimasD_data_frame) {
  
  # get file as a dataframe
  fubar <- read.table(fubar_file, header = TRUE)
  
  # get individual ID from fubar name
  individual_id <- str_extract(fubar_file, "(?<=data/)(.*?)(?=_rev2miss\\.aln\\.FUBAR_simple\\.txt$)|(?<=data/)(.*?)(?=\\.aln\\.FUBAR_simple\\.txt$)")
  print(individual_id)
  
  # get TajimasD excel sheet as a dataframe
  d_file <- tajimasD_data_frame[[individual_id]]
  
  # make new dataframe from TajimasD containing midpoint and D
  tajimasD <- d_file[, c("Midpoint", "D")]
  
  # make sure to trim white space, convert to numeric type
  tajimasD$Midpoint <- as.numeric(str_trim(tajimasD$Midpoint))
  
  # get codon site by taking Midpoint divided by 3
  tajimasD$CodonSite <- tajimasD$Midpoint / 3
  
  return(list(fubar, individual_id, tajimasD))
}


get_plot <- function(individual_id, fubar, tajimasD, consensus_annotation) {
  
  # initialize some variables for plotting 
  dNdS_upper_limit <- max(fubar$dN.dS) + 20
  TajimasD_upper_limit <- max(tajimasD$D) + 0.5
  
  # Group the data by unique annotations and get start and stop positions
  Regions1_df <- aggregate(Consensus.Aligned.to.HXB2.Position ~ Regions1, data = consensus_annotation, FUN = function(x) c(min(x), max(x)))
  Regions2_df <- aggregate(Consensus.Aligned.to.HXB2.Position ~ Regions2, data = consensus_annotation, FUN = function(x) c(min(x), max(x)))
  Regions3_df <- aggregate(Consensus.Aligned.to.HXB2.Position ~ Regions3, data = consensus_annotation, FUN = function(x) c(min(x), max(x)))
  
  # get color palettes for each annotation type
  Region1_palette <- wes_palette(length(Regions1_df$Regions1), name = "AsteroidCity1", type = "continuous")
  Region2_palette <- rainbow(length(Regions2_df$Regions2))
  Region3_palette <- wes_palette(length(Regions3_df$Regions3), name = "Darjeeling1", type = "continuous")
  
  # open file to store plot
  png(glue('results/plots/{individual_id}.png'), width = 1500, height = 700)
  
  # create layout to store the two plots; 4 rows with second plot being 3x bigger than first
  layout(matrix(c(1,2,2,2), nrow = 4, ncol = 1, byrow = TRUE))
  
  # plot rectangles in the first row
  plot.new()
  plot.window(xlim = range(fubar$site), ylim = c(0, 5))
  
  # Regions 1 annotations
  for (i in 1:nrow(Regions1_df)) {
    # get codon positions for start and stop based on the aligned consensus
    start_codon <- Regions1_df$Consensus.Aligned.to.HXB2.Position[i, 1]
    stop_codon <- Regions1_df$Consensus.Aligned.to.HXB2.Position[i, 2]
    
    # get codon positions for start and stop based on the original consensus
    start_pos <- which(consensus_annotation$Original.Consensus.Position == start_codon)
    stop_pos <- which(consensus_annotation$Original.Consensus.Position == stop_codon)
    
    # Draw the rectangle
    rect(start_pos, 5, stop_pos, 4, col = Region1_palette[i], border = NA)
  }
  
  # Regions 2 annotations
  for (i in 1:nrow(Regions2_df)) {
    # get codon positions for start and stop based on the aligned consensus
    start_codon <- Regions2_df$Consensus.Aligned.to.HXB2.Position[i, 1]
    stop_codon <- Regions2_df$Consensus.Aligned.to.HXB2.Position[i, 2]
    
    # get codon positions for start and stop based on the original consensus
    start_pos <- which(consensus_annotation$Original.Consensus.Position == start_codon)
    stop_pos <- which(consensus_annotation$Original.Consensus.Position == stop_codon)
    
    # Draw the rectangle
    rect(start_pos, 3, stop_pos, 2, col = Region2_palette[i], border = NA)
  }
  
  # Regions 3 annotations
  for (i in 1:nrow(Regions3_df)) {
    # get codon positions for start and stop based on the aligned consensus
    start_codon <- Regions3_df$Consensus.Aligned.to.HXB2.Position[i, 1]
    stop_codon <- Regions3_df$Consensus.Aligned.to.HXB2.Position[i, 2]
    
    # get codon positions for start and stop based on the original consensus
    start_pos <- which(consensus_annotation$Original.Consensus.Position == start_codon)
    stop_pos <- which(consensus_annotation$Original.Consensus.Position == stop_codon)
    
    # Draw the rectangle
    rect(start_pos, 1, stop_pos, 0, col = Region3_palette[i], border = NA)
  }
  
  # add legend to the first plot
  legend('topright', legend=c(Regions1_df$Regions1, Regions2_df$Regions2, Regions3_df$Regions3),
         fill = c(Region1_palette, Region2_palette, Region3_palette),
         bty = 'n', cex = .8, inset = c(-.5, 0), xpd = TRUE)
  
  # plot Tajimas D and dN/dS vs. codon site in the second row
  twoord.plot(fubar$site, fubar$dN.dS, tajimasD$CodonSite, tajimasD$D, 
              main = glue("Individual {individual_id}"), 
              xlab = "Codon Site", ylab = "dN/dS", rylab = "Tajima's D", 
              lcol = 1, rcol = 4,
              type=c("p","l"),
              lytickpos = seq(0, dNdS_upper_limit, by = 20),
              lylim = c(0, dNdS_upper_limit))
              # rylim = (c(-(TajimasD_upper_limit), TajimasD_upper_limit)))
              # ylab.at = (dNdS_upper_limit / 2), 
              # rylab.at = 0)

  
  # finish plot
  dev.off()
}


main()
