#!/usr/bin/env Rscript

##### Script to get HIV env gene annotations mapped to codon position for each patient consensus sequence

# imports
library(readxl)
library(tidyverse)
library(glue)


# get command line arguments
args <- commandArgs(trailingOnly = TRUE)
xl_dir <- args[[2]]
results_dir <- args[[3]]

# read in excel workbook of annotations for HIV envelope gene
xl_annotation <- read_xlsx(glue('{xl_dir}/Env_features[8422].xlsx'), col_names = TRUE)

# get new dataframe containing just the annotations and HXB2 codon position
HXB2_annotations <- xl_annotation[10:866, c('Regions_1: signal peptide, cleavage sites, disulfied bonds, hypervariable regions,  integrin α4β7 binding, Lectin DC-sign binding, gp41 regions: Leucine/Isoleucine zipper',
                                            'Regions_2: Variable loops, gp41 regions:  Kennedy epitope, GxxxG motif, Cytoplasmic tail LLP1, LLP2, LLP3; RxKR motif',
                                            'Regions_3: gp120 Regions related to the CD4 binding site, gp41 regions: fusion peptide, Leucine/Isoleucine zipper, fusion domain, YXXL motif',
                                            'HXB2 position')]

# rename columns of dataframe for easier access
HXB2_annotations <- HXB2_annotations %>%
  rename(
    Regions_1 = `Regions_1: signal peptide, cleavage sites, disulfied bonds, hypervariable regions,  integrin α4β7 binding, Lectin DC-sign binding, gp41 regions: Leucine/Isoleucine zipper`,
    Regions_2 = `Regions_2: Variable loops, gp41 regions:  Kennedy epitope, GxxxG motif, Cytoplasmic tail LLP1, LLP2, LLP3; RxKR motif`,
    Regions_3 = `Regions_3: gp120 Regions related to the CD4 binding site, gp41 regions: fusion peptide, Leucine/Isoleucine zipper, fusion domain, YXXL motif`
  )

# get all position mapping files as a list
mapped_pos_files <- list.files(path = results_dir, 
                               pattern = '*_position_mappings.csv')

# for every position mapping file,
for (i in 1:length(mapped_pos_files)) {
  
  # get the individual ID from the file name
  individual_id <- str_replace(mapped_pos_files[i], '_position_mappings.csv', '')
  
  # get file as dataframe
  pos_map_file <- read.csv(glue('{results_dir}/{mapped_pos_files[i]}'))
  
  # create a new dataframe to store mapped codon positions (and, later, annotations)
  mapping_df <- pos_map_file[, c('HXB2.Position', 
                                 'Consensus.Aligned.to.HXB2.Position', 
                                 'Original.Consensus.Position')]
  
  # get new columns to store annotations
  mapping_df$Regions1 <- NA
  mapping_df$Regions2 <- NA
  mapping_df$Regions3 <- NA
  
  
  # for each HXB2 codon position in annotation dataframe, 
  for (p in 1:length(HXB2_annotations$`HXB2 position`)) {
    
    # check if Regions1 has an associated annotation
    if (is.na(HXB2_annotations$Regions_1[p]) == FALSE) {
      # if yes, add that annotation to the index of mapping dataframe where that HXB2 codon is
      position <- which(mapping_df$HXB2.Position == HXB2_annotations$`HXB2 position`[p])
      mapping_df$Regions1[position] <- HXB2_annotations$Regions_1[p]
    }
    
    # check if Regions2 has an associated annotation
    if (is.na(HXB2_annotations$Regions_2[p]) == FALSE) {
      # if yes, add that annotation to the index of mapping dataframe where that HXB2 codon is
      position <- which(mapping_df$HXB2.Position == HXB2_annotations$`HXB2 position`[p])
      mapping_df$Regions2[position] <- HXB2_annotations$Regions_2[p]
    }
    
    # check if Regions3 has an associated annotation
    if (is.na(HXB2_annotations$Regions_3[p]) == FALSE) {
      # if yes, add that annotation to the index of mapping dataframe where that HXB2 codon is
      position <- which(mapping_df$HXB2.Position == HXB2_annotations$`HXB2 position`[p])
      mapping_df$Regions3[position] <- HXB2_annotations$Regions_3[p]
    }
    
    # if no annotation available, next
    else {
      next
    }
  }
  # save dataframe as csv file
  write.csv(mapping_df, file = glue('{results_dir}{individual_id}_consensus_annotations.csv'))
}


