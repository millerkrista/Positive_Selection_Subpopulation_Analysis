#!/usr/bin/env Rscript

#### Script to open TajimasD workbooks and find every midpoint for every sample that 
#### has a D value greater than or equal to 1.5
#### Output is a json file with participant ID: [midpoints...], participant ID: [midpoints....]


# import packages
library(readxl)
library(dplyr)
library(collections)
library(jsonlite)


main <- function() {
  
  # get lists of dataframes from excel workbooks for 01 and 607
  tajimasD_01 <- read_excel_workbook('./data/VRC601_TajimasD_1.xlsx')
  tajimasD_607 <- read_excel_workbook('./data/TajimasD.xlsx')
  tajimasD_mexico <- read_excel_workbook('./data/mexico_TajimasD.xlsx')
  
  # concatenate the lists of tajimasD data
  tajimasD_data_frame <- c(tajimasD_01, tajimasD_607, tajimasD_mexico)
  
  # initialize dictionary to contain IDs and significant midpoints
  sig_midpoints <- dict()
  
  # go through each participant's data
  for (sheet in 1:length(tajimasD_data_frame)) {
    
    # get ID
    name <- names(tajimasD_data_frame[sheet])
 
    # filter data to get midpoints where Tajimas D >= 1.5
    tajimasD_data_frame[[sheet]] <- dplyr::filter(tajimasD_data_frame[[sheet]], D >= 1.5)
    
    # add ID and midpoints to dictionary
    sig_midpoints$set(name, as.numeric(trimws(tajimasD_data_frame[[name]]$Midpoint, whitespace = "[\\h\\v]")))
  }
  
  # convert dictionary object to list
  sig_midpoints_list <- sig_midpoints$as_list()
  
  # convert list to JSON
  sig_midpoints_json <- toJSON(sig_midpoints_list)
  
  # save as JSON file
  write(sig_midpoints_json, file='results/sig_midpoints.json')
}


read_excel_workbook <- function(tajimasD_path) {
  
  # open Tajimas D workbook
  # getting names of all excel sheets in TajimasD workbook
  sheets <- readxl::excel_sheets(tajimasD_path) 
  
  # reading every sheet of TajimasD workbook into list
  df_list <- lapply(sheets, function(x) readxl::read_excel(tajimasD_path, sheet = x)) 
  
  # get each sheet from list as a dataframe
  tajimasD_data_frame <- lapply(df_list, as.data.frame) 
  
  # assigning names to data frames 
  names(tajimasD_data_frame) <- sheets
  
  return(tajimasD_data_frame)
}


main()



