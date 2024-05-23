#===============================================================================

#R Shiny application for automated m/z domain fragment identification
#of proteins and nucleic acids

# MBL 14 March 2024

#Data read-in Module

#This module handles data input; either a Thermo .raw file or a .csv or .xlsx
#centroid list, including associated metadata.


#Packages=======================================================================

require(data.table)
require(magrittr)
require(stringr)
require(rawrr)
require(tools)


#Functions======================================================================

read_rawfile <- function(rawfile_path){
  
  spectrum <- readSpectrum(rawfile_path, 1)
  
  isolation_width <- spectrum[[1]]$`MS2 Isolation Width` %>% as.numeric()
  isolation_center <- spectrum[[1]]$scanType %>% 
    str_extract("ms2 \\d{3,4}\\.\\d+\\@") %>% 
    str_extract("\\d{3,4}\\.\\d+") %>% 
    as.numeric()
  
  if(str_detect(spectrum[[1]]$scanType, "\\@")){#Check if not an MS1
    
    simple_spectrum <- list(
      
      centroid_data = data.table(
        centroid_mz = spectrum[[1]]$centroid.mZ,
        centroid_intensity = spectrum[[1]]$centroid.intensity
      )[centroid_intensity > 10
      ][, centroid_index := 1:length(centroid_mz)],
      
      profile_data = data.table(
        mz = spectrum[[1]]$mZ,
        intensity = spectrum[[1]]$intensity
      ),
      
      polarity = str_extract(spectrum[[1]]$scanType, "MS [+-]") %>% 
        str_remove("MS "),
      
      mz_range = str_extract(spectrum[[1]]$scanType, "\\d+\\.\\d+\\-\\d+\\.\\d+") %>% str_split("-", simplify = T) %>% as.numeric(),
      
      activation = str_extract(spectrum[[1]]$scanType, "\\@(hcd)|(uvpd)|(cid)") %>% 
        str_remove("\\@"),
      
      precursor = data.table(
        precursor_mz = spectrum[[1]]$centroid.mZ,
        precursor_intensity = spectrum[[1]]$centroid.intensity
      )[(precursor_mz > isolation_center - isolation_width/2) & (precursor_mz < isolation_center + isolation_width/2)]
      
    )
    
  } else {
    
    simple_spectrum <- list(
      centroid_data = data.table(
        centroid_mz = spectrum[[1]]$centroid.mZ,
        centroid_intensity = spectrum[[1]]$centroid.intensity
      )[centroid_intensity > 10
      ][, centroid_index := 1:length(centroid_mz)],
      
      profile_data = data.table(
        mz = spectrum[[1]]$mZ,
        intensity = spectrum[[1]]$intensity
      ),
      
      polarity = str_extract(spectrum[[1]]$scanType, "MS [+-]") %>% 
        str_remove("MS "),
      
      mz_range = str_extract(spectrum[[1]]$scanType, "\\d+\\.\\d+\\-\\d+\\.\\d+") %>% str_split("-", simplify = T) %>% as.numeric()
      
    )
    
  }
  
  
  return(simple_spectrum)
  
}
#Reads in centroids and profile data from a rawfile stored in a list of data tables with some metadata


#Module Server==================================================================

data_readin_server <- function(input, output, session, user_inputs){
  
  spectrum <- reactive({
    
    req(user_inputs$file_upload())

    ext <- tools::file_ext(user_inputs$file_upload()$name)

    switch(ext,
      raw = read_rawfile(user_inputs$file_upload()$datapath),
      RAW = read_rawfile(user_inputs$file_upload()$datapath),
      validate("Invalid file; Please upload a .raw file")
    )
    
    #read_rawfile(input$file_upload$datapath)
    
  })
  
  output$spectrum_centroid_table <- renderTable(
    spectrum()$centroid_data
  )
  
  return(
    list(
      spectrum = reactive({ spectrum() })
    )
    # list(
    #   centroid_data = switch(ext(),
    #     csv = ,
    #     tsv = ,
    #     raw = ,
    #     RAW = 
    #   ),
    #   profile_data = switch(ext(),
    #     csv = ,
    #     tsv = ,
    #     raw = ,
    #     RAW =              
    #   )
    # )
  )
  
}


#Module UI======================================================================

data_readin_UI <- function(id){
  
  ns <- NS(id)
  

  return(
    tagList(
      tableOutput(ns("spectrum_centroid_table"))
    )
    # list(
    #   sidebar = tagList(
    #     fileInput(ns("file_upload"),
    #       "Rawfile or Mass List"
    #     )
    #   ),
    #   main_panel = tagList(
    #     tableOutput(ns("spectrum_centroid_table"))
    #   )
    # )
  )
  
}