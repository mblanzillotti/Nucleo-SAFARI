#===============================================================================

#R Shiny application for automated m/z domain fragment identification
#of proteins and nucleic acids

# MBL 14 March 2024

#Fragment identification Module

#This module handles data input; either a Thermo .raw file or a .csv or .xlsx
#centroid list, including associated metadata.


#Packages=======================================================================

require(data.table)
require(magrittr)


#Functions======================================================================

ID_centroid <- function(theoretical, centroid_data, tolerance = 10){
  
  identification <- centroid_data[(centroid_mz < theoretical*(1+tolerance/1e6)) & (centroid_mz > theoretical*(1-tolerance/1e6))
  ][, ppm_error := (centroid_mz-theoretical)*1e6
  ][order(-abs(ppm_error))
  ][, .(centroid_index)
  ][[1]][1] #Give back a single integer value with the lowest ppm error, not all identified centroids
  
  if(identical(identification, integer(0))){
    
    return(NA)
    
  } else {
    
    return(
      unlist(identification)
    )
    
  }
  
}


identify_precursor <- function(precursor_mzs, spectrum, tolerance = 10){
  
  precursor_identifications <- precursor_mzs[, !"mass"
    ][, centroid_index := lapply(theoretical_mz, ID_centroid, spectrum$centroid_data, tolerance = tolerance) %>% unlist()
    ][!is.na(centroid_index)]
  
  precursor_identifications <- merge(
    precursor_identifications, spectrum$centroid_data, by = "centroid_index"
  )
  
  
  return(precursor_identifications)
  
}


identify_fragments <- function(fragment_mzs, spectrum, mass_tolerance = 10){
  
  fragment_identifications <- fragment_mzs[, !"mass"
  ][, centroid_index := lapply(theoretical_mz, ID_centroid, spectrum$centroid_data, tolerance = mass_tolerance) %>% unlist()
  ][!is.na(centroid_index)
  ][, sequence := c(0, cumsum(diff(nuclide)>1)), by = .(ion_type, position, charge)
  ][, sequence_prob := sum(prob), by = .(ion_type, charge, position, sequence)
  ][, .SD[max(sequence_prob) == sequence_prob], by = .(ion_type, charge, position)
  ][sequence_prob > 0.6
  ][, .(ion_type, position, charge, nuclide, prob, theoretical_mz, centroid_index)]
  
  hydrogen_shifts <-data.table(
    ion_type = c("a", "a-B", "b", "b-B", "c", "c-B", "d", "d-B", "w", "w-B", "x", "x-B", "y", "y-B", "z", "z-B"),
    shift = list(c(-1, -2), c(-1, -2, -3), -1, -1, -1, -1, c(-1, -2), -1, -1, -1, -1, -1, -1, -1, -1, -1)
  )
  
  
  #might as well search for all the fragments? Given the minimal additional processing time with the new framework
  hydrogen_shift_identifications <- fragment_mzs[unique(fragment_identifications[, .(ion_type, position, charge)]), on = .(ion_type, position, charge)
  ][, .(ion_type, position, charge, nuclide, prob, theoretical_mz)
  ][hydrogen_shifts, on = .(ion_type)
  ][, .(shift = unlist(shift)), by = .(ion_type, position, charge, nuclide, prob, theoretical_mz)
  ][, theoretical_mz := theoretical_mz + 1.00783*shift/abs(charge)
  ][, ion_type := paste0(ion_type, shift)
  ][, !c("shift")
  ][, centroid_index := lapply(theoretical_mz, ID_centroid, spectrum$centroid_data, tolerance = mass_tolerance) %>% unlist()
  ][, .SD[any(nuclide == 0)], by = .(ion_type, position, charge)#Is this what I need to only keep ones where nuclide 0 is found?
  ][!is.na(centroid_index)
  ][, sequence := c(0, cumsum(diff(nuclide)>1)), by = .(ion_type, position, charge)
  ][, sequence_prob := sum(prob), by = .(ion_type, charge, position, sequence)
  ][, .SD[max(sequence_prob) == sequence_prob], by = .(ion_type, charge, position)
  ][sequence_prob > 0.6
  ][, .(ion_type, position, charge, nuclide, prob, theoretical_mz, centroid_index)]
  
  
  fragment_identifications <- rbind(
    fragment_identifications, hydrogen_shift_identifications
  )[, parent_type := str_extract(ion_type, "[abcdwxyz]((\\-B)|(\\-H2O))?")] 
  
  
  return(fragment_identifications)
  
}


fit_predicted_intensity <- function(fragment_data){
  
  
  fragment_data_dcast <- dcast(fragment_data,
                               centroid_intensity ~ ion_type, value.var = "prob", fill = 0
  )
  
  solved_intensity <- fragment_data_dcast[, !"centroid_intensity"] %>% as.matrix() %>%
    nnls::nnls(b = fragment_data_dcast$centroid_intensity) %>% 
    .$x
  
  
  return(
    fragment_data[data.table(ion_type = names(fragment_data_dcast)[-1], contributed_intensity = solved_intensity), on = .(ion_type)
    ][, .(contributed_intensity)
    ][[1]]
  )
  
}


validate_fragment_identifications <- function(fragment_identifications, intensity_tolerance = 25, spectrum){
  
  validated_fragments <- spectrum$centroid_data[fragment_identifications, on = .(centroid_index)
  ][, contributed_intensity := fit_predicted_intensity(.SD), by = .(parent_type, position, charge), .SDcols = c("ion_type", "prob", "centroid_index", "centroid_intensity")
  ][contributed_intensity > 0
  ][, predicted_intensity := prob*contributed_intensity
  ][, intensity_error := (centroid_intensity - sum(predicted_intensity))/sum(predicted_intensity), by = .(parent_type, position, charge, centroid_index)
  ][, .SD[weighted.mean(abs(intensity_error), prob*contributed_intensity) <= intensity_tolerance/100], by = .(parent_type, position, charge)
  ][, ppm_error := (theoretical_mz-centroid_mz)/theoretical_mz*1e6]
  
  
  return(validated_fragments)
  
}


#Module Server==================================================================

identification_server <- function(input, output, session, user_inputs, mz, data_readin){
  
  precursor_identifications <- reactive({
    
    req(data_readin$spectrum())
    req(user_inputs$input_sequence())#these shouldn't be necessary right?
    
    identify_precursor(
      mz$precursor_mz_data(),
      data_readin$spectrum(),
      user_inputs$ppm_tolerance()
    )
    
  })
  
  fragment_identifications <- reactive({
    
    req(data_readin$spectrum())
    req(user_inputs$input_sequence())
    
    identify_fragments(
      mz$fragment_mz_data(),
      data_readin$spectrum(),
      mass_tolerance = user_inputs$ppm_tolerance()#,
      #biomolecule = user_inputs$biomolecule()
    ) 
    
  })
  
  validated_identifications <- reactive({
    
    req(data_readin$spectrum())
    req(user_inputs$input_sequence())
    
    validate_fragment_identifications(
      fragment_identifications(),
      user_inputs$intensity_tolerance(),
      data_readin$spectrum()
    )
    
  })
  
  output$precursor_identifications_table <- renderTable(
    
    precursor_identifications()
    
  )
  
  output$fragment_identifications_table <- renderTable(
    
    validated_identifications()
    
  )
  
  
  return(
    list(
      precursor_identifications = reactive({ precursor_identifications() }),
      fragment_identifications = reactive({ validated_identifications() })
    )
  )
  
}


#Module UI======================================================================

identification_UI <- function(id){
  
  ns <- NS(id)
  
  
  # return(
  
  # )
  
}