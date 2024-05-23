#===============================================================================

#R Shiny application for automated m/z domain fragment identification
#of proteins and nucleic acids

# MBL 14 March 2024

#m/z Generation Module

#This module consists of two connected parts: generation of monoisotopic and
#isotopic masses. Mass generation is dictated by the "deconvoluted_toggle"
#from the user inputs module. Monoisotopic mass calculation is apart of isotopic
#generation, however the converse is untrue and thus toggled off when searching
#deconvoluted data.

#Mass generation is agnostic to the type of biomolecule selected (protein or
#nucleic acid), and is thus applied to the precursor and fragment compositions
#indiscriminately by position and ion type (where applicable)


#Packages=======================================================================

require(data.table)
require(IsoSpecR)
require(magrittr)


#Functions======================================================================

calculate_precursor_mzs <- function(precursor_isotopes, precursor_charge, adduct = 1.00783, mz_range = c(200, 2000)){
  
  precursor_mzs <- precursor_isotopes[, theoretical_mz := (mass + precursor_charge*adduct)/abs(precursor_charge)]
  
  
  return(precursor_mzs)
  
}

calculate_fragment_mzs <- function(fragment_isotopes, precursor_charge, n_charge_states = 3L, adduct = 1.00783, mz_range = c(200, 2000)){
  
  fragment_mzs <- fragment_isotopes[, .(position, ion_type)] %>% unique()
  
  fragment_mzs <- fragment_mzs[, adj_position := position
  ][str_detect(ion_type, "[wxyz]"), adj_position := max(position)-position
  ][, ratio := precursor_charge/max(adj_position)
  ][, charge := round(adj_position*ratio)
  ][, lapply(.SD, `+`, 1L:n_charge_states - ceiling(n_charge_states/2)), by = .(ion_type, position), .SDcols = c("charge")#wow ok that's easy
  ][sign(charge) == sign(precursor_charge)
  ][abs(charge) <= abs(precursor_charge)
  ][, .(ion_type, position, charge)
  ][fragment_isotopes, on = .(ion_type, position), allow.cartesian = T
  #][, rbindlist(isotopes), by = .(ion_type, position, charge)
  ][, theoretical_mz := (mass + adduct*charge)/abs(charge)
  ][, .SD[all(theoretical_mz %between% mz_range)], by = .(ion_type, position, charge), .SDcols = c("nuclide", "prob", "mass", "theoretical_mz")]
  
  
  return(fragment_mzs)
  
}


#Module Server==================================================================

mz_server <- function(input, output, session, user_inputs, mass){
  
  precursor_mz_data <- reactive({
    
    req(user_inputs$deconvoluted_toggle() == F)
    
    calculate_precursor_mzs(
      mass$precursor_mass_data(),
      user_inputs$precursor_charge()
    )
    
  })
  
  output$precursor_average_mz <- renderText(
    paste0(
      "Precursor Average m/z: ", 
      round(precursor_mz_data()[, .(average_mz = weighted.mean(theoretical_mz, prob))][[1]], 3)
    )
  )
  
  output$precursor_mz_table <- renderTable(
    precursor_mz_data()
  )
  
  
  fragment_mz_data <- reactive({
    
    req(user_inputs$deconvoluted_toggle() == F)
    
    calculate_fragment_mzs(
      mass$fragment_mass_data(),
      user_inputs$precursor_charge(),
      user_inputs$n_charge_states()#,
      #mass_range = data_readin$spectrum$mass_range()
    )#add an if statement here to worry about mass range if it's a mass list?
    
  })
  
  output$fragment_mz_download <- downloadHandler(
    filename = function(){
      "Theoretical Fragment mz Data.csv"
    },
    content = function(file){
      write.csv(fragment_mz_data(), file)
    }
  )
  
  
  return(
    list(
      precursor_mz_data = reactive({ precursor_mz_data() }),
      fragment_mz_data = reactive({ fragment_mz_data() })
    )
  )
  
}


#Module UI======================================================================

mz_UI <- function(id){
  
  ns <- NS(id)
  
  
  return(
    tagList(
      textOutput(ns("precursor_average_mz")),
      downloadButton(ns("fragment_mz_download"),
        "Download Theoretical Fragment m/zs"
      ),
      tableOutput(ns("precursor_mz_table"))
    )
  )
  
}