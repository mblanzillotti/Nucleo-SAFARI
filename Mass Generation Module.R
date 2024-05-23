#===============================================================================

#R Shiny application for automated m/z domain fragment identification
#of proteins and nucleic acids

# MBL 14 March 2024

#Mass Generation Module

#This module consists of two connected parts: generation of monoisotopic and
#isotopic masses. Mass generation is dictated by the "deconvoluted_toggle"
#from the user inputs module. Monoisotopic mass calculation is apart of isotopic
#generation, however the converse is untrue and thus toggled off when searching
#deconvoluted data.


#Packages=======================================================================

require(data.table)
require(IsoSpecR)
require(magrittr)


#Functions======================================================================

#Mass calculation helpers
calculate_monoisotopic <- function(element, amount){
  
  data.table(
    element = element,
    amount = amount
  )[element_masses, on = .(element)
  ][amount != 0L
  ][, .(mass = amount*mono_mass)
  ][, .(sum(mass))
  ][[1]]
  
}

IsoSpec_dist <- function(element, amount){
  
  monoisotopic_mass <- calculate_monoisotopic(element, amount)
  
  composition <- data.table(
    element = element,
    amount = amount
  )[, .(element = list(element), amount = list(amount))
  ][, .(IsoSpec = purrr::pmap(list(amount, element), `names<-`))
  ][, IsoSpec
  ]
  
  nominal_dist <- IsoSpecR::IsoSpecify(composition[[1]], 0.97) %>% 
    as.data.table()
  
  nominal_dist <- nominal_dist[, nuclide := mass-monoisotopic_mass
  ][, nuclide := round(nuclide) %>% as.integer()
  ][, .(prob, nuclide)
  ][, .(prob = sum(prob)), by = .(nuclide)
  ][order(nuclide)
  ][, mass := nuclide*1.002185 + monoisotopic_mass]
  
  return(nominal_dist)
  
}


#Precursor mass calculations
calculate_precursor_monoisotopic <- function(precursor_composition){
  
  precursor_monoisotopic <- precursor_composition[, .(monoisotopic_mass = calculate_monoisotopic(element, amount))
  ][[1]]
  
  
  return(precursor_monoisotopic)
  
}

calculate_precursor_isotopes <- function(precursor_composition){
  
  precursor_isotopes <- precursor_composition[, .(element, amount)
  ][, .(amount = sum(amount)), by = .(element)
  ][amount != 0L
  ][, .(isotopes = list(IsoSpec_dist(element, amount)))
  ][, rbindlist(isotopes)]
  
  
  return(precursor_isotopes)

}


#Fragment mass calculations
calculate_fragment_monoisotopics <- function(fragment_compositions){
  
  fragment_monoisotopics <- fragment_compositions[, .(monoisotopic_mass = calculate_monoisotopic(element, amount)), by = .(position, ion_type)]
  
  
  return(fragment_monoisotopics)
  
}

calculate_fragment_isotopes <- function(fragment_compositions){
  
  fragment_isotopes <- fragment_compositions[, monoisotopic_mass := calculate_monoisotopic(element, amount), by = .(position, ion_type)
  ][, .(isotopes = list(IsoSpec_dist(element, amount))), by = .(position, ion_type)
  ][, rbindlist(isotopes), by = .(ion_type, position)]
  
  
  return(fragment_isotopes)
  
}

#Module Server==================================================================

mass_server <- function(input, output, session, user_inputs, parsing, compositions){
  
  #Precursor Data
  precursor_mass_data <- reactive({
    
    req(user_inputs$input_sequence())
    
    if(user_inputs$deconvoluted_toggle()){
      calculate_precursor_monoisotopic(compositions$precursor_composition())
    } else {
      calculate_precursor_isotopes(compositions$precursor_composition())
    }
    
  })
  
  output$precursor_mass_summary <- renderText(
    
    if(user_inputs$deconvoluted_toggle()){
      paste0(
        "Precursor Monoisotopic Mass: ", 
        round(precursor_mass_data(), 3),
        " Da"
      )
    } else {
      paste0(
        "Precursor Average Mass: ", 
        round(precursor_mass_data()[, .(average_mass = weighted.mean(mass, prob))][[1]], 3),
        " Da"
      )
    }
    
  )
  
  #Fragment Data
  fragment_mass_data <- reactive({
    
    req(user_inputs$input_sequence())
    
    if(user_inputs$deconvoluted_toggle()){
      calculate_fragment_monoisotopics(compositions$fragment_compositions())
    } else {
      calculate_fragment_isotopes(compositions$fragment_compositions())
    }
    
  })
  
  output$fragment_mass_table <- renderTable(
    
    if(user_inputs$deconvoluted_toggle()){
      fragment_mass_data()
    } else {
      fragment_mass_data()[, .(average_mass = weighted.mean(mass, prob)), by = .(ion_type, position)]
    },
    
    digits = 3,
    
  )
  
  output$fragment_mass_download <- downloadHandler(
    filename = function(){
      paste0("Theoretical Fragment ",
        ifelse(user_inputs$deconvoluted_toggle(), "Monoisotopic ", "Isotopic "),
        "Mass Data.csv"
      )
    },
    content = function(file){
      write.csv(fragment_mass_data(), file)
    }
  )
  
  
  return(
    list(
      precursor_mass_data = reactive({ precursor_mass_data() }),
      fragment_mass_data = reactive({ fragment_mass_data() })
    )
  )
  
}


#Module UI======================================================================

mass_UI <- function(id){
  
  ns <- NS(id)
  
  
  return(
    tagList(
      textOutput(ns("precursor_mass_summary")),
      downloadButton(ns("fragment_mass_download"), "Download Theoretical Fragment Masses"),
      tableOutput(ns("fragment_mass_table"))
    )
  )
  
}