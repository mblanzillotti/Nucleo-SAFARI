#===============================================================================

#R Shiny application for automated m/z domain fragment identification
#of proteins and nucleic acids

# MBL 14 March 2024

#User Inputs Module

#This module handles reactive user inputs and ferries them to the appropriate
#modules downstream for parsing, composition determination, fragment matching,
#and visualization.


#UI=============================================================================

user_inputs_UI <- function(id){
  
  ns <- NS(id)
  
 
  sidebar_layout <- tagList(
    textInput(ns("input_sequence"),
      "Sequence Input"
    ),
    # selectInput(ns("biomolecule"),
    #   "Biomolecule",
    #   choices = c(
    #     "Nucleic Acid" = "nucleic_acid",
    #     "Protein" = "protein"
    #   ),
    #   selected = "nucleic_acid"
    # ),
    tags$hr(),
    checkboxInput(ns("neutral_loss_toggle"),
      "Neutral Losses",
      value = T
    ),
    checkboxInput(ns("deconvoluted_toggle"),
      "Deconvoluted Data",
      value = F
    ),
    tags$hr(),
    numericInput(ns("precursor_charge"),
      "Precursor Charge State",
      -1L,
      min = -50L,
      max = 50L,
      step = 1L
    ),
    numericInput(ns("n_charge_states"),
      "Number of Fragment Charge States",
      3L,
      min = 1L,
      max = 8L,
      step = 1L
    ),
    tags$hr(),
    fileInput(ns("file_upload"),
      "Upload Rawfile or Mass List"
    ),
    numericInput(ns("ppm_tolerance"),
      "m/z Tolerance (ppm)",
      value = 10L,
      min = 1L,
      max = 50L,
      step = 1L
    ),
    numericInput(ns("intensity_tolerance"),
      "Intensity Tolerance (%)",
      value = 25L,
      min = 5L,
      max = 75L,
      step = 5L
    ),
    checkboxInput(ns("hydrogen_shift_toggle"),
      "Search Hydrogen Shifts",
      value =T
    )
  )
  
   
  return(sidebar_layout)
  
}

#Server=========================================================================

user_inputs_server <- function(input, output, session, which_tab){
  
  
  # observeEvent(input$biomolecule, {
  #   updateNumericInput(
  #     inputId = "precursor_charge", 
  #     value = switch(input$biomolecule,
  #       nucleic_acid = -1,
  #       protein = 1
  #     )
  #   )
  # })
  
  
  return(
    list(
      # Sequence parsing
      input_sequence = reactive({ input$input_sequence }),
      # biomolecule = reactive({ input$biomolecule }),
      # Mass generation
      neutral_loss_toggle = reactive({ input$neutral_loss_toggle }),
      deconvoluted_toggle = reactive({ input$deconvoluted_toggle  }),
      # m/z generation
      precursor_charge = reactive({ input$precursor_charge }),
      n_charge_states = reactive({ input$n_charge_states }),
      # Data read-in
      file_upload = reactive({ input$file_upload }),
      # Search parameters
      hydrogen_shift_toggle = reactive({ input$hydrogen_shift_toggle }),
      ppm_tolerance = reactive({ input$ppm_tolerance }),
      intensity_tolerance = reactive({ input$intensity_tolerance })
    )
  )
  
}