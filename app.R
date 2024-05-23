#Fragment Identification App

#===============================================================================

#R Shiny application for automated m/z domain fragment identification
#of proteins and nucleic acids

# MBL 11 March 2024

#Application Body

#This file contains the server and UI combining several modules to make the app run.
#current version as of 11 March 2024 is in testing, and only parses sequences

#MBL 14 March 2024

#Reorganized UI elements out of user_inputs so inputs live within the modules
#where they are first needed.
  #Note - this does cause "sequence" do disappear when changing tabs. The data
  #is retained, (I can see mass properly) but it ain't it
  #Might need to revert back to an inputs module, which won't be too bad, but hey worth a shot.

  #Will input be carried across the modules? do I need to specify that argument in callModule?

#Added protein functionality, which does require proofreading.

#MBL 15 March 2024

#Reinstituted user_inputs module with non-dynamic sidebar. 

#MBL 27 March 2024

#Implemented new visualizations (charge site analysis, abundance summaries)
#Adjusted visualization layouts - fragment map and interactive annotated spectrum
# now live in the "identification" tab, separate from other visualizations
#Adjusted fragment map sizing, implemented flag plotting for protein spectra
#Interested in adding an interactive element to the fragment map, allowing selection
# of a particular ion to visualize - a sort of shortcut for the residues you're interested in
#Considering adding color-coding to symbol encodings for modifications to include in
# fragment map without taking up too much space.

#MBL 09 April 2024

# Quick updates to fragment map plotting for readability and legend structure

# Additionally, began separation of protein and nucleic acid apps.
# Rather than one unified app working for both, we have decided to break them apart
# to make their use cases clearer. This is now the app.R file for the newly
# minted Oligo-SAFARI (name pending)

#Packages=======================================================================

require(shiny)
require(data.table)
require(magrittr)
require(stringr)
require(IsoSpecR)
require(rawrr)
require(tools)
require(vroom)
require(nnls)
require(ggthemes)
require(grid)
require(scales)
require(ggplot2)
require(plotly)


#Modules========================================================================

#setwd("D:/UT Austin/Research/Shiny App/20240328")
source("Prerequisites and Helpers.R")
source("User Inputs Module.R")
source("Sequence Parsing Module.R")
source("Chemical Compositions Module.R")
source("Mass Generation Module.R")
source("mz Generation Module.R")
source("Data Readin Module.R")
source("Fragment Identification Module.R")
source("Visualizations Module.R")


#User Interface=================================================================

main_panel_layouts <- tabsetPanel(
  type = "tabs",
  id = "which_tab",
  
  tabPanel("Sequence Parsing", value = "parsing_tab",
    parsing_UI("parsing")
  ),
  tabPanel("Mass Data", value = "mass_tab",
    mass_UI("mass")
  ),
  tabPanel("m/z Data", value = "mz_tab",
    mz_UI("mz")
  ),
  tabPanel("Spectrum Read-In", value = "data_readin_tab",
    data_readin_UI("data_readin")
  ),
  tabPanel("Identifications", value = "identificaiton_tab",
    visualization_UI("visualization")$identification_tab
  ),
  tabPanel("Visualization", value = "visualization_tab",
    visualization_UI("visualization")$visualization_tab
  )#,
  # tabPanel("Citation", value = "citation_tab",
  #   citation_UI("visualization")
  # )
)


ui <- fluidPage(
  
  #titlePanel("Oligo - SAFARI"),
  titlePanel("Nucleo-SAFARI"),
  
  sidebarLayout(
    
    sidebarPanel(width = 3,
      user_inputs_UI("user_inputs")
    ),
    
    mainPanel(
      main_panel_layouts
    ) 
    
  )
  
  
)


#Server=========================================================================

server <- function(input, output, session){
  
  #Modules
  user_inputs_module <- callModule(user_inputs_server, "user_inputs"
  )
  
  parsing_module <- callModule(parsing_server, "parsing",
    user_inputs = user_inputs_module
  )
  
  composition_module <- callModule(composition_server, "composition",
    user_inputs = user_inputs_module,
    parsing = parsing_module
  )

  mass_module <- callModule(mass_server, "mass",
    user_inputs = user_inputs_module,
    parsing = parsing_module,
    compositions = composition_module
  )

  mz_module <- callModule(mz_server, "mz",
    user_inputs = user_inputs_module,
    mass = mass_module
  )

  data_readin_module <- callModule(data_readin_server, "data_readin",
    user_inputs = user_inputs_module
  )
  
  identification_module <- callModule(identification_server, "identification",
    user_inputs = user_inputs_module,
    mz = mz_module,
    data_readin = data_readin_module
  )

  visualization_module <- callModule(visualization_server, "visualization",
    user_inputs = user_inputs_module,
    parsing = parsing_module,
    data_readin = data_readin_module,
    identification = identification_module
  )
  
  # citation_module <- callModule(citation_server, "citation",
  # )
  
  ##Dynamic UI Components
  # observeEvent(input$which_tab, {
  #   updateTabsetPanel(inputId = "which_sidebar", selected = input$which_tab)
  # })
  #For now limit UI dynamism to avoid "removing" user inputs when tab changes.
  #I think it has to do with namespacing and calling the same module UI in multiple spots
  
  
}


#RUN APP========================================================================

shinyApp(ui = ui, server = server)

#End Document===================================================================