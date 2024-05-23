#===============================================================================

#R Shiny application for automated m/z domain fragment identification
#of proteins and nucleic acids

# MBL 14 March 2024

#Sequence Parsing Module

#This module contains two modules which interpret the protein or nucleic acid
#sequence inputs provided by the app user and parses them by amide-N/residue/carbonyl
#or sugar/base/babckbone, respectively. This format in necessary for composition
#and mass determination downstream

#A user-selected input of "biomolecule" dictates which module is called to parse the 
#user-input sequence.

#These parsing functions require both stringr (for string matching) and 
#the prerequisite file for accurate interpretation


#Packages=======================================================================

require(data.table)
require(stringr)
require(magrittr)


#Functions======================================================================

#note, any additionally defined single-symbols for bases MUST be added to the initial
#regular expression AND to the nitbase and nitbase mod sections of the code.
#Please play nice with the regexp
parse_nucleic_input <- function(input_seq){
  
  parsed <- str_extract_all(input_seq, "[rmd\\+]?(\\([12345]?'?[\\+\\-]?[:alnum:]+[.-]?[:alnum:]*\\)){0,5}[ACGTUIYSDLH]{1}(\\([\\+\\-]?[:alnum:]+[.-]?[:alnum:]+\\))?[p*]?(\\([\\+\\-]?[:alnum:]+[.-]?[:alnum:]+\\))?") %>% 
    unlist()
  
  parsed <- data.table(
    IDT_parsed = parsed,
    position = 1:length(parsed),
    sugar = str_extract(parsed, "^[rmd\\+]"),
    sugar_mod = str_extract(parsed, "[rmd\\+](\\([12345]?'?[\\+\\-]?[:alnum:]+[.-]?[:alnum:]+\\)){1,5}") %>%
      str_remove("[rmd\\+]"),
    nitbase = parsed %>%
      str_remove("[rmd\\+]?(\\([12345]?'?[\\+\\-]?[:alnum:]+[.-]?[:alnum:]*\\)){0,5}") %>%
      str_remove("(\\([\\+\\-]?[:alnum:]+[.-]?[:alnum:]+\\))?[p*]?(\\([\\+\\-]?[:alnum:]+[.-]?[:alnum:]+\\))") %>%
      str_extract("[ACGTUYSIDHL]{1}"),
    nitbase_mod = str_extract(parsed, "[ACGTUYISDLH]{1}(\\([\\+\\-]?[:alnum:]+[.-]?[:alnum:]+\\))") %>%
      str_remove("[ACGTUIYSDLH]") %>% 
      str_remove_all("[\\(\\)]"),
    backbone = str_extract(parsed, "[p*]"),
    backbone_mod = str_extract(parsed, "[p*](\\([\\+\\-]?[:alnum:]+[\\].?-][:alnum:]+\\))") %>%
      str_remove("[p*]") %>% 
      str_remove_all("[\\(\\)]")
  ) %>%
    .[is.na(sugar), sugar := "d"] %>%
    .[is.na(backbone), backbone := "p"] %>%
    .[nrow(.), backbone := NA] %>%
    .[, `1'` := str_extract(sugar_mod, "\\(1'[\\+\\-]?[:alnum:]+[.-]?[:alnum:]+\\)") %>% str_remove("1'") %>% str_remove_all("[\\(\\)]")] %>%
    #maybe it's possible to do this with lapply?
    .[, `2'` := str_extract(sugar_mod, "\\(2'[\\+\\-]?[:alnum:]+[.-]?[:alnum:]+\\)") %>% str_remove("2'") %>% str_remove_all("[\\(\\)]")] %>%
    .[, `3'` := str_extract(sugar_mod, "\\(3'[\\+\\-]?[:alnum:]+[.-]?[:alnum:]+\\)") %>% str_remove("3'") %>% str_remove_all("[\\(\\)]")] %>%
    .[is.na(`3'`), `3'` := "O1"] %>%
    .[`3'` == "O1" & position == nrow(.), `3'` := "O1H1"] %>%
    .[, `4'` := str_extract(sugar_mod, "\\(4'[\\+\\-]?[:alnum:]+[.-]?[:alnum:]+\\)") %>% str_remove("4'") %>% str_remove_all("[\\(\\)]")] %>%
    .[, `5'` := str_extract(sugar_mod, "\\(5'[\\+\\-]?[:alnum:]+[.-]?[:alnum:]+\\)") %>% str_remove("5'") %>% str_remove_all("[\\(\\)]")] %>%
    .[is.na(`5'`), `5'` := "O1"] %>%
    .[`5'` == "O1" & position == 1, `5'` := "O1H1"] %>%
    .[, !c("IDT_parsed", "sugar_mod")] %>%
    melt(
      measure.vars = c("sugar", "nitbase", "nitbase_mod", "backbone", "backbone_mod", "1'", "2'", "3'", "4'", "5'"),
      variable.name = "item", variable.factor = F,
      value.name = "symbol"
    ) %>%
    .[, item := factor(item, levels = c("5'", "sugar", "1'", "2'", "4'", "nitbase", "nitbase_mod", "3'", "backbone", "backbone_mod"))] %>%
    .[order(position, item)] %>%
    .[!is.na(symbol)]
  
  
  return(parsed)
  
}


# parse_protein_input <- function(input_sequence){
#   
#   parsed <- str_extract_all(input_sequence, "(n(\\([\\+\\-]?[:alnum:]+[\\.]?[\\d]{0,}\\)))?[:upper:](\\([\\+\\-]?[:alnum:]+[\\.]?[\\d]{0,}\\))?(c(\\([\\+\\-]?[:alnum:]+[\\.]?[\\d]{0,}\\)))?") %>% 
#     unlist()
#   
#   parsed <- data.table(
#     protein_parsed = parsed,
#     position = 1:length(parsed),
#     N_term = ifelse(str_detect(parsed, "n"), str_extract(parsed, "n(\\([\\+\\-]?[:alnum:]+[\\.]?[\\d]{0,}\\))"), "n"),
#     #N_term_mod = str_extract(N_term, "(\\([\\+\\-]?[:alnum:]+[\\.]?[\\d]{0,}\\))") %>% str_remove_all("[\\(\\)]"),
#     residue = str_remove(parsed, "(n(\\([\\+\\-]?[:alnum:]+[\\.]?[\\d]{0,}\\)))?") %>%
#       str_extract("[:upper:]((\\([\\+\\-]?[:alnum:]+[\\.]?[\\d]{0,}\\)))?"),
#     #residue_mod = str_extract(residue, "(\\([\\+\\-]?[:alnum:]+[\\.]?[\\d]{0,}\\))") %>% str_remove_all("[\\(\\)]"),
#     C_term = ifelse(str_detect(parsed, "c"), str_extract(parsed, "c((\\([\\+\\-]?[:alnum:]+[\\.]?[\\d]{0,}\\))?)"), "c")#,
#     #C_term_mod = str_extract(C_term, "(\\([\\+\\-]?[:alnum:]+[\\.]?[\\d]{0,}\\))") %>% str_remove_all("[\\(\\)]")
#   )[, N_term_mod := str_extract(N_term, "(\\([\\+\\-]?[:alnum:]+[\\.]?[\\d]{0,}\\))") %>% str_remove_all("[\\(\\)]")
#     ][, N_term := str_extract(N_term, "n")
#     ][, C_term_mod := str_extract(C_term, "(\\([\\+\\-]?[:alnum:]+[\\.]?[\\d]{0,}\\))") %>% str_remove_all("[\\(\\)]")
#     ][, C_term := str_extract(C_term, "c")
#     ][, residue_mod := str_extract(residue, "(\\([\\+\\-]?[:alnum:]+[\\.]?[\\d]{0,}\\))") %>% str_remove_all("[\\(\\)]")
#     ][, residue := str_extract(residue, "[:upper:]")
#     ][position == 1 & is.na(N_term_mod), N_term_mod := "H1",
#     ][position == max(position) & is.na(C_term_mod), C_term_mod := "O1H1"
#     ] %>% 
#     melt(
#       measure.vars = c("N_term", "N_term_mod", "residue", "residue_mod", "C_term", "C_term_mod"),
#       variable.name = "item", variable.factor = F,
#       value.name = "symbol"
#     )
#     
#   parsed <- parsed[, item := factor(item, levels = c("N_term", "N_term_mod", "residue", "residue_mod", "C_term", "C_term_mod"))
#     ][order(position, item)
#     ][!is.na(symbol)
#     ][, .(position, item, symbol)]
#   
#   
#   return(parsed)
#   
# }

#Module Server==================================================================

parsing_server <- function(input, output, session, user_inputs){
  
  precursor_parsed <- reactive({
    
    parse_nucleic_input(user_inputs$input_sequence())
    
    # switch(user_inputs$biomolecule(),
    #   nucleic_acid = parse_nucleic_input(user_inputs$input_sequence()),
    #   protein = parse_protein_input(user_inputs$input_sequence())
    # )
    
  })
    
  output$precursor_parsed_table <- renderTable({
    
    req(user_inputs$input_sequence())
    
    precursor_parsed() %>% 
      dcast(position~item)
    
  })
  

  return(
    list(
      precursor_parsed = reactive({ precursor_parsed() })
    )
  )
  
}

#Module UI======================================================================

parsing_UI <- function(id){
  
  ns <- NS(id)
  
  
  return(
    tagList(
      tableOutput(ns("precursor_parsed_table"))
    )
  )

}


