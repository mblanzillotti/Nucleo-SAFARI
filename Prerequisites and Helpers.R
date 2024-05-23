#===============================================================================

#R Shiny application for automated m/z domain fragment identification
#of proteins and nucleic acids

# MBL 14 March 2024

#Prerequisites and Helpers Module

#This module contains element and symbol data as well as helpers for sequence parsing
#Used downstream for mass determination of protein and nucleic acid sequences.
#As of 11 March 2024, digit modifications are not allowed, but will be reenabled.

#Packages=======================================================================

require(data.table)
require(magrittr)


#Prerequisites==================================================================

#Note: prerequisite variable (e.g. element_masses, etc.) may need to be
#included in each module separately. Test in the small scale with
#generation of monoisotopic masses

#Element Masses: contains the monoisotopic masses of elements that will be considered
#for mass calculation. Missing elements will not be considered. Likewise, elements
#will need to be added as an additional column in the symvols tables for each compound

#I apologize they are not particularly readable like tribbles would be,
#might be smart to put these in a separate file in the future, I just don't know how to do that yet

element_masses <- data.table(
  element = c("C", "H", "N", "O", "P", "S", "Na", "K", "Fe"),
  mono_mass = c(12.00000, 1.00783, 14.00307, 15.99492, 30.97376, 31.97207, 22.98977, 38.96370, 55.934940)
)


#Protein symbols containds predefined symbols (including amino acids and common modifications)
#for use in sequecnce parsing an mass calculation.
#The protein backbone is notated in this app by three symbols of the form: nXc
#where "n" corresponds to the amide N, "X" corresponds to the alpha carbon and side chain,
#and "c" corresponds to the carbonyl carbon. Modifications can be added to any
#of these elements by adding the mod in parentheses "()" after the appropriate symbol.

#e.g., a methylated lysine could be: K(Methyl), K(C1H2), or K(+14.01)
#Note that digits are a work in progress as of 11 Mar 2024

protein_symbols <- data.table(
  symbol =  c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R",  "S", "T", "V", "W", "Y", "Acetyl", "Carbamyl", "Methyl", "Deoxid", "Heme", "n", "c"),
  C =  c(2, 2, 3, 4, 8, 1, 5, 5,  5,  5,  4, 3, 4, 4, 5,  2, 3, 4, 10, 8, 2, 1, 1, 0,  34, 0, 1),
  H =  c(4, 4, 4, 6, 8, 2, 6, 10, 11, 10, 8, 5, 6, 7, 11, 4, 6, 8, 9,  8, 2, 2, 2, 0,  32, 1, 0),
  N =  c(0, 0, 0, 0, 0, 0, 2, 0,  1,  0,  0, 1, 0, 1, 3,  0, 0, 0, 1,  0, 0, 1, 0, 0,  4,  1, 0),
  O =  c(0, 0, 2, 2, 0, 0, 0, 0,  0,  0,  0, 1, 0, 1, 0,  1, 1, 0, 0,  1, 1, 1, 0, -1, 4,  0, 1),
  S =  c(0, 1, 0, 0, 0, 0, 0, 0,  0,  0,  1, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0,  0,  0, 0),
  Fe = c(0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0,  1,  0, 0)
) %>% 
  melt(
    measure.vars = c("C", "H", "N", "O", "S", "Fe"),
    variable.name = "element", variable.factor = F,
    value.name = "amount"
  ) %>% .[, amount := as.integer(amount)]


#Nucleic symbols contains a separate set of symbols corresponding to nucleic acid
#nomenclature and symbol parsing.
#The phosphodiester backbone is encoded in this app of the for sBb, where
#the sugar and any modifications to the 1', 2', 3', 4', or 5' sugar is encoded in "s",
#the nitrogenous base is encoded in "B", and the corresponding backbone is encoded in "b".

#Modifications to any of these three elements or likewise to any carbons can be
#encoded by adding the modification in parentheses after the residue. Some residues
#are implied, such as deoxyribose and phosphate backbones by default.

#modifications to the 3' and 5' oxygens (default) should be specified using the
#form: s(5'Phos), etc. which will override the default behavior of oxygen interpolation.

nucleic_symbols <- data.table(
  symbol =  c("r", "d", "m", "+", "A", "C", "G", "H", "T", "I", "L", "U", "D", "S", "Y",  "*", "p", "Acetyl", "Methyl", "Dimethyl", "Phos", "MethylPhos", "DimethylPhos", "DiPhos", "TriPhos"),
  C =  c(5, 5, 6, 6, 5, 4, 5, 6, 5, 5, 0, 4, 4, 4, 16, 0, 0, 2, 1, 2, 0, 1, 2, 0, 0),
  H =  c(7, 7, 9, 7, 4, 4, 4, 6, 5, 3, 0, 3, 5, 3, 19, 1, 1, 2, 2, 4, 2, 4, 6, 3, 4),
  N =  c(0, 0, 0, 0, 5, 3, 5, 5, 2, 4, 0, 2, 2, 2, 6,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  O =  c(2, 1, 2, 2, 0, 1, 1, 1, 2, 1, 0, 2, 2, 1, 5,  1, 2, 1, 0, 0, 4, 4, 4, 7, 10),
  P =  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 0, 0, 0, 1, 1, 1, 2, 3),
  S =  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  Na = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  K =  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
) %>% 
  melt(
    measure.vars = c("C", "H", "N", "O", "P", "S", "Na", "K"),
    variable.name = "element", variable.factor = F,
    value.name = "amount"
  ) %>% .[, amount := as.integer(amount)]
#H N7 methyl G, where there is a fixed positive charge
#(and a net +3 H rather than typical 2, but considered as 2 for charge accounting)


#Functions======================================================================

interpret_formula <- function(encoding){
  data.table(
    element = str_extract_all(encoding, "[:alpha:]{1}[:lower:]?\\-?\\d+") %>% 
      unlist() %>% 
      str_extract("[:alpha:]{1}[:lower:]?"),
    amount = str_extract_all(encoding, "[:upper:]{1}[:lower:]?\\-?\\d+") %>% 
      unlist() %>% 
      str_extract("\\-?\\d+") %>% 
      as.integer()
  )[, .(amount = sum(amount)), by = element
  ][element_masses[, element], on = .(element)
  ][is.na(amount), amount := 0L]
}
#Helper which generates a tibble of the format [element, amount] from an input chemical formula


interpret_digits <- function(encoding){
  str_remove_all(encoding, "[\\(\\)]") %>% 
    str_extract("[\\-\\+]?\\d*\\.?\\d*") %>% 
    as.numeric()
}
#Helper that extracts a net loss or gain in mass from a modification


interpret_symbol <- function(encoding, symbol_composition){
  symbol_composition[symbol == encoding] %>% 
    .[, !"symbol"]
}
#Helper that extracts the associated chemical formula from a symbol in "symbol_compositions"


determine_composition <- function(encoding, symbol_composition){
  
  if(str_detect(encoding, "[:alpha:]") && str_detect(encoding, "[:digit:]")){
    #Letters and numbers means formula
    composition <- interpret_formula(encoding)
    
  } else if(str_detect(encoding, "[:graph:]") && !str_detect(encoding, "[:digit:]")){
    
    composition <- interpret_symbol(encoding, symbol_composition)
    
  } else{ 
    #If it's digits only, that will be handled by determine_digits separately
    composition <- NA
    
  }
  
  return(composition)
  
}
#Function to interpret and return the composition for either a symbol or formal composition


# determine_digits <- function(encoding){
#   
#   if(!str_detect(encoding, "[:alpha:]") && str_detect(encoding, "[:digit:]")){
#     
#     composition <- interpret_digits(encoding)
#     
#   } else {
#     
#     composition <- 0L
#     
#   }
#   
#   return(composition)
#   
# }
#Function to interpret and return digit modifications. Temporarily disabled