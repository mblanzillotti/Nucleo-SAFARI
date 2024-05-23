#===============================================================================

#R Shiny application for automated m/z domain fragment identification
#of proteins and nucleic acids

# MBL 14 March 2024

#Fragment Composition Module

#Generate fragments based on precursor composition, biomolecule, and neutral loss
#toggle.


#Packages=======================================================================

require(data.table)
require(stringr)


#Functions======================================================================

generate_precursor_composition <- function(parsed_sequence, symbol_composition){
  
  precursor_composition <- parsed_sequence[, composition := lapply(symbol, determine_composition, symbol_composition = symbol_composition)
  ][, rbindlist(composition), by = .(position, symbol, item)
  ][, filter_key := (amount == 0), by = .(element)
  ][, filter_key := all(filter_key), by = .(element)
  ][filter_key == F
  ][, .(position, symbol, item, element, amount)]
  
  
  return(precursor_composition)
  
}


generate_nucleic_fragment_compositions <- function(precursor_composition, base_losses = T){
  
  fragment_compositions <- precursor_composition 
  fragment_compositions <- fragment_compositions[str_detect(item, "[124]|(sugar)|(nitbase)"), ion_type := "nside"
  ][is.na(ion_type), ion_type := item
  ][, .(amount = sum(amount)), by = .(position, ion_type, element)
    #][str_detect(ion_type, "nside"), amount := sum(amount), by = .(element, position)
  ][, acc_5 := cumsum(amount), by = .(element)
  ][, acc_3 := rev(cumsum(rev(amount))), by = .(element)
  ][, .(position, ion_type, acc_5, acc_3, element)] %>% 
    melt(
      measure.vars = c("acc_5", "acc_3"),
      variable.name = "terminus",
      value.name = "amount"
    )
  
    #5'-terminal fragments (a/b/c/d)
  fragment_compositions <- fragment_compositions[str_detect(ion_type, "5'") & str_detect(terminus, "acc_5"), ion_type := "d"
  ][str_detect(ion_type, "nside") & str_detect(terminus, "acc_5"), ion_type := "a"
  ][str_detect(ion_type, "3'") & str_detect(terminus, "acc_5"), ion_type := "b"
  ][str_detect(ion_type, "backbone") & str_detect(terminus, "acc_5"), ion_type := "c"
    #3'-terminal fragments (w/x/y/z)
  ][str_detect(ion_type, "5'") & str_detect(terminus, "acc_3"), ion_type := "y"
  ][str_detect(ion_type, "nside") & str_detect(terminus, "acc_3"), ion_type := "z"
  ][str_detect(ion_type, "3'") & str_detect(terminus, "acc_3"), ion_type := "w"
  ][str_detect(ion_type, "backbone") & str_detect(terminus, "acc_3"), ion_type := "x"
    #Position corrections
  ][str_detect(ion_type, "d"), position := position-1
  ][str_detect(ion_type, "[yz]"), position := position-1
  ][, .(position, ion_type, element, amount)
  ][position > 0]
  
  
  if(base_losses){
  #base loss fragments
    base_loss_fragment_comps <- fragment_compositions
    base_loss_fragment_comps <- base_loss_fragment_comps[, adj_position := position + 1*str_detect(ion_type, "[wxyz]")
                                                       #The above := modification alters fragment_comps in place, likely because b_l_f_c is only a reference at this point
    ][precursor_composition[item == "nitbase"], on = .(adj_position = position, element)
    ][, amount := amount - i.amount
    ][, .(position, ion_type, element, amount, adj_position)
    ][, ion_type := paste0(ion_type, "-B")]
  
    fragment_compositions <- rbind(fragment_compositions, base_loss_fragment_comps)[, .(position, ion_type, element, amount)]
  
  }
  
  hydrogen_adjustments <- data.table(
    ion_type = c("a", "a-B", "b", "b-B", "c", "c-B", "d", "d-B", "w", "w-B", "x", "x-B", "y", "y-B", "z", "z-B"),
    element = c("H"),
    amount = c(0, -2, 1, 1, 0, 0, 2, 2, 1, 1, -1, -1, 1, 1, 2, 2)
  )
  
  fragment_compositions <- merge(
    fragment_compositions, hydrogen_adjustments, by = c("ion_type", "element"), all.x = T
  )[!is.na(amount.y), amount.x := amount.x + amount.y
  ][, amount := amount.x
  ][, !c("amount.x", "amount.y")
  ][order(ion_type, position, element)]
  
  
  return(fragment_compositions)
  
}

# 
# generate_protein_fragment_compositions <- function(precursor_composition, neutral_losses = T){
#   
#   fragment_compositions <- precursor_composition[, ion_type := str_extract(item, "(N_term)|(residue)|(C_term)")
#     ][, .(amount = sum(amount)), by = .(position, ion_type, element)
#     ][, acc_N := cumsum(amount), by = .(element)
#     ][, acc_C := rev(cumsum(rev(amount))), by = .(element)
#     ][, .(position, ion_type, acc_N, acc_C, element)] %>% 
#     melt(
#       measure.vars = c("acc_N", "acc_C"),
#       variable.name = "terminus",
#       value.name = "amount"
#     )
#   
#     #N-terminal Fragments (a/b/c)
#   fragment_compositions <- fragment_compositions[str_detect(ion_type, "N_term") & str_detect(terminus, "acc_N"), ion_type := "c"
#   ][str_detect(ion_type, "residue") & str_detect(terminus, "acc_N"), ion_type := "a"
#   ][str_detect(ion_type, "C_term") & str_detect(terminus, "acc_N"), ion_type := "b"
#     #C-terminal fragments (x/y/z)
#   ][str_detect(ion_type, "N_term") & str_detect(terminus, "acc_C"), ion_type := "y"
#   ][str_detect(ion_type, "residue") & str_detect(terminus, "acc_C"), ion_type := "z"
#   ][str_detect(ion_type, "C_term") & str_detect(terminus, "acc_C"), ion_type := "x"
#     #Position corrections
#   ][str_detect(ion_type, "[bc]"), position := position-1
#   ][str_detect(ion_type, "[yz]"), position := position-1#maybe, not sure about this w/o testing
#   ][, .(position, ion_type, element, amount)
#   ][position > 0]
#   
#   if(neutral_losses){
#     #neutral loss fragments
#     neutral_losses <- data.table(
#       ion_type = c("b", "b", "y", "y"),
#       element = c("H", "O", "H", "O"),
#       amount = c(2, 1, 2, 1)
#     )
#     
#     neutral_loss_fragments <- merge(fragment_compositions[str_detect(ion_type, "[by]")],
#       neutral_losses,
#       by = c("element", "ion_type"), allow.cartesian = T, all.x = T
#     )[is.na(amount.y), amount.y := 0
#       ][, amount := amount.x-amount.y
#       ][, .(ion_type, position, element, amount)
#       ][order(ion_type, position, element)
#       ][, ion_type := paste0(ion_type, "-H2O") #This will break down with more neutral losses...
#       ][, .SD[all(amount >= 0)], by = .(ion_type, position), .SDcols = c("amount", "element")]
# 
#     fragment_compositions <- rbind(fragment_compositions, neutral_loss_fragments)
#     
#   }
#   
#   hydrogen_adjustments <- data.table(
#     ion_type = c("a", "b", "b-H2O", "c", "x", "y", "y-H2O", "z"),
#     element = c("H"),
#     amount = c(-2, 0, 0, 1, 0, 2, 2, 1)
#   )
#   
#   fragment_compositions <- merge(
#     fragment_compositions, hydrogen_adjustments, by = c("ion_type", "element"), all.x = T
#   )[!is.na(amount.y), amount.x := amount.x + amount.y
#   ][, amount := amount.x
#   ][, !c("amount.x", "amount.y")
#   ][order(ion_type, position, element)
#   ][position %in% unique(precursor_composition[symbol == "P"]$position) & ion_type == "y" & element == "H", amount := amount - 1
#   ][position != max(precursor_composition$position)]
#   
#   
#   return(fragment_compositions)
#   
# }


#Module Server==================================================================

composition_server <- function(input, output, session, user_inputs, parsing){
  
  precursor_composition <- reactive({
    
    req(user_inputs$input_sequence())
    
    generate_precursor_composition(
      parsing$precursor_parsed(),
      nucleic_symbols
    )
    
    # switch(user_inputs$biomolecule(),
    #   nucleic_acid = generate_precursor_composition(
    #     parsing$precursor_parsed(),
    #     nucleic_symbols
    #   ),
    #   protein = generate_precursor_composition(
    #     parsing$precursor_parsed(), 
    #     protein_symbols
    #   )
    # )
    
  })
  
  fragment_compositions <- reactive({
    
    req(user_inputs$input_sequence())
    
    generate_nucleic_fragment_compositions(
      precursor_composition(),
      user_inputs$neutral_loss_toggle()
    )
    
    # switch(user_inputs$biomolecule(),
    #   nucleic_acid = generate_nucleic_fragment_compositions(
    #     precursor_composition(),
    #     user_inputs$neutral_loss_toggle()
    #   ),
    #   protein = generate_protein_fragment_compositions(
    #     precursor_composition(),
    #     user_inputs$neutral_loss_toggle()
    #   )
    # )
    
  })
  
  
  return(
    list(
      precursor_composition = reactive({ precursor_composition() }),
      fragment_compositions = reactive({ fragment_compositions() })
    )
  )
  
}


#Module UI======================================================================

composition_UI <- function(id){
  
  ns <- NS(id)
  
  #Don't think there are any outputs or UI elements for this module...
  
  
}