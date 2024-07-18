#===============================================================================

#R Shiny application for automated m/z domain fragment identification
#of proteins and nucleic acids

# MBL 14 March 2024

#Visualization Module

#Visualize fragment identification data


#Packages=======================================================================

require(data.table)
require(magrittr)
require(ggplot2)
require(plotly)
require(grid)
require(scales)
require(ggthemes)

#Functions======================================================================


scale_color_publication <- function(...){
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

scale_fill_publication <- function(...){
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

#Protein Colors
# protein_colors <- c(
#   "a" = "#B4D496",
#   "b" = "#A0BEEE",
#   "c" = "#FF8A8A",
#   "x" = "#00B050",
#   "y" = "#3764AA",
#   "z" = "#C00000"
# )

#Nucleic Acid Colors
nucleic_colors <- c(
  "a"   = "#33889F",
  "a-B" = "#4BACC6",
  "b"   = "#B22112",
  "b-B" = "#D42816",
  "c"   = "#9BBB59",
  "c-B" = "#BFDF37",
  "d"   = "#8064A2",
  "d-B" = "#8E81D8",
  "w"   = "#2B4C73",
  "w-B" = "#1F3651",
  "x"   = "#A13B39",
  "x-B" = "#722A28",
  "y"   = "#6B8335",
  "y-B" = "#55682A",
  "z"   = "#584470",
  "z-B" = "#3E304E"
)

#Theme
theme_publication <- function(base_size=14) {
  (theme_foundation(base_size=base_size)
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(1.2), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.title.y.right = element_text(vjust=2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(), 
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           axis.ticks.length = unit(-0.1, "cm"),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "bottom",
           legend.direction = "horizontal",
           legend.key.size= unit(0.2, "cm"),
           legend.spacing = unit(0, "cm"),
           legend.title = element_text(face="italic"),
           plot.margin = unit(c(2,2,2,2),"mm"),
           strip.background = element_blank(),
           strip.text = element_text(face="bold")
   ))
  
}


generate_ion_label <- function(ion_type, position, charge){
  
  paste0(
    ion_type, 
    "[", position, "]",
    "^{\n", 
    rep("~", ceiling(log10(position+1))+1) %>% paste0(collapse = ""),
    abs(charge),
    " * \"", ifelse(charge > 0, "+", "-"), "\"\n}"
  )
  
}


#Plots
plot_ppm_error <- function(fragment_identifications){
  
  plot <- ggplot()+
    geom_point(data = fragment_identifications,
      aes(x = centroid_mz, y = ppm_error, color = parent_type), alpha = 0.25
    )+
    scale_color_manual(values = nucleic_colors, limits = names(nucleic_colors), drop = F)+
    theme_publication(base_size = 12)+
    theme(
      panel.grid.major = element_blank()
    )+
    guides(color = guide_legend(nrow = 2, position = "top", override.aes = list(alpha = 1)))+
    labs(
      x = expression("Experimental"~italic("m/z")),
      y = expression(italic("m/z")~"Error (ppm)"),
      color = "Ion Type"
    )
  
  
  return(plot)
  
}


plot_abundance_summary <- function(fragment_identifications, spectrum){
  
  TIC <- sum(spectrum$centroid_data$centroid_intensity)
  
  intensity_data <- fragment_identifications[, .(total_intensity = sum(centroid_intensity)), by = .(parent_type)]
  count_data <- fragment_identifications[, .(parent_type, position)] %>% 
                unique() %>% 
                .[, .(N_fragments = .N), by = .(parent_type)]
  
  figure_data <- merge(
    intensity_data, count_data,
    by = c("parent_type")
  )
  
  plot <- ggplot()+
    geom_col(data = figure_data,
      aes(x = parent_type, y = total_intensity/TIC*100, fill = parent_type)
    )+
    geom_text(data = figure_data,
      aes(x = parent_type, y = total_intensity/TIC*100+0.05*max(total_intensity)/TIC*100, label = N_fragments),
      color = "darkgrey", size = 4
    )+
    scale_fill_manual(values = nucleic_colors, limits = names(nucleic_colors), drop = F)+
    theme_publication(base_size = 12)+
    theme(
      panel.grid.major = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1.0),
      legend.position = "none"
    )+
    labs(
      x = "Ion Type",
      y = "Ion Current (%)"
    )
    
  
  return(plot)
  
}


plot_abundance_by_position <- function(fragment_identifications, spectrum){
  
  TIC <- sum(spectrum$centroid_data$centroid_intensity)
  
  plot <- ggplot()+
    geom_col(data = fragment_identifications[, .(pos_sum_intensity = sum(predicted_intensity)), by = .(position, parent_type)
                      ][, adj_position := position
                      ][str_detect(parent_type, "[wxyz]"), adj_position := position + 1],
      aes(x = adj_position, y = pos_sum_intensity/TIC*100, fill = parent_type)
    )+
    scale_fill_manual(values = nucleic_colors, limits = names(nucleic_colors), drop = F)+
    theme_publication(base_size = 12)+
    theme(
      panel.grid.major = element_blank()
    )+
    guides(fill = guide_legend(nrow = 2, position = "top", override.aes = list(drop = F)))+
    labs(
      x = "Position",
      y = "Ion Current (%)",
      fill = "Ion Type"
    )
  
    
  return(plot)
    
}


plot_annotated_spectrum <- function(fragment_identifications, spectrum, mass_range = c(200, 2000)){
  
  scalar <- spectrum$profile[intensity == max(intensity)]$intensity[1] %>% log10() %>% floor()
  
  plot <- ggplot()+
    geom_line(data = spectrum$profile_data,
      aes(x = mz, y = intensity/10^scalar), linewidth = 0.5, color = "black"
    )+
    geom_line(data = fragment_identifications[, .(predicted_intensity = sum(predicted_intensity)), by = .(parent_type, position, charge, centroid_mz)
                                              ][, TPZ := paste(parent_type, position, charge, sep = "|")],
      aes(x = centroid_mz, y = predicted_intensity/10^scalar, group = TPZ, color = parent_type)
    )+
    geom_point(data = fragment_identifications,
      aes(x = theoretical_mz, y = predicted_intensity/10^scalar, color = parent_type), size = 0.75
    )+
    # geom_text(data = fragment_identifications[, .SD[max(centroid_intensity) == centroid_intensity], by = .(parent_type, position, charge), .SDcols = c("centroid_mz", "centroid_intensity")
    #                    ][, ion_label := paste0(parent_type, "[", position, "]", "^{\n", abs(charge), " * \"-\"\n}")],
    #   aes(x = centroid_mz, y = centroid_intensity*1.05/10^scalar, color = parent_type, label = ion_label),
    #   parse = T, size = 5, hjust = 0.5, check_overlap = T, vjust = 0.0
    # )+
    # scale_color_manual(
    #   values = switch(biomolecule,
    #     nucleic_acid = nucleic_colors,
    #     protein = protein_colors
    #   )
    # )+
    scale_color_manual(values = nucleic_colors, limits = names(nucleic_colors), drop = F)+
    theme_publication(base_size = 14)+
    theme(
      panel.grid.major = element_blank(),
      legend.position = "none"
    )+
    labs(
      x = "m/z",
      y = "Intensity"
    )
    # labs(
    #   x = expr(italic("m/z")),
    #   y = expr("Intensity"~"(10"^(!!scalar)*")")
    # )
  #note these axis labels won't work in a plotly context
  #Nor will the nice parsed fragment ion labels, which is a shame.
  
  
  return(plot)
  
}


plot_fragment_map <- function(fragment_identifications, precursor_parsed){
  
  coordinate_scaffold <- data.table(
      ordering = rep(c(1, 2, 3), 8),
      x_coord = c(rep(c(0, 0, -0.15), 4), rep(c(0, 0, 0.15), 4)),
      y_coord = c(rep(c(0, 0.25, 0.25), 4), rep(c(0, -0.25, -0.25), 4)),
      parent_type = lapply(c("a", "b", "c", "d", "w", "x", "y", "z"), rep, times = 3) %>% unlist()
    ) %>% 
      rbind(
        data.table(
          ordering = rep(c(1, 2), 8),
          x_coord = c(rep(c(-0.16, -0.24), 4), rep(c(0.16, 0.24), 4)),
          y_coord = c(rep(c(0.25, 0.25), 4), rep(c(-0.25, -0.25), 4)),
          parent_type = lapply(c("a-B", "b-B", "c-B", "d-B", "w-B", "x-B", "y-B", "z-B"), rep, times = 2) %>% unlist()
        )
      )
  
  ion_coordinates <- data.table(
      x_adj = c(0, 0, 1, 1, 2, 2, 3, 3, 0, 0, 1, 1, 2, 2, 3, 3)*0.1 + 0.5-0.15,
      y_adj = (rep(c(0, 0, 0.6, 0.6, 1.2, 1.2, 1.8, 1.8), 2) - 0.9)*0.1,
      parent_type = c("a", "a-B", "b", "b-B", "c", "c-B", "d", "d-B", "w", "w-B", "x", "x-B", "y", "y-B", "z", "z-B")
    ) %>% 
      merge(coordinate_scaffold, by = "parent_type")
  
  
  ion_coordinates <- ion_coordinates[order(parent_type, ordering)
  ][, x_coord := x_coord + x_adj
  ][, y_coord := y_coord + y_adj]
  
  
  plot_data <- fragment_identifications[, .(position, parent_type)] %>% unique()
  
  plot_data <- ion_coordinates[plot_data, on = .(parent_type), allow.cartesian = T
  ][, position_x := (position %% 20) - 1
  ][position_x == -1 & str_detect(parent_type, "[abcd]"), position_x := 19
  ][, position_y := - (position %/% 20)
  ][position_x == 19 & str_detect(parent_type, "[abcd]"), position_y := position_y + 1
  ][, x_coord := x_coord + position_x
  ][, y_coord := y_coord + position_y]
  
  
  plot <- ggplot()+
    geom_path(data = plot_data[order(parent_type, position, ordering), T_P := paste0(parent_type, position)],
      aes(x = x_coord, y = y_coord, group = T_P, color = parent_type), linewidth = 0.5, linejoin = "mitre"
    )+
    geom_text(data = precursor_parsed[item == "sugar" | item == "nitbase", .(glyph = paste0(symbol, collapse = "")), by = position
        ][, x_pos := (position - 1) %% 20
        ][, y_pos := -(position) %/% 20 + 1],
      aes(x = x_pos, y = y_pos, label = glyph),
      size = 3.5
    )+
    geom_text(data = precursor_parsed[position %% 20 == 1, .(position)
      ][, y_pos := -position %/% 20 + 1
      ][, x_pos := -1] %>% unique(),
      aes(x = x_pos, y = y_pos, label = position),
      color = "grey", size = 3, vjust = 0.5, hjust = 1.0
    )+
    geom_text(data = precursor_parsed[position %% 20 == 0, .(position)
      ][, y_pos := -position %/% 20 + 1
      ][, x_pos := 20] %>% unique(),
      aes(x = x_pos, y = y_pos, label = position),
      color = "grey", size = 3, vjust = 0.5, hjust = 0.0
    )+
    scale_color_manual(values = nucleic_colors, limits = names(nucleic_colors), drop = F)+
    theme_void()+
    theme(
      panel.grid.major = element_blank(),
      axis.line = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.key.size= unit(0.2, "cm"),
      legend.spacing = unit(0, "cm"),
      legend.title = element_text(face="italic"),
      plot.margin = unit(c(2,2,2,2),"mm"),
      legend.text = element_text(size = 10)
    )+
    guides(color = guide_legend(nrow = 2, position = "top", override.aes = list(linewidth = 2, drop = F)))+
    coord_fixed(ratio = 1.1/1)+
    labs(
      color = "Ion Type"
    )
  
  
  return(plot)
  
}


plot_chargesite_localization <- function(fragment_identifications, precursor_parsed){
  
  charge_data <- fragment_identifications[, .(position, parent_type, centroid_intensity, predicted_intensity, charge)]
  
  charge_data <- charge_data[, terminal_type := ifelse(str_detect(parent_type, "[abcd]"), "L", "R")
  ][, .(predicted_intensity = sum(predicted_intensity)), by = .(terminal_type, position, charge)
  ][, charge_normalized_intensity := predicted_intensity/sum(predicted_intensity), by = .(terminal_type, position)
  ][str_detect(terminal_type, "R"), adj_position := position + 1
  ][str_detect(terminal_type, "L"), adj_position := position]
  
  charge_data <- precursor_parsed[, .(position)] %>% unique() %>% .[charge_data, on = .(position), all = T]
  #This does overwrite validated_fragments oddly with terminal_type, gotta properly make a copy somehow
  
  min_charge <- charge_data$charge %>% unique() %>% abs() %>% min() * sign(charge_data$charge %>% unique() %>% min())
  max_charge <- charge_data$charge %>% unique() %>% abs() %>% max() * sign(charge_data$charge %>% unique() %>% max())
    
  facet_labels <- c(
    "L" = "5' (a/b/c/d)",
    "R" = "3' (w/x/y/z)"
  )
  
  plot <- ggplot()+
    geom_col(data = charge_data,
             aes(x = adj_position, y = charge_normalized_intensity*100, fill = charge)
    )+
    facet_grid(
      terminal_type~., 
      labeller = labeller(terminal_type = facet_labels)
    )+
    theme_publication(base_size = 12)+
    theme(
      #axis.text.x = element_blank(),
      panel.grid.major = element_blank(),
      legend.key.width = unit(0.75, "in")
    )+
    scale_fill_gradientn(colours = rainbow(8), breaks = c(min_charge, max_charge))+
    scale_x_continuous(
      breaks = c(
        seq(1, max(precursor_parsed$position), 10),
        max(precursor_parsed$position)
      )
    )+
    scale_y_continuous(n.breaks = 3)+
    labs(
      x = "Position",
      y = "%",
      sec.y = "Fragment Ion Types",
      fill = "Charge State "
    )
  
  
  return(plot)
  
}


#Scoring
calculate_sequence_coverage <- function(fragment_identifications, precursor_parsed){
  
  n_IDs <- fragment_identifications[position != max(precursor_parsed$position)
    ][, .(position)] %>% unique() %>% nrow()

  total_IDs <- precursor_parsed[, .(position)
    ][position == max(position)] %>% unique() %>% .[[1]] - 1
  
  
  return(
    round(n_IDs / total_IDs * 100, 1)
  )
  
}


calculate_score <- function(fragment_identifications){
  
  average_monomer_mass <- 280
  
  n_residues <- 4
  
  n_identifications <- fragment_identifications[, .(position, parent_type)
    ] %>% unique() %>% nrow()
  
  
  ID_sequence_tags <- fragment_identifications[, terminus := ifelse(str_detect(parent_type, "[abcd]"), "L", "R")
    ][, .(terminus, position, charge)] %>% unique()
  
  ID_sequence_tags <- ID_sequence_tags[order(terminus, position, charge)
    ][, tag_ID := c(0, cumsum(diff(position) > 1)), by = terminus
    ][, 
      .(
        mean_position = mean(position),
        mean_charge = mean(charge),
        tag_length = unique(position) %>% length() #not sure if the .N sentinel will work the way I want here
      ), 
      by = .(terminus, tag_ID)
    ][, score := ((1/n_residues)*(1/average_monomer_mass*abs(mean_charge)))^tag_length %>% log() %>% -.
    ][, rand_prob := ((1/n_residues)*(1/average_monomer_mass*abs(mean_charge)))^tag_length
    ][, .(tag_score = round(sum(score) + n_identifications), P_score = prod(rand_prob))]
  
  
  return(
    list(
      `P Score` = ID_sequence_tags$P_score,
      `Tag Score` = ID_sequence_tags$tag_score
    )
  )
  
}


#Module Server==================================================================

visualization_server <- function(input, output, session, user_inputs, parsing, data_readin, identification){
  
  #Scores
  sequence_coverage <- reactive({
    
    req(identification$fragment_identifications())
    
    calculate_sequence_coverage(
      identification$fragment_identifications(),
      parsing$precursor_parsed()
    )
    
  })
  
  output$sequence_coverage_text <- renderText(
    
    paste0("Sequence Coverage: ", round(sequence_coverage(), 1), " %")
    
  )
  
  characterization_scores <- reactive({
    
    req(identification$fragment_identifications())
    
    calculate_score(
      identification$fragment_identifications()# ,
      # user_inputs$biomolecule()
    )
    
  })
  
  output$p_score_text <- renderText(
    
    paste0("P-Score: ", formatC(characterization_scores()$`P Score`), format = "e", digits = 3)
    
  )
  
  
  #Visualizations
  fragment_map <- reactive({
    
    req(identification$fragment_identifications())
    
    plot_fragment_map(
      identification$fragment_identifications(),
      parsing$precursor_parsed()# ,
      # user_inputs$biomolecule()
    )
    
  })
  
  
  output$fragment_map_plot <- renderPlot(
    fragment_map(),
    width = 650,
    height = reactive({(max(parsing$precursor_parsed()$position) %/% 10 + 1)*50+25}),
    res = 75
  )
  

  output$fragment_identifications_download <- downloadHandler(
    filename = function(){
      "Identified Fragments.csv"
    },
    content = function(file){
      write.csv(identification$fragment_identifications(), file)
    }
  )
  
  annotated_spectrum <- reactive({
    
    plot_annotated_spectrum(
      identification$fragment_identifications(),
      data_readin$spectrum()# ,
      # user_inputs$biomolecule()
    )
    
  })
  
  output$annotated_spectrum_plotly <- renderPlotly(
    annotated_spectrum()
  )
  
  
  selected_visualization <- reactive({
    
    switch(input$which_plot,
      summarized_abundance = plot_abundance_summary(
        identification$fragment_identifications(), 
        data_readin$spectrum()# ,
        # user_inputs$biomolecule()
      ),
      position_abundance = plot_abundance_by_position(
        identification$fragment_identifications(), 
        data_readin$spectrum()# ,
        # user_inputs$biomolecule()
      ),
      ppm_error = plot_ppm_error(
        identification$fragment_identifications()# ,
        # user_inputs$biomolecule()
      ),
      charge_site = plot_chargesite_localization(
        identification$fragment_identifications(),
        parsing$precursor_parsed()# , 
        # user_inputs$biomolecule()
      )
    )
    
  })
  
  output$selected_visualization_plot <- renderPlot(
    selected_visualization(),
    #width = 1000,
    height = 400,
    res = 75
  )
  
  
}


#Module UI======================================================================

visualization_UI <- function(id){
  
  ns <- NS(id)
  

  return(
    list(
      identification_tab = tagList(
        downloadButton(ns("fragment_identifications_download"), 
          "Download Fragment Identifications"
        ),
        plotOutput(ns("fragment_map_plot"), inline = T),
        textOutput(ns("sequence_coverage_text")),
        textOutput(ns("p_score_text")),
        plotlyOutput(ns("annotated_spectrum_plotly"))
      ),
      visualization_tab = tagList(
        selectInput(ns("which_plot"),
          "Visualization",
          choices = c(
            "Summarized Abundance" = "summarized_abundance",
            "Abundance by Position" = "position_abundance",
            "PPM Error" = "ppm_error",
            "Charge Site Analysis" = "charge_site"
          ),
          selected = "summarised_abundance"
        ),
        plotOutput(ns("selected_visualization_plot"))
      )
    )
  )
  
}