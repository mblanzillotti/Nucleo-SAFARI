#Packages=======================================================================

require(shiny)
require(data.table)
require(magrittr)
require(stringr)
require(IsoSpecR)
require(rawrr)
require(tools)
require(vroom)
require(ggthemes)
require(grid)
require(scales)


#Modules========================================================================

#setwd("D:/UT Austin/Research/Shiny App/20240314")
source("Prerequisites and Helpers.R")
source("User Inputs Module.R")
source("Sequence Parsing Module.R")
source("Chemical Compositions Module.R")
source("Mass Generation Module.R")
source("mz Generation Module.R")
source("Data Readin Module.R")
source("Fragment Identification Module.R")
source("Visualizations Module.R")


#Nucleic Acid Testing===========================================================

date()

parsed_sequence <- parse_nucleic_input("r(5'Phos)GrCrGrGrArUrUrUrArG(Methyl)rCrUrCrArGrDrDrGrGrGrArGrArGrCrG(Dimethyl)rCrCrArGrAmCrUmGrArArYrArUrC(Methyl)rUrGrGrArGrHrUrCrC(Methyl)rUrGrUrGrTrUrCrGrA(Methyl)rUrCrCrArCrArGrArArUrUrCrGrCrArCrCrA")
# parsed_sequence <- parse_nucleic_input("rU*rC*rA*rC*rU*rU*rU*rC*rA*rU*rA*rA*rU*rG*rC*rU*rG*rG*")

precursor_composition <- generate_precursor_composition(parsed_sequence, nucleic_symbols)
precursor_monoisotopic <- calculate_precursor_monoisotopic(precursor_composition)
precursor_isotopes <- calculate_precursor_isotopes(precursor_composition)

date()

fragment_compositions <- generate_nucleic_fragment_compositions(precursor_composition)
fragment_isotopes <- calculate_fragment_isotopes(fragment_compositions)
fragment_mzs <- calculate_fragment_mzs(fragment_isotopes, -28, n_charge_states = 5)

date()

spectrum <- read_rawfile("20231020_OBE_Desalt_tRNAF_890_CID_q14_NCE25_10uM_2_XIC.raw")
#spectrum <- read_rawfile("20231020_OBE_Desalt_tRNAF_890_UVPD_20ms_10uM_2_XIC.raw")

#spectrum <- read_rawfile("05122021_SpinrazaPT_5-_1184_193_2ms_1mJ_5e5_120k_L1.raw")

date()

fragment_identifications <- identify_fragments(fragment_mzs, spectrum, 10, "nucleic_acid")

date()

validated_identifications <- validate_fragment_identifications(fragment_identifications, spectrum)

date()


plot_abundance_summary(validated_identifications, spectrum, "nucleic_acid")
plot_abundance_by_position(validated_identifications, spectrum, "nucleic_acid")
plot_ppm_error(validated_identifications, "nucleic_acid")
plot_fragment_map(validated_identifications, parsed_sequence, "nucleic_acid")
plot_chargesite_localization(validated_identifications, parsed_sequence, "nucleic_acid")
plot_annotated_spectrum(validated_identifications, spectrum, biomolecule = "nucleic_acid")

calculate_sequence_coverage(parsed_sequence, validated_identifications)
calculate_score(validated_identifications, "nucleic_acid")


#Protein Testing================================================================

date()

parsed_sequence <- parse_protein_input("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG")

precursor_composition <- generate_precursor_composition(parsed_sequence, protein_symbols)
precursor_monoisotopic <- calculate_precursor_monoisotopic(precursor_composition)
precursor_isotopes <- calculate_precursor_isotopes(precursor_composition)

date()

fragment_compositions <- generate_protein_fragment_compositions(precursor_composition)
fragment_isotopes <- calculate_fragment_isotopes(fragment_compositions)
fragment_mzs <- calculate_fragment_mzs(fragment_isotopes, 10, n_charge_states = 5)

date()

spectrum <- read_rawfile("07252022_UBQ_10+_857_193nm_2ms_1mJ_8e5_240k.raw")

date()

fragment_identifications <- identify_fragments(fragment_mzs, spectrum, 10, "protein")

date()

validated_identifications <- validate_fragment_identifications(fragment_identifications, 25, spectrum)

date()

plot_abundance_summary(validated_identifications, spectrum, "protein")
plot_abundance_by_position(validated_identifications, spectrum, "protein")
plot_ppm_error(validated_identifications, "protein")
plot_fragment_map(validated_identifications, parsed_sequence, "protein")
plot_chargesite_localization(validated_identifications, parsed_sequence, "protein")
ggplotly(plot_annotated_spectrum(validated_identifications, spectrum, biomolecule = "protein"))

calculate_sequence_coverage(validated_identifications, parsed_sequence)
calculate_score(validated_identifications, "protein")
