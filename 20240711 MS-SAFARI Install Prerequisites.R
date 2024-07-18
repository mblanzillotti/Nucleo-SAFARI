#Package Requirements for Nucleo-SAFARI=========================================

#This script installs prerequisite packages for use with Nucleo-SAFARI
#by M Lanzillotti as of 08 May 2024. This script assumes
#R (4.3.2 or later) and R Studio (2023.06.0 Build 421) are installed, alongside
#.NET framework (4.8.1 or later) for use with Thermo Rawfile Reader. RTools will
#be installed alongside the following packages:

# .NET Framework: https://dotnet.microsoft.com/en-us/download
# R and RStudio: https://posit.co/download/rstudio-desktop/

#Packages=======================================================================

#CRAN packages

install.packages("shiny")
install.packages("data.table")
install.packages("magrittr")
install.packages("stringr")
install.packages("IsoSpecR")
install.packages("vroom")
install.packages("nnls")
install.packages("ggthemes")
install.packages("scales")
install.packages("ggplot2")
install.packages("plotly")


#Bioconductor packages
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
  
BiocManager::install("rawrr", force = T)
rawrr::installRawFileReaderDLLs()
rawrr::installRawrrExe()
  #Requires consent to RawFileReader license agreement via Thermo


