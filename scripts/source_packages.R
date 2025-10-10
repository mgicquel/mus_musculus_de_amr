# Source file for libraries and functions used in the project


mylibraries <- c("dplyr", "tidyr",  "data.table", "units", "stringr", "phyloseq", "abind", "here", "forcats",
                 "ggplot2", "corrplot", "viridis", "ggpubr", "ggsankey", "RColorBrewer", "officedown", "flextable", "extrafont", "cowplot",
                 "ggcorrplot", "ggrepel", "ggspatial",
                 "ggraph", "igraph", "shades",
                 "tmap", "ggmap", "terra", "tidyterra", "sf", "ggdensity", "rnaturalearth",
                 "dismo", "Hmsc", "MCMCvis", "tictoc")

for (i in 1:length(mylibraries)) {
  if(mylibraries[i] %in% rownames(installed.packages()) == FALSE) {install.packages(mylibraries[i])}
}

lapply(mylibraries, require, character.only = TRUE)

