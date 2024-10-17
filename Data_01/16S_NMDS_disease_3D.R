# This script runs the plot of 3D NMDS on Bray-Curtis distances of bacterial community 
# composition between healthy and diseased colonies

library(tidyverse)
library(readxl)
library(rgl)

nmds <- read_tsv(file="raw_data/melipona16s.subsample.braycurtis.3dnmds.axes")
metadata <- read_excel(path="raw_data/metadata.xlsx")
metadata_nmds <- inner_join(metadata, nmds, by=c('sample'='group')) %>%
  mutate(disease = factor(disease,
                        levels=c("FALSE", "TRUE")),
        disease_color = case_when(disease == "TRUE" ~ "red",
                                   disease == "FALSE" ~ "blue"
                                )
        )

plot3d(x=metadata_nmds$axis1, y=metadata_nmds$axis2, z=metadata_nmds$axis3,
       xlab = "NMDS Axis 1", ylab = "NMDS Axis 2", zlab = "NMDS Axis 3",
       col=metadata_nmds$disease_color, type="s", size=1)

snapshot3d("16S_NMDS_3D_disease")


