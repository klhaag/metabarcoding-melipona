# This script runs the plot of 3D NMDS on Bray-Curtis distances of bacterial community 
# composition between months

library(tidyverse)
library(readxl)
library(rgl)

nmds <- read_tsv(file="raw_data_fungi/ITS_Fungi_only.agc.unique_list.0.05.filter.0.05.subsample.braycurtis.3dnmds.axes")
metadata <- read_excel(path="raw_data_fungi/metadata.xlsx")
metadata_nmds <- inner_join(metadata, nmds, by=c('sample'='group')) %>%
  mutate(month = factor(month,
                        levels=c("Jan", "Feb", "Mar", "Apr")),
         month_color = case_when(month == "Jan" ~ "blue",
                                 month == "Feb" ~ "red",
                                 month == "Mar" ~ "green",
                                 month == "Apr" ~ "yellow"
                                 
         )
  )

plot3d(x=metadata_nmds$axis1, y=metadata_nmds$axis2, z=metadata_nmds$axis3,
       xlab = "NMDS Axis 1", ylab = "NMDS Axis 2", zlab = "NMDS Axis 3",
       col=metadata_nmds$month_color, type="s", size=1)

snapshot3d("ITS_fungi_NMDS_3D_month")


