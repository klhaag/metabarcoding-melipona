# This script runs the plot of NMDS on Bray-Curtis distances of fungal community 
# composition across months

library(tidyverse)
library(readxl)
library(ggtext)

source("16S_adonis.R")

nmds <- read_tsv(file="raw_data_fungi/ITS_Fungi_only.agc.unique_list.0.05.filter.0.05.subsample.braycurtis.0.05.square.nmds.axes")
metadata <- read_excel(path="raw_data_fungi/metadata.xlsx")
metadata_nmds <- inner_join(metadata, nmds, by=c('sample'='group')) %>%
  mutate(month = factor(month,
                        levels=c("Jan",
                                 "Feb",
                                 "Mar",
                                 "Apr")
  )
  )

my_legend <- tibble(x = c(-0.85, -0.7, 0.52),
                    y = c(0.75, -0.7, 0.83),
                    color = c("Jan", "Feb", "Mar", "Apr"),
                    label = c("<strong>Jan</strong>", "<strong>Feb</strong>", 
                              "<strong>Mar</strong>", "<strong>Apr</strong>"))

ggplot(metadata_nmds,
       aes(x=axis1, y=axis2, color=month, fill=month)) +
  stat_ellipse(geom="polygon",type="norm", level=0.75, alpha=0.2, show.legend=F) +
  geom_point(show.legend=TRUE) +
  coord_cartesian(xlim=c(-0.8, 0.8), ylim=c(-0.8, 0.8)) +
  labs(title=NULL,
       caption=NULL) +
  scale_color_manual(name=NULL,
                     values=c("darkblue", "darkred", "darkgreen", "orange"),
                     breaks=c("Jan", "Feb", "Mar", "Apr"),
                     labels=c("JAN", "FEB", "MAR", "APR")) +
  scale_fill_manual(name=NULL,
                    values=c("darkblue", "darkred", "darkgreen", "orange"),
                    breaks=c("Jan", "Feb", "Mar", "Apr"),
                    labels=c("JAN", "FEB", "MAR", "APR")) +
  theme_classic() +
  theme(
    legend.key.size = unit(0.25, "cm"),
    legend.position = c(1.0, 0.95),
    legend.background = element_rect(fill="white", color="black"),
    legend.margin = margin(t=-2, r=3, b=3, l=3),
    legend.text = element_markdown(),
    plot.margin = margin(l=1, r=4, unit="lines"),
    plot.title.position = "plot",
    plot.title = element_markdown(margin = margin(b=1, unit="lines")),
    axis.title.y = element_text(hjust = 1),
    axis.title.x = element_text(hjust = 0),
    plot.caption = element_text(hjust=0),
    plot.caption.position = "plot"
  )
