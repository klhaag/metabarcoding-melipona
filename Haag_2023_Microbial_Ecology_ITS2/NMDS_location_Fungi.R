# This script runs the plot of NMDS on Bray-Curtis distances of fungal community 
# composition between locations

library(tidyverse)
library(readxl)
library(ggtext)

nmds <- read_tsv(file="raw_data_fungi/ITS_Fungi_only.agc.unique_list.0.05.filter.0.05.subsample.braycurtis.0.05.square.nmds.axes")
metadata <- read_excel(path="raw_data_fungi/metadata.xlsx")
metadata_nmds <- inner_join(metadata, nmds, by=c('sample'='group')) %>%
  mutate(location = factor(location,
                        levels=c("BP", "PA")
  ))

ggplot(metadata_nmds,
       aes(x=axis1, y=axis2, color=location, fill=location)) +
  stat_ellipse(geom="polygon",type="norm", level=0.75, alpha=0.2, show.legend=F) +
  geom_point(show.legend=TRUE) +
  # geom_richtext(data=my_legend,
  #           aes(x=x, y=y, label=label, color=color),
  #           inherit.aes=FALSE, show.legend = FALSE, hjust=0, fill=NA, label.color=NA) +
  coord_cartesian(xlim=c(-0.8, 0.8), ylim=c(-0.8, 0.8)) +
  labs(title=NULL,
       caption=NULL) +
  scale_color_manual(name=NULL,
                     values=c("darkred", "darkblue"),
                     breaks=c("BP", "PA"),
                     labels=c("BP", "POA")) +
  scale_fill_manual(name=NULL,
                    values=c("darkred", "darkblue"),
                    breaks=c("BP", "PA"),
                    labels=c("BP", "POA")) +
  theme_classic() +
  theme(
    legend.key.size = unit(0.25, "cm"),
    legend.position = c(1.0, 0.95),
    legend.background = element_rect(fill="white", color="black"),
    legend.text = element_markdown(),
    plot.margin = margin(l=1, r=4, unit="lines"),
    plot.title.position = "plot",
    plot.title = element_markdown(margin = margin(b=1, unit="lines")),
    axis.title.y = element_text(hjust = 1),
    axis.title.x = element_text(hjust = 0),
    plot.caption = element_text(hjust=0),
    plot.caption.position = "plot"
  )