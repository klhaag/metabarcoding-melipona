library(tidyverse)
library(readxl)

rarefy <- read_tsv(file="raw_data_plants/ITS_plants.rarefaction") %>%
  select(-contains("lci-"), -contains("hci-")) %>%
  pivot_longer(cols=c(-numsampled), names_to='sample', values_to='sobs') %>%
  mutate(sample=str_replace_all(sample, pattern="0.05-", replacement="")) %>%
  drop_na() 

metadata <- read_excel(path="raw_data_plants/metadata.xlsx", col_names=TRUE)

metadata_rarefy <- inner_join(metadata, rarefy, c("sample"))

ggplot(metadata_rarefy, aes(x=numsampled, y=sobs, color=colony, group=sample)) +
  geom_line() +
  scale_color_manual(name=NULL,
                     values=c("darkblue", "blue", "red", "green", "yellow", "orange"),
                     breaks=c("BP1", "BP4", "BP5", "PA1", "PA4", "PA5"),
                     labels=c("BP1", "BP4", "BP5", "PA1", "PA4", "PA5")) +
  labs(x="Number of Sequences Sampled per Bee",
       y="Number of OTUs per Bee") +
  theme_classic()