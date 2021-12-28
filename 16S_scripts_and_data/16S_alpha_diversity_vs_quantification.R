#This script contains the analysis of bacterial quantification 
#versus alpha diversity and other metadata

library(tidyverse)
library(broom)
library(readxl)

get_metadata <- function(){
  metadata <- read_excel(path="raw_data/metadata.xlsx",
                         col_types=c(location = "text", colony = "text", sister = "text", disease = "logical",
                                     forager = "text", month = "text", weight = "numeric", sample = "text")
  )
  return(metadata)
}

quant <- read_csv(file="raw_data/16S_quantification.csv",
                  col_types=cols(group = col_character())) %>%
  select(group, relquant)
alpha <- read_tsv(file="raw_data/melipona16s.phylip.an.0.03.filter.groups.ave-std.summary",
                  col_types=cols(group = col_character())) %>%
  filter(method=='ave') %>%
  select(group, sobs, shannon, shannoneven, invsimpson, coverage)
alpha_quant <- inner_join(alpha, quant, by=c('group'='group'))

metadata <- get_metadata()
meta_quant <- inner_join(metadata, quant, by=c('sample'='group'))

# statistics

cor.test(alpha_quant$relquant, alpha_quant$shannoneven)

cor.test(alpha_quant$relquant, alpha_quant$shannoneven, method = "spearman")

cor.test(meta_quant$relquant, meta_quant$weight, method = "spearman")


meta_quant %>%
  mutate(disease = as.character(disease)) %>%
  select(sample, relquant, disease, location, month) %>%
  pivot_longer(cols=c(disease, location, month), names_to="characteristic", values_to="value") %>%
  drop_na() %>%
  nest(data = -characteristic) %>%
  mutate(tests = map(data, ~tidy(kruskal.test(relquant ~ value, data=.x)))) %>%
  unnest(cols=tests) %>%
  select(-data) %>%
  mutate(p.value.adj = p.adjust(p.value, method="BH"))

pairwise.wilcox.test(g=meta_quant$month, x=meta_quant$relquant, p.adjust.method="BH")