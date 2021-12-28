#This script contains the analysis of bacterial quantification and forager weight per month, location and health status
#Diseased colonies are those that manifested any symptom of disease during the course of our study

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
metadata <- get_metadata()
meta_quant <- inner_join(metadata, quant, by=c('sample'='group'))

# plotting 16S quantification x month x colony

ggplot(meta_quant, aes(x=month, y=relquant, color=colony)) +
  geom_jitter(shape=19, size=2, width=0.2) +
    scale_color_manual(name=NULL,
                     values=c("black", "blue", "red", "green", "yellow", "pink"),
                     breaks=c("BP1", "BP4", "BP5", "PA1", "PA4", "PA5"),
                     labels=c("BP1", "BP4", "BP5", "PA1", "PA4", "PA5")) +
  scale_x_discrete(limits=c("Jan", "Feb", "Mar", "Apr"),
                   labels=c("JAN", "FEB", "MAR", "APR")) +
  labs(title="Bacterial quantification in all samples",
       x=NULL,
       y="Bacteria/Bee Gene Copies") +
  theme_classic()

#plotting 16S quantification x month x disease

ggplot(meta_quant, aes(x=month, y=relquant, color=disease)) +
  geom_jitter(shape=19, size=2, width=0.2) +
  scale_color_manual(name=NULL,
                     values=c("blue", "red"),
                     breaks=c("FALSE", "TRUE"),
                     labels=c("Healty", "Diseased")) +
  scale_x_discrete(limits=c("Jan", "Feb", "Mar", "Apr"),
                   labels=c("JAN", "FEB", "MAR", "APR")) +
  labs(title="Bacterial quantification in relation to disease",
       x=NULL,
       y="Bacteria/Bee Gene Copies") +
  theme_classic()

#box plot of 16S quantification x colony health

ggplot(meta_quant, aes(x=disease, y=relquant, color=disease)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(shape=19)+
  scale_color_manual(name=NULL,
                     values=c("blue", "red"),
                     breaks=c("FALSE", "TRUE"),
                     labels=c("Healthy", "Diseased")) +
  scale_x_discrete(limits=c("FALSE", "TRUE"),
                   labels=c("Healthy", "Diseased")) +
  labs(title="Relationship between bacterial quantities and disease",
       x=NULL,
       y="Bacteria/Bee Gene Copies") +
  theme_classic()

#box plot of 16S quantification x month

ggplot(meta_quant, aes(x=month, y=relquant, color=month)) +
  ylim(0,20) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(shape=19)+
  scale_color_manual(name=NULL,
                     values=c("darkblue", "darkred", "darkgreen", "orange"),
                     breaks=c("Jan", "Feb", "Mar", "Apr"),
                     labels=c("JAN", "FEB", "MAR", "APR")) +
  scale_x_discrete(limits=c("Jan", "Feb", "Mar", "Apr"),
                   labels=c("JAN", "FEB", "MAR", "APR")) +
  labs(title="Temporal variation in bacterial quantities",
       x=NULL,
       y="Bacteria/Bee Gene Copies") +
  theme_classic()

# statistics

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

cor.test(meta_quant$relquant, meta_quant$weight)

cor.test(meta_quant$relquant, meta_quant$weight, method = "spearman")

lm_relquant_weight <- lm(relquant~weight + month, data=meta_quant)
summary(lm_relquant_weight)

ggplot(meta_quant, aes(x=weight, y=relquant, color=month)) +
  geom_point() +
  geom_smooth(method="lm", show.legend=FALSE) +
  scale_color_manual(name=NULL,
                     values=c("darkblue", "darkred", "darkgreen", "orange"),
                     breaks=c("Jan", "Feb", "Mar", "Apr"),
                     labels=c("JAN", "FEB", "MAR", "APR")) +
  labs(title=NULL,
       x="Bee Weight (g)",
       y="Bacteria/Bee Gene Copies") +
  theme_classic()
