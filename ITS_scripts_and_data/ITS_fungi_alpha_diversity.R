#This script contains the analysis of fungal alpha diversity and forager weight per month, location and health status
#Diseased colonies are those that manifested any symptom of disease during the course of our study

library(tidyverse)
library(broom)
library(readxl)

get_metadata <- function(){
  metadata <- read_excel(path="raw_data_fungi/metadata.xlsx",
                         col_types=c(location = "text", colony = "text", sister = "text", disease = "logical",
                                     forager = "text", month = "text", weight = "numeric", sample = "text")
  )
  return(metadata)
}

alpha <- read_tsv(file="raw_data_fungi/ITS_Fungi_only.agc.unique_list.groups.ave-std.summary",
                  col_types=cols(group = col_character())) %>%
  filter(method=='ave') %>%
  select(group, sobs, shannon, shannoneven, invsimpson, coverage)
metadata <- get_metadata()
meta_alpha <- inner_join(metadata, alpha, by=c('sample'='group'))

# plotting shannon index x month x colony

ggplot(meta_alpha, aes(x=month, y=shannon, color=colony)) +
  geom_jitter(shape=19, size=2, width=0.2) +
  scale_color_manual(name=NULL,
                     values=c("black", "blue", "red", "green", "yellow", "pink"),
                     breaks=c("BP1", "BP4", "BP5", "PA1", "PA4", "PA5"),
                     labels=c("BP1", "BP4", "BP5", "PA1", "PA4", "PA5")) +
  scale_x_discrete(limits=c("Jan", "Feb", "Mar", "Apr"),
                   labels=c("JAN", "FEB", "MAR", "APR")) +
  labs(title="Shannon diversity in all samples",
       x=NULL,
       y="Fungi SDI") +
  theme_classic()

#plotting shannon index x month x disease

ggplot(meta_alpha, aes(x=month, y=shannon, color=disease)) +
  geom_jitter(shape=19, size=2, width=0.2) +
  scale_color_manual(name=NULL,
                     values=c("blue", "red"),
                     breaks=c("FALSE", "TRUE"),
                     labels=c("Healty", "Diseased")) +
  scale_x_discrete(limits=c("Jan", "Feb", "Mar", "Apr"),
                   labels=c("JAN", "FEB", "MAR", "APR")) +
  labs(title="Shannon diversity in relation to disease",
       x=NULL,
       y="Fungi SDI") +
  theme_classic()

#box plot of shannon index x colony health

ggplot(meta_alpha, aes(x=disease, y=shannon, color=disease)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(shape=19)+
  scale_color_manual(name=NULL,
                     values=c("blue", "red"),
                     breaks=c("FALSE", "TRUE"),
                     labels=c("Healthy", "Diseased")) +
  scale_x_discrete(limits=c("FALSE", "TRUE"),
                   labels=c("Healthy", "Diseased")) +
  labs(title="Relationship between fungal Shannon Index and disease",
       x=NULL,
       y="Fungi SDI") +
  theme_classic()

#box plot of eveness x colony health

ggplot(meta_alpha, aes(x=disease, y=shannoneven, color=disease)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(shape=19)+
  scale_color_manual(name=NULL,
                     values=c("blue", "red"),
                     breaks=c("FALSE", "TRUE"),
                     labels=c("Healthy", "Diseased")) +
  scale_x_discrete(limits=c("FALSE", "TRUE"),
                   labels=c("Healthy", "Diseased")) +
  labs(title="Relationship between fungal eveness and disease",
       x=NULL,
       y="Fungi Eveness") +
  theme_classic()

#box plot of eveness x month

ggplot(meta_alpha, aes(x=month, y=shannoneven, color=month)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(shape=19)+
  scale_color_manual(name=NULL,
                     values=c("darkblue", "darkred", "darkgreen", "orange"),
                     breaks=c("Jan", "Feb", "Mar", "Apr"),
                     labels=c("JAN", "FEB", "MAR", "APR")) +
  scale_x_discrete(limits=c("Jan", "Feb", "Mar", "Apr"),
                   labels=c("JAN", "FEB", "MAR", "APR")) +
  labs(title="Temporal changes in fungal eveness",
       x=NULL,
       y="Fungi Eveness") +
  theme_classic()

#box plot of shannon index x month

ggplot(meta_alpha, aes(x=month, y=shannon, color=month)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(shape=19)+
  scale_color_manual(name=NULL,
                     values=c("darkblue", "darkred", "darkgreen", "orange"),
                     breaks=c("Jan", "Feb", "Mar", "Apr"),
                     labels=c("JAN", "FEB", "MAR", "APR")) +
  scale_x_discrete(limits=c("Jan", "Feb", "Mar", "Apr"),
                   labels=c("JAN", "FEB", "MAR", "APR")) +
  labs(title="Temporal changes in fungal SDI",
       x=NULL,
       y="Fungi SDI") +
  theme_classic()

#box plot of weight x month

ggplot(meta_alpha, aes(x=month, y=weight, color=month)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(shape=19)+
  scale_color_manual(name=NULL,
                     values=c("darkblue", "darkred", "darkgreen", "orange"),
                     breaks=c("Jan", "Feb", "Mar", "Apr"),
                     labels=c("JAN", "FEB", "MAR", "APR")) +
  scale_x_discrete(limits=c("Jan", "Feb", "Mar", "Apr"),
                   labels=c("JAN", "FEB", "MAR", "APR")) +
  labs(title="Temporal changes in forager weight",
       x=NULL,
       y="Forager Weight (g)") +
  theme_classic()

#box plot of weight x location

ggplot(meta_alpha, aes(x=location, y=weight, color=location)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(shape=19)+
  scale_color_manual(name=NULL,
                     values=c("darkblue", "darkred"),
                     breaks=c("BP", "PA"),
                     labels=c("BP", "PA")) +
  scale_x_discrete(limits=c("BP", "PA"),
                   labels=c("BP", "PA")) +
  labs(title="Relationship between forager weight and location",
       x=NULL,
       y="forager weight (g)") +
  theme_classic()

#box plot of weight x disease

ggplot(meta_alpha, aes(x=disease, y=weight, color=disease)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(shape=19)+
  scale_color_manual(name=NULL,
                     values=c("blue", "red"),
                     breaks=c("FALSE", "TRUE"),
                     labels=c("Healthy", "Diseased")) +
  scale_x_discrete(limits=c("FALSE", "TRUE"),
                   labels=c("Healthy", "Diseased")) +
  labs(title="Relationship between forager weight and disease",
       x=NULL,
       y="forager weight (g)") +
  theme_classic()

# statistics

meta_alpha %>%
  mutate(disease = as.character(disease)) %>%
  select(sample, shannon, disease, location, month) %>%
  pivot_longer(cols=c(disease, location, month), names_to="characteristic", values_to="value") %>%
  drop_na() %>%
  nest(data = -characteristic) %>%
  mutate(tests = map(data, ~tidy(kruskal.test(shannon ~ value, data=.x)))) %>%
  unnest(cols=tests) %>%
  select(-data) %>%
  mutate(p.value.adj = p.adjust(p.value, method="BH"))

pairwise.wilcox.test(g=meta_alpha$month, x=meta_alpha$shannon, p.adjust.method="BH")

meta_alpha %>%
  mutate(disease = as.character(disease)) %>%
  select(sample, shannoneven, disease, location, month) %>%
  pivot_longer(cols=c(disease, location, month), names_to="characteristic", values_to="value") %>%
  drop_na() %>%
  nest(data = -characteristic) %>%
  mutate(tests = map(data, ~tidy(kruskal.test(shannoneven ~ value, data=.x)))) %>%
  unnest(cols=tests) %>%
  select(-data) %>%
  mutate(p.value.adj = p.adjust(p.value, method="BH"))

pairwise.wilcox.test(g=meta_alpha$month, x=meta_alpha$shannoneven, p.adjust.method="BH")

meta_alpha %>%
  mutate(disease = as.character(disease)) %>%
  select(sample, weight, disease, location, month) %>%
  pivot_longer(cols=c(disease, location, month), names_to="characteristic", values_to="value") %>%
  drop_na() %>%
  nest(data = -characteristic) %>%
  mutate(tests = map(data, ~tidy(kruskal.test(weight ~ value, data=.x)))) %>%
  unnest(cols=tests) %>%
  select(-data) %>%
  mutate(p.value.adj = p.adjust(p.value, method="BH"))

pairwise.wilcox.test(g=meta_alpha$month, x=meta_alpha$weight, p.adjust.method="BH")

cor.test(meta_alpha$shannon, meta_alpha$weight)

cor.test(meta_alpha$shannon, meta_alpha$weight, method = "spearman")

lm_shannon_weight <- lm(shannon~weight + month, data=meta_alpha)
summary(lm_shannon_weight)

ggplot(meta_alpha, aes(x=weight, y=shannon, color=month)) +
  geom_point() +
  geom_smooth(method="lm", show.legend=FALSE) +
  scale_color_manual(name=NULL,
                     values=c("darkblue", "darkred", "darkgreen", "orange"),
                     breaks=c("Jan", "Feb", "Mar", "Apr"),
                     labels=c("JAN", "FEB", "MAR", "APR")) +
  labs(title="Linear correlation of weight and fungal SDI",
       x="Bee Weight (g)",
       y="Fungi SDI") +
  theme_classic()


