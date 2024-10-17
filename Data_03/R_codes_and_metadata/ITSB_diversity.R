library(tidyverse)
library(broom)
library(readxl)


get_metadata <- function(){
  metadata <- read_excel("metadata_supplement.xlsx",
                         col_types=c(sample = "text", colony = "text", period = "text", treatment = "text", pertreatment = "text")
  )
  return(metadata)
}

alpha <- read_tsv(file="ITSB.filtered.agc.groups.ave-std.summary",
                  col_types=cols(group = col_character())) %>%
  filter(method=='ave') %>%
  select(group, sobs, shannon, shannoneven, invsimpson, coverage)
metadata <- get_metadata()
meta_alpha <- inner_join(metadata, alpha, by=c('sample'='group'))

# plotting shannon index x treatment x period

ggplot(meta_alpha, aes(x=treatment, y=shannon, color=period, fill=period)) +
  geom_boxplot(alpha=0.9, outlier.shape=NA) +
  geom_jitter(shape=19)+
  scale_color_manual(name=NULL,
                     values=c("darkblue","darkred"),
                     breaks=c("A","D"),
                     labels=c("Before","After")) +
  scale_fill_manual(name=NULL,
                    values=c("lightblue","pink"),
                    breaks=c("A","D"),
                    labels=c("Before","After")) +
  scale_x_discrete(limits=c("AL","SA","CT"),
                   labels=c("AP+SA","SA","CT")) +
  labs(title=NULL,
       x=NULL,
       y="Shannon Index") +
  theme_classic()


#box plot of shannon index x period

ggplot(meta_alpha, aes(x=period, y=shannon, color=period, fill=period)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(shape=19)+
  scale_color_manual(name=NULL,
                     values=c("darkblue", "darkred"),
                     breaks=c("A", "D"),
                     labels=c("Before", "After")) +
  scale_x_discrete(limits=c("A", "D"),
                   labels=c("Before", "After")) +
  scale_fill_manual(name=NULL,
                    values=c("lightblue","pink"),
                    breaks=c("A","D"),
                    labels=c("Before","After")) +
  labs(title=NULL,
       x=NULL,
       y="Shannon Index") +
  theme_classic()


# statistics

meta_alpha %>%
  select(sample, shannon, colony, period, treatment, pertreatment) %>%
  pivot_longer(cols=c(colony, period, treatment), names_to="characteristic", values_to="value") %>%
  drop_na() %>%
  nest(data = -characteristic) %>%
  mutate(tests = map(data, ~tidy(kruskal.test(shannon ~ value, data=.x)))) %>%
  unnest(cols=tests) %>%
  select(-data) %>%
  mutate(p.value.adj = p.adjust(p.value, method="BH"))

pairwise.wilcox.test(g=meta_alpha$pertreatment, x=meta_alpha$shannon, p.adjust.method="BH")

meta_alpha %>%
  select(sample, invsimpson, colony, period, treatment, pertreatment) %>%
  pivot_longer(cols=c(colony, period, treatment), names_to="characteristic", values_to="value") %>%
  drop_na() %>%
  nest(data = -characteristic) %>%
  mutate(tests = map(data, ~tidy(kruskal.test(invsimpson ~ value, data=.x)))) %>%
  unnest(cols=tests) %>%
  select(-data) %>%
  mutate(p.value.adj = p.adjust(p.value, method="BH"))

pairwise.wilcox.test(g=meta_alpha$pertreatment, x=meta_alpha$invsimpson, p.adjust.method="BH")
