library(tidyverse)
library(readxl)
library(broom)
library(purrr)

#COMPARING PERIODS

metadata <- read_excel("metadata_supplement.xlsx", range = cell_cols(1:5),
                       col_types=c(sample = "text", colony = "text", period = "text", treatment = "text", pertreatment = "text"))

taxonomy <- read_tsv(file="ITSB.taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep=";")

otu_data <- read_tsv("ITSB.subsample.filter.shared", 
                     col_types=cols(Group=col_character(),
                                    .default = col_double())) %>%
  select(Group, starts_with("Otu")) %>%
  pivot_longer(cols=-Group, names_to="otu", values_to="count") %>%
  mutate(rel_abund=count/10000)

otu_data %>%
  group_by(Group) %>%
  summarise(N = sum(count), .groups="drop") %>%
  summarise(mean_N = mean(N), sd_N = sd(N))

agg_genus_data <- inner_join(otu_data, taxonomy, by="otu") %>%
  group_by(Group, genus) %>%
  summarize(agg_rel_abund = sum(rel_abund))

agg_genus_metadata <- inner_join(metadata, agg_genus_data, by=c("sample" = "Group"))

agg_genus_metadata %>%
  group_by(genus) %>%
  summarize(median=median(agg_rel_abund)) %>%
  arrange((desc(median)))

top_taxa <- agg_genus_metadata %>%
  group_by(genus) %>%
  summarize(median=median(agg_rel_abund)) %>%
  arrange((desc(median))) %>% 
  top_n(6, median) %>%
  pull(genus)

agg_genus_metadata %>%
  filter(genus %in% top_taxa) %>%
  mutate(agg_rel_abund = agg_rel_abund + 1/21000) %>%
  mutate(genus=factor(genus, levels=top_taxa)) %>%
  ggplot(aes(x=genus, y=agg_rel_abund, color=period, fill=period)) +
  geom_boxplot(alpha=0.9) +
  scale_color_manual(name=NULL,
                     values=c("darkblue", "darkred"),
                     breaks=c("A", "D"),
                     labels=c("Before", "After")) +
  scale_fill_manual(name=NULL,
                    values=c("lightblue", "pink"),
                    breaks=c("A", "D"),
                    labels=c("Before", "After")) +
  labs(x=NULL,
       y="Relative abundance %") +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100)) +
  coord_flip() +
  theme_classic() +
  theme()

#COMPARING PERTREATMENT

metadata <- read_excel("metadata_supplement.xlsx", range = cell_cols(1:5),
                       col_types=c(sample = "text", colony = "text", period = "text", treatment = "text", pertreatment = "text"))

taxonomy <- read_tsv(file="ITSB.taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep=";")

otu_data <- read_tsv("ITSB.subsample.filter.shared", 
                     col_types=cols(Group=col_character(),
                                    .default = col_double())) %>%
  select(Group, starts_with("Otu")) %>%
  pivot_longer(cols=-Group, names_to="otu", values_to="count") %>%
  mutate(rel_abund=count/10000)

otu_data %>%
  group_by(Group) %>%
  summarise(N = sum(count), .groups="drop") %>%
  summarise(mean_N = mean(N), sd_N = sd(N))

agg_genus_data <- inner_join(otu_data, taxonomy, by="otu") %>%
  group_by(Group, genus) %>%
  summarize(agg_rel_abund = sum(rel_abund))

agg_genus_metadata <- inner_join(metadata, agg_genus_data, by=c("sample" = "Group"))

agg_genus_metadata %>%
  group_by(genus) %>%
  summarize(median=median(agg_rel_abund)) %>%
  arrange((desc(median)))

top_taxa <- agg_genus_metadata %>%
  group_by(genus) %>%
  summarize(median=median(agg_rel_abund)) %>%
  arrange((desc(median))) %>% 
  top_n(4, median) %>%
  pull(genus)

agg_genus_metadata %>%
  filter(genus %in% top_taxa) %>%
  mutate(agg_rel_abund = agg_rel_abund + 1/21000) %>%
  mutate(genus=factor(genus, levels=top_taxa)) %>%
  mutate(pertreatment=factor(pertreatment, levels=c("DCT", "ACT", "DSA", "ASA", "DAL", "AAL"))) %>%
  ggplot(aes(x=genus, y=agg_rel_abund, color=pertreatment, fill=pertreatment)) +
  geom_boxplot(alpha=0.9, outlier.shape=NA) +
  scale_color_manual(name=NULL,
                     values=c("darkblue", "darkred", "darkgreen", "darkorange", "darkblue", "darkred"),
                     breaks=c("AAL", "DAL", "ASA", "DSA", "ACT", "DCT"),
                     labels=c("AP+SA - Before", "AP+SA - After", "SA - Before", "SA - After", "CT - Before", "CT - After")) +
  scale_fill_manual(name=NULL,
                    values=c("blue", "red", "green", "yellow", "lightblue", "pink"),
                    breaks=c("AAL", "DAL", "ASA", "DSA", "ACT", "DCT"),
                    labels=c("AP+SA - Before", "AP+SA - After", "SA - Before", "SA - After", "CT - Before", "CT - After")) +
  labs(x=NULL,
       y="Relative abundance %") +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100)) +
  coord_flip() +
  theme_classic() +
  theme()

#statistics

genus_tests1 <- agg_genus_metadata %>%
  nest(sample_data = c(-genus)) %>%
  mutate(test=map(sample_data, ~tidy(kruskal.test(agg_rel_abund~pertreatment, data=.)))) %>%
  unnest(test) %>%
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
  arrange(p.value.adj)

genus_tests2 <- agg_genus_metadata %>%
  nest(sample_data = c(-genus)) %>%
  mutate(test=map(sample_data, ~tidy(kruskal.test(agg_rel_abund~period, data=.)))) %>%
  unnest(test) %>%
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
  arrange(p.value.adj)

