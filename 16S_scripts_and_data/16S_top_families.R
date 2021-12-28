library(tidyverse)
library(readxl)
library(broom)
library(purrr)

metadata <- read_excel(path="raw_data/metadata.xlsx", range = cell_cols(1:8),
                       col_types=c(location = "text", colony = "text", sister = "text", disease = "logical", 
                                   forager = "text", month = "text", weight = "numeric", sample = "text"))

metadata <- mutate(metadata, month = factor(month, levels=c("Apr", "Mar", "Feb", "Jan")))

taxonomy <- read_tsv(file="raw_data/melipona16s.otu.taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep=";")

otu_data <- read_tsv("raw_data/melipona16s.phylip.an.0.03.filter.0.03.subsample.shared", 
                     col_types=cols(Group=col_character(),
                                    .default = col_double())) %>%
  select(Group, starts_with("Otu")) %>%
  pivot_longer(cols=-Group, names_to="otu", values_to="count") %>%
  mutate(rel_abund=count/10000)

otu_data %>%
  group_by(Group) %>%
  summarise(N = sum(count), .groups="drop") %>%
  summarise(mean_N = mean(N), sd_N = sd(N))

agg_family_data <- inner_join(otu_data, taxonomy, by="otu") %>%
  group_by(Group, family) %>%
  summarize(agg_rel_abund = sum(rel_abund))

agg_family_metadata <- inner_join(metadata, agg_family_data, by=c("sample" = "Group"))

agg_family_metadata %>%
  group_by(family) %>%
  summarize(median=median(agg_rel_abund)) %>%
  arrange((desc(median)))

top_taxa <- agg_family_metadata %>%
  group_by(family) %>%
  summarize(median=median(agg_rel_abund)) %>%
  arrange((desc(median))) %>% 
  top_n(4, median) %>%
  pull(family)

agg_family_metadata %>%
  filter(family %in% top_taxa) %>%
  mutate(agg_rel_abund = agg_rel_abund + 1/21000) %>%
  mutate(genus=factor(family, levels=top_taxa)) %>%
  ggplot(aes(x=family, y=agg_rel_abund, color=disease, fill=disease)) +
  geom_boxplot(alpha=0.5) +
    scale_color_manual(name=NULL,
                     values=c("darkred", "darkblue"),
                     breaks=c("TRUE", "FALSE"),
                     labels=c("Diseased", "Healthy")) +
    scale_fill_manual(name=NULL,
                      values=c("darkred", "darkblue"),
                      breaks=c("TRUE", "FALSE"),
                      labels=c("Diseased", "Healthy")) +
    labs(x=NULL,
       y="Relative abundance %") +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100)) +
  coord_flip() +
  theme_classic() +
  theme()

agg_family_metadata %>%
  filter(family %in% top_taxa) %>%
  mutate(agg_rel_abund = agg_rel_abund + 1/21000) %>%
  mutate(genus=factor(family, levels=top_taxa)) %>%
  ggplot(aes(x=family, y=agg_rel_abund, color=month, fill=month)) +
  geom_boxplot(alpha=0.5) +
  scale_color_manual(name=NULL,
                     values=c("darkred", "darkblue", "darkgreen", "orange"),
                     breaks=c("Jan", "Feb", "Mar", "Apr"),
                     labels=c("JAN", "FEB", "MAR", "APR")) +
  scale_fill_manual(name=NULL,
                     values=c("darkred", "darkblue", "darkgreen", "orange"),
                     breaks=c("Jan", "Feb", "Mar", "Apr"),
                     labels=c("JAN", "FEB", "MAR", "APR")) +
  labs(x=NULL,
       y="Relative abundance %") +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100)) +
  coord_flip() +
  theme_classic() +
  theme()

agg_family_metadata %>%
  filter(family %in% top_taxa) %>%
  mutate(agg_rel_abund = agg_rel_abund + 1/21000) %>%
  mutate(genus=factor(family, levels=top_taxa)) %>%
  ggplot(aes(x=family, y=agg_rel_abund, color=location, fill=location)) +
  geom_boxplot(alpha=0.5) +
  scale_color_manual(name=NULL,
                     values=c("darkred", "darkblue"),
                     breaks=c("BP", "PA"),
                     labels=c("BP", "POA")) +
  scale_fill_manual(name=NULL,
                     values=c("darkred", "darkblue"),
                     breaks=c("BP", "PA"),
                     labels=c("BP", "POA")) +
  labs(x=NULL,
       y="Relative abundance %") +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100)) +
  coord_flip() +
  theme_classic() +
  theme()

#statistics

# this tests the significance of the differences in relative abundances of
# bacterial top families between colonies with different health status,
# correcting for multiple comparisons

family_tests1 <- agg_family_metadata %>%
  nest(sample_data = c(-family)) %>%
  mutate(test=map(sample_data, ~tidy(kruskal.test(agg_rel_abund~disease, data=.)))) %>%
  unnest(test) %>%
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
  arrange(p.value.adj)

# this tests the significance of the differences in relative abundances of
# bacterial top families between months, correcting for multiple comparisons

family_tests2 <- agg_family_metadata %>%
  nest(sample_data = c(-family)) %>%
  mutate(test=map(sample_data, ~tidy(kruskal.test(agg_rel_abund~month, data=.)))) %>%
  unnest(test) %>%
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
  arrange(p.value.adj)

# this tests the significance of the differences in relative abundances of
# bacterial top families between different locations, correcting for multiple comparisons

family_tests3 <- agg_family_metadata %>%
  nest(sample_data = c(-family)) %>%
  mutate(test=map(sample_data, ~tidy(kruskal.test(agg_rel_abund~location, data=.)))) %>%
  unnest(test) %>%
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
  arrange(p.value.adj)