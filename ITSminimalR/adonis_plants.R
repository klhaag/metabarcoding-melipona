library(tidyverse)
library(readxl)
library(vegan)

set.seed(19962606)
permutations <- 10000

metadata <- read_excel(path="raw_data_plants/metadata.xlsx")

distance <- read_tsv("raw_data_plants/ITS_plants_only.agc.unique_list.0.05.filter.0.05.subsample.braycurtis.0.05.square.dist", 
                     skip=1, col_names=FALSE)

colnames(distance) <- c("sample", distance$X1)

meta_distance <- inner_join(metadata, distance, by="sample")

all_dist <- meta_distance %>%
  select(all_of(.[["sample"]])) %>%
  as.dist()

adonis(all_dist~location, 
       data=meta_distance, 
       permutations = permutations)

adonis(all_dist~location*disease, 
       data=meta_distance, 
       permutations = permutations)

adonis(all_dist~month, 
       data=meta_distance, 
       permutations = permutations)

adonis(all_dist~month*disease, 
       data=meta_distance, 
       permutations = permutations)

adonis(all_dist~month*disease, 
       data=meta_distance, 
       permutations = permutations)

adonis(all_dist~sister, 
       data=meta_distance, 
       permutations = permutations)

pairwise_p <- numeric()

# Jan vs Feb
Jan_Feb <- meta_distance %>%
  filter(month == "Jan" | month == "Feb")

Jan_Feb_dist <- Jan_Feb %>%
  select(all_of(.[["sample"]])) %>%
  as.dist()

Jan_Feb_test <- adonis(Jan_Feb_dist~month,
                       data=Jan_Feb,
                       permutations=permutations)

pairwise_p["Jan_Feb"] <- Jan_Feb_test[["aov.tab"]][["Pr(>F)"]][1]

# Jan vs Mar
Jan_Mar <- meta_distance %>%
  filter(month == "Jan" | month == "Mar")

Jan_Mar_dist <- Jan_Mar %>%
  select(all_of(.[["sample"]])) %>%
  as.dist()

Jan_Mar_test <- adonis(Jan_Mar_dist~month,
                       data=Jan_Mar,
                       permutations=permutations)

pairwise_p["Jan_Mar"] <- Jan_Mar_test[["aov.tab"]][["Pr(>F)"]][1]

# Jan vs Apr
Jan_Apr <- meta_distance %>%
  filter(month == "Jan" | month == "Apr")

Jan_Apr_dist <- Jan_Apr %>%
  select(all_of(.[["sample"]])) %>%
  as.dist()

Jan_Apr_test <- adonis(Jan_Apr_dist~month,
                       data=Jan_Apr,
                       permutations=permutations)

pairwise_p["Jan_Apr"] <- Jan_Apr_test[["aov.tab"]][["Pr(>F)"]][1]

# Feb vs Mar
Feb_Mar <- meta_distance %>%
  filter(month == "Feb" | month == "Mar")

Feb_Mar_dist <- Feb_Mar %>%
  select(all_of(.[["sample"]])) %>%
  as.dist()

Feb_Mar_test <- adonis(Feb_Mar_dist~month,
                       data=Feb_Mar,
                       permutations=permutations)

pairwise_p["Feb_Mar"] <- Feb_Mar_test[["aov.tab"]][["Pr(>F)"]][1]

# Feb vs Apr
Feb_Apr <- meta_distance %>%
  filter(month == "Feb" | month == "Apr")

Feb_Apr_dist <- Feb_Apr %>%
  select(all_of(.[["sample"]])) %>%
  as.dist()

Feb_Apr_test <- adonis(Feb_Apr_dist~month,
                       data=Feb_Apr,
                       permutations=permutations)

pairwise_p["Feb_Apr"] <- Feb_Apr_test[["aov.tab"]][["Pr(>F)"]][1]

# Mar vs Apr
Mar_Apr <- meta_distance %>%
  filter(month == "Mar" | month == "Apr")

Mar_Apr_dist <- Mar_Apr %>%
  select(all_of(.[["sample"]])) %>%
  as.dist()

Mar_Apr_test <- adonis(Mar_Apr_dist~month,
                       data=Mar_Apr,
                       permutations=permutations)

pairwise_p["Mar_Apr"] <- Mar_Apr_test[["aov.tab"]][["Pr(>F)"]][1]

p.adjust(pairwise_p, method="BH")