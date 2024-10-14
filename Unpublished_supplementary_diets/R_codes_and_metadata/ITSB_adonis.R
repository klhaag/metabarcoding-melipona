library(tidyverse)
library(readxl)
library(vegan)

set.seed(19962606)
permutations <- 10000

metadata <- read_excel("metadata_supplement.xlsx")

distance <- read_tsv("ITSB.braycurtis.square.dist", 
                     skip=1, col_names=FALSE)

colnames(distance) <- c("sample", distance$X1)

meta_distance <- inner_join(metadata, distance, by="sample")

all_dist <- meta_distance %>%
  select(all_of(.[["sample"]])) %>%
  as.dist()

all_test <- adonis2(all_dist~period*treatment, 
                    data=meta_distance, 
                    permutations = permutations)

# ASA vs DSA
ASA_DSA <- meta_distance %>%
  filter(pertreatment == "ASA" | pertreatment == "DSA")

ASA_DSA_dist <- ASA_DSA %>%
  select(all_of(.[["sample"]])) %>%
  as.dist()

ASA_DSA_test <- adonis2(ASA_DSA_dist~pertreatment,
                        data=ASA_DSA,
                        permutations=permutations)

# DAL vs DCT
DAL_DCT <- meta_distance %>%
  filter(pertreatment == "DAL" | pertreatment == "DCT")

DAL_DCT_dist <- DAL_DCT %>%
  select(all_of(.[["sample"]])) %>%
  as.dist()

DAL_DCT_test <- adonis2(DAL_DCT_dist~pertreatment,
                        data=DAL_DCT,
                        permutations=permutations)

# DAL vs DSA
DAL_DSA <- meta_distance %>%
  filter(pertreatment == "DAL" | pertreatment == "DSA")

DAL_DSA_dist <- DAL_DSA %>%
  select(all_of(.[["sample"]])) %>%
  as.dist()

DAL_DSA_test <- adonis2(DAL_DSA_dist~pertreatment,
                        data=DAL_DSA,
                        permutations=permutations)

# DSA vs DCT
DSA_DCT <- meta_distance %>%
  filter(pertreatment == "DAL" | pertreatment == "DCT")

DSA_DCT_dist <- DSA_DCT %>%
  select(all_of(.[["sample"]])) %>%
  as.dist()

DSA_DCT_test <- adonis2(DSA_DCT_dist~pertreatment,
                        data=DSA_DCT,
                        permutations=permutations)



# AAL vs ACT
AAL_ACT <- meta_distance %>%
  filter(pertreatment == "AAL" | pertreatment == "ACT")

AAL_ACT_dist <- AAL_ACT %>%
  select(all_of(.[["sample"]])) %>%
  as.dist()

AAL_ACT_test <- adonis2(AAL_ACT_dist~pertreatment,
                        data=AAL_ACT,
                        permutations=permutations)

# ASA vs ACT
ASA_ACT <- meta_distance %>%
  filter(pertreatment == "ASA" | pertreatment == "ACT")

ASA_ACT_dist <- ASA_ACT %>%
  select(all_of(.[["sample"]])) %>%
  as.dist()

ASA_ACT_test <- adonis2(ASA_ACT_dist~pertreatment,
                        data=ASA_ACT,
                        permutations=permutations)

# AAL vs ASA
AAL_ASA <- meta_distance %>%
  filter(pertreatment == "AAL" | pertreatment == "ASA")

AAL_ASA_dist <- AAL_ASA %>%
  select(all_of(.[["sample"]])) %>%
  as.dist()

AAL_ASA_test <- adonis2(AAL_ASA_dist~pertreatment,
                        data=AAL_ASA,
                        permutations=permutations)

# AAL vs DSA
AAL_DSA <- meta_distance %>%
  filter(pertreatment == "AAL" | pertreatment == "DSA")

AAL_DSA_dist <- AAL_DSA %>%
  select(all_of(.[["sample"]])) %>%
  as.dist()

AAL_DSA_test <- adonis2(AAL_DSA_dist~pertreatment,
                        data=AAL_DSA,
                        permutations=permutations)

# AAL vs DCT
AAL_DCT <- meta_distance %>%
  filter(pertreatment == "AAL" | pertreatment == "DCT")

AAL_DCT_dist <- AAL_DCT %>%
  select(all_of(.[["sample"]])) %>%
  as.dist()

AAL_DCT_test <- adonis2(AAL_DCT_dist~pertreatment,
                        data=AAL_DCT,
                        permutations=permutations)

# AAL vs DAL
AAL_DAL <- meta_distance %>%
  filter(pertreatment == "AAL" | pertreatment == "DAL")

AAL_DAL_dist <- AAL_DAL %>%
  select(all_of(.[["sample"]])) %>%
  as.dist()

AAL_DAL_test <- adonis2(AAL_DAL_dist~pertreatment,
                        data=AAL_DAL,
                        permutations=permutations)

# ASA vs DAL
ASA_DAL <- meta_distance %>%
  filter(pertreatment == "ASA" | pertreatment == "DAL")

ASA_DAL_dist <- ASA_DAL %>%
  select(all_of(.[["sample"]])) %>%
  as.dist()

ASA_DAL_test <- adonis2(ASA_DAL_dist~pertreatment,
                        data=ASA_DAL,
                        permutations=permutations)

# ASA vs DCT
ASA_DCT <- meta_distance %>%
  filter(pertreatment == "ASA" | pertreatment == "DCT")

ASA_DCT_dist <- ASA_DCT %>%
  select(all_of(.[["sample"]])) %>%
  as.dist()

ASA_DCT_test <- adonis2(ASA_DCT_dist~pertreatment,
                        data=ASA_DCT,
                        permutations=permutations)

# ACT vs DAL
ACT_DAL <- meta_distance %>%
  filter(pertreatment == "ACT" | pertreatment == "DAL")

ACT_DAL_dist <- ACT_DAL %>%
  select(all_of(.[["sample"]])) %>%
  as.dist()

ACT_DAL_test <- adonis2(ACT_DAL_dist~pertreatment,
                        data=ACT_DAL,
                        permutations=permutations)

# ACT vs DSA
ACT_DSA <- meta_distance %>%
  filter(pertreatment == "ACT" | pertreatment == "DSA")

ACT_DSA_dist <- ACT_DSA %>%
  select(all_of(.[["sample"]])) %>%
  as.dist()

ACT_DSA_test <- adonis2(ACT_DSA_dist~pertreatment,
                        data=ACT_DSA,
                        permutations=permutations)

# ACT vs DCT
ACT_DCT <- meta_distance %>%
  filter(pertreatment == "ACT" | pertreatment == "DCT")

ACT_DCT_dist <- ACT_DCT %>%
  select(all_of(.[["sample"]])) %>%
  as.dist()

ACT_DCT_test <- adonis2(ACT_DCT_dist~pertreatment,
                        data=ACT_DCT,
                        permutations=permutations)




