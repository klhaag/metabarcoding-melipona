library(readxl)
library(psych)
library(corrplot)

corrdata <- read_excel(path="raw_data/combined_bacteria_plants.xlsx", 
                       col_types=c(Group = "text", Otu0001 = "numeric", Otu0002 = "numeric", Otu0003 = "numeric", Otu0004 = "numeric",
                                   Otu0005 = "numeric", Otu0006 = "numeric", Otu0007 = "numeric", Otu0008 = "numeric", Otu0009 = "numeric",
                                   Otu0010 = "numeric", Otu0011 = "numeric", Otu0012 = "numeric", Otu0013 = "numeric", Otu0015 = "numeric",
                                   Otu0021 = "numeric", Otu0024 = "numeric", Otu0026 = "numeric", Otu0029 = "numeric", Otu001 = "numeric",
                                   Otu002 = "numeric", Otu003 = "numeric", Otu004 = "numeric", Otu005 = "numeric", Otu006 = "numeric", 
                                   Otu007 = "numeric", Otu008 = "numeric", Otu009 = "numeric"))

corrresults <- cor(corrdata[,unlist(lapply(corrdata, is.numeric))], method = "spearman", use = "pairwise.complete.obs")
test <- corr.test(corrdata[2:28], method = "spearman")
corrplot(corrresults, type="upper")




