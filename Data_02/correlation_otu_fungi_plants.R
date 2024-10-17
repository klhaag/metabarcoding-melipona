library(readxl)
library(psych)
library(corrplot)

corrdata <- read_excel(path="raw_data_fungi/combined_fungi_plants.xlsx", 
                       col_types=c(Group = "text", OtuFG01 = "numeric", OtuFG02 = "numeric", OtuFG03 = "numeric", OtuFG04 = "numeric",
                                   OtuFG05 = "numeric", OtuFG06 = "numeric", OtuFG07 = "numeric", OtuFG08 = "numeric", OtuFG09 = "numeric",
                                   OtuFG10 = "numeric", OtuFG11 = "numeric", OtuFG12 = "numeric", OtuFG13 = "numeric", OtuFG14 = "numeric",
                                   OtuFG15 = "numeric", OtuFG16 = "numeric", Otu001 = "numeric", Otu002 = "numeric", Otu003 = "numeric",
                                   Otu004 = "numeric", Otu005 = "numeric", Otu006 = "numeric", Otu007 = "numeric", Otu008 = "numeric", 
                                   Otu009 = "numeric"))
                                   
corrresults <- cor(corrdata[,unlist(lapply(corrdata, is.numeric))], method = "spearman", use = "pairwise.complete.obs")
test <- corr.test(corrdata[2:26], method = "spearman")
corrplot(corrresults, type="upper")
p_values <- test[["p"]]