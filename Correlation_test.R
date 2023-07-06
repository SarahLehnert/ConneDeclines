#examining correlation among variables

setwd("C://Users/LehnertS/Desktop/Manuscripts/ConneDeclines/RandomForest/RScript/Update_July6_2023/")
library(randomForest)
library(ggplot2)
library(corrplot)
library(rfPermute)

#Read data
env<-read.table("All_data_1987_2021_ConneRiver.csv", sep=",", header=T)

#Manipulate/prepare data for analysis
rownames(env)<-env[,1] #Add rownames (years)


#Scale environmental data - start at column 5 
scaleenv<-as.data.frame(apply(env[,5:ncol(env)], 2, function(x) scale(x)))
rownames(scaleenv)<-rownames(env) #add year row names
dim(scaleenv)


# Evaluate correlations between variables - could consider removing highly correlated variables (>0.7)
corrplot(cor(scaleenv, use = "complete.obs"), method = "ellipse",
         title = "method = 'ellipse'",
         type="lower",
         tl.col = "black", # Labels color
         mar = c(2, 1, 3, 1)) 

#Look at values of correlations
write.table(as.data.frame(cor(scaleenv, use = "complete.obs")), file = "Correlation_matrix1987_2021.txt", quote = F, row.names = T, col.names = T, sep="\t")

#Things to remove - 
#Use a cut-off of 0.75 here

#Note Climate index is highly correlated wtih  Sea Ice, CIL area, St27 Summer CIL, St27 Temp  - Keep Climate Index only
#Smolt weight and length are highly correalted - Keep length only which should be less prone to error than weight
#Smolt weight and length were estimated for 2021 based on data from 2016-2022 (average)
#Spring and Summer SST on South coast highly correlated - just keep Spring (may be more relavant to smolts)
#Salinity was removed as it was missing data for 2020 and is likely less relevant for salmon  

#Remove these variables and created a new dataframe
