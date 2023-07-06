#This script below is from Brieuc et al. 2018
#See: Brieuc, M.S., Waters, C.D., Drinan, D.P. and Naish, K.A., 2018. A practical introduction to Random Forest for genetic association studies in ecology and evolution. Molecular ecology resources, 18(4), pp.755-766


######
# Introductory tutorial for conducting Random Forest on data with a continuous response variable 
# (i.e. Random Forest using regression trees)

# The code relies on the package randomForest (Liaw and Wiener 2002)
library(randomForest)

# Import the  data set, which includes annual information on salmon abundance and environmental/anthropogenic variables
# The objective is to identify which environmental/anthropogenic variables are associated with salmon abundance at Conne River

#set directory
setwd("C://Users/LehnertS/Desktop/Manuscripts/ConneDeclines/RandomForest/RScript/Update_July6_2023/")

#Read data
env<-read.table("All_data_1987_2021_ConneRiver_nonCorrelated_july6_2023.csv", sep=",", header=T)

#Manipulate/prepare data for analysis
rownames(env)<-env[,1] #Add rownames (years)

#Scale variables - start at column 5 (Smolt length)
#Note the excluded columns include Year, Small return, Large return, and Total returns (these are not predictors)
#In addition, these data were already assessed for correlation among variables - and highly correlated variabels were removed
scaleenv<-as.data.frame(apply(env[,5:ncol(env)], 2, function(x) scale(x)))
rownames(scaleenv)<-rownames(env) #add year row names

#Check dimensions of scaled dataset
dim(scaleenv) #This includes 35 rows (years) and 12 variables 

#First, explore the overall distribution of small returns
#This represents the scaled data
hist(as.data.frame(scale(env$Return_Small))[,1])   
hist(env$Return_Small) #non scaled   

###########################################################################################################################################
###########################################################################################################################################
# Now run Random Forest analysis. Since this is a continuous variable, we need to conduct a regression RF

# First, we need to optimize mtry by running different values of mtry at different values of ntree. 

# We will run mtry values - since we have 12 variables - these will be run from 4 to 12 
# We will initially run each of these mtry values at ntree=100 to 5000 (by increments of 100). 
# We are looking for a plateau where the proportion variation explained (PVE) stops increasing with larger values of ntree
# Once we reach the plateau, we will choose the mtry value that maximizes PVE.

#This optimization uses all 12 scaled variables with the response being 'small salmon returns' (scaled as well)
results_optimization <- matrix(data=NA , nrow = 0, ncol = 3)

for (i in seq(from = 100, to = 5000 , by = 100)){  # values of ntree
  print(i)
  for (j in c(4:12)){    #values of mtry 
    rf_ij <- randomForest(x = scaleenv, y = as.data.frame(scale(env$Return_Small))[,1],
                          importance=TRUE ,proximity=TRUE, ntree=i, mtry=j)
    results_optimization <- rbind(results_optimization, c(i,j,tail(rf_ij$rsq,1)))
  }
}

# Clean up the file format
results_optimization<-as.data.frame(results_optimization)
colnames(results_optimization)<-c("ntree", "mtry","PVE")

# Now plot results to see if there's a plateau
plot(results_optimization$ntree[results_optimization$mtry == 4],results_optimization$PVE[results_optimization$mtry == 4], type="l", col="black", xlab="ntree",ylab="PVE",ylim=c(0,1))
lines(results_optimization$ntree[results_optimization$mtry == 6],results_optimization$PVE[results_optimization$mtry == 6], col="blue")
lines(results_optimization$ntree[results_optimization$mtry == 8],results_optimization$PVE[results_optimization$mtry == 8], col="green")
lines(results_optimization$ntree[results_optimization$mtry == 10],results_optimization$PVE[results_optimization$mtry == 10], col="purple")
lines(results_optimization$ntree[results_optimization$mtry == 11],results_optimization$PVE[results_optimization$mtry == 11], col="orange")
lines(results_optimization$ntree[results_optimization$mtry == 12],results_optimization$PVE[results_optimization$mtry == 12], col="red")

#Generally PVE (percent variance explained) did not improve much beyond 8 mtry - thus we will use 8 for mtry in our analysis.


###########################################################################################################################################
###########################################################################################################################################
# Now begin the full Random Forest analyses

# Recall that even though we optimized mtry, we must now run a larger number of trees in order to achieve convergence of importance values between forests.
# As a starting point, we will grow 1000 trees and increase if necessary. We do not need to worry about this increase in ntree affecting our mtry optimization,
# Did not see much change in PVE across different number of trees here - so using 1000.

rf_all_1 = randomForest(x = scaleenv, y = as.data.frame(scale(env$Return_Small))[,1], importance=TRUE ,proximity=TRUE, mtry=8, ntree=1000)
save(rf_all_1,file="rf_all_1_jul6.Rdata")

rf_all_2 = randomForest(x = scaleenv, y = as.data.frame(scale(env$Return_Small))[,1], importance=TRUE ,proximity=TRUE, mtry=8, ntree=1000)
save(rf_all_2,file="rf_all_2_jul6.Rdata")

#Check correlation of variable importance values between forests 
importance_rf_all_1<-data.frame(importance(rf_all_1,type=1)) # Type 1 is mean decrease in accuracy 
colnames(importance_rf_all_1)<-c("importance")
importance_rf_all_2<-data.frame(importance(rf_all_2,type=1))
colnames(importance_rf_all_2)<-c("importance")

cor(importance_rf_all_1,importance_rf_all_2)  # A correlation of 0.99 for variable importance values between forests is extremely good, so we'll use 1000 trees for the remaining forests

rf_all_3 = randomForest(x = scaleenv, y = as.data.frame(scale(env$Return_Small))[,1], importance=TRUE ,proximity=TRUE, mtry=8, ntree=1000)
save(rf_all_3,file="rf_all_3_jul6.Rdata")

importance_rf_all_3<-data.frame(importance(rf_all_3,type=1))
colnames(importance_rf_all_3)<-c("importance")

############################################################################################################################################
############################################################################################################################################

# The predictive ability of regression trees is measured by the proportion of variation explained.  A value is calculated for each tree within a forest.
# For each forest, we will take the proportion variation explained from the last tree (after convergence)

rf_all_1_rsq <- tail(rf_all_1$rsq,1)
rf_all_2_rsq <- tail(rf_all_2$rsq,1)
rf_all_3_rsq <- tail(rf_all_3$rsq,1)

#Combine importance (mean decrease in accuracy) values of each locus across the three forests
importance_rf_all <-cbind(rownames(importance_rf_all_1),importance_rf_all_1,importance_rf_all_2, importance_rf_all_3)
colnames(importance_rf_all)<-c("Variable","Importance1","Importance2", "Importance3")

cor(importance_rf_all[,2:4]) #results are highly correlated r>0.99

# Export importance values for future reference
write.csv(importance_rf_all,file="rf_importance_values_SmallRetunrs_1987-2021_jul6.csv",row.names=FALSE)

############################################################################################################################################
############################################################################################################################################

#much of the remaining script for Brieuc et al. 2018 is on identifying groups of genetic loci that may be predictive of phenotype 
#in their study there are 1000s of loci, here we have only 12 predictors, and we did not explore this further. 

