#Script for r unning random forest analyses with final dataset

#Set directory and open libraries
setwd("C://Users/LehnertS/Desktop/Manuscripts/ConneDeclines/RandomForest/RScript/Update_July6_2023/")
library(randomForest)
library(ggplot2)
library(corrplot)
library(rfPermute)

#Read data
env<-read.table("All_data_1987_2021_ConneRiver_nonCorrelated_july6_2023.csv", sep=",", header=T)

#Manipulate/prepare data for analysis
rownames(env)<-env[,1] #Add rownames (years)


#Scale environmental data - start at column 5 
scaleenv<-as.data.frame(apply(env[,5:ncol(env)], 2, function(x) scale(x)))
rownames(scaleenv)<-rownames(env) #add year row names
dim(scaleenv)


# Evaluate correlations between variables - could consider removing highly correlated variables (>0.7)
corrplot(cor(scaleenv), method = "ellipse",
         title = "method = 'ellipse'",
         type="lower",
         tl.col = "black", # Labels color
         mar = c(2, 1, 3, 1)) 

# all correlations are |r|<0.75 in this final dataset


#Create datasets to save results (Note change nrows for different datasets)
MDA<-as.data.frame(matrix(nrow=ncol(scaleenv), ncol=1)) #Number of rows refers to number of enviro data
colnames(MDA)<-"abundance"

#tested different values of mtry and ntree to improve PVE - 
#did not improve much beyond 8 variables, and all trees produced similar PVE - used 1000
#See other script for selection info

#Run Random forest
set.seed(4563)
rf<-randomForest(scaleenv, as.data.frame(scale(env$Return_Small))[,1],
                 ntree=1000, mtry=8, importance=TRUE)

print(rf)

MDA[,1]<-as.data.frame(importance(rf, type = 1))

#Add rownames for rf$importance
rownames(MDA)<-rownames(rf$importance)

#Write table with results
write.table(MDA, file="meanMDA_SmallReturns_July6_2023.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

##Open data and make plot 
MDA_data=read.table("meanMDA_SmallReturns_July6_2023.txt", header=T)
MDA_data$Measure<- rownames(MDA_data) 
MDA_data$MDA<- MDA_data$abundance

#get correlations among variables with abundance
data<- data.frame("variable"= colnames(scaleenv), "cor" = rep(NA, 12))
for(i in 1:ncol(scaleenv)){
  
  lm1 <-cor(as.data.frame(scale(env$Return_Small))[,1],scaleenv[,i])
  data[i,2] <- lm1 
}
MDA_dat2<- cbind(MDA_data,data )
MDA_dat2$col <- rep("negative")
MDA_dat2$col[which(MDA_dat2$cor >= 0)] <- "positive"

#Plot RF results for small returns
ggplot(data=MDA_dat2, aes(x=reorder(Measure, MDA), y=MDA, fill=col))+
  scale_fill_manual(values=c("firebrick","dodgerblue4"), name="Association with\nabundance (r)")+
  geom_bar(stat="identity")+theme_bw()+
  coord_flip()+
  geom_hline(yintercept = 0)+
  scale_y_continuous(lim=c(-5,30),expand = c(0, 0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #note significance of factors is determined below - added astrerisk here for significant ones
  geom_text(aes(x=12, y=27), label="*", size=8)+
  geom_text(aes(x=11, y=18), label="*", size=8)+
  geom_text(aes(x=10, y=8), label="*", size=8)+
  ylab("Mean decrease in acccuracy")+ xlab("Variable") 

#save for manuscript

#Determine significance of variables

set.seed(443)
rf_sig<-rfPermute(scaleenv, as.data.frame(scale(env$Return_Small))[,1],
                  ntree=1000, mtry=8, importance=TRUE, num.rep = 1000)

predict(rf_sig)
summary(rf_sig)
plotImportance(rf_sig, scale = TRUE)
importance(rf_sig, scale = TRUE)

####################Run same model wtih total retunrs

#Create datasets to save results (Note change nrows for different datasets)
MDA_rf_All<-as.data.frame(matrix(nrow=ncol(scaleenv), ncol=1)) #Number of rows refers to number of enviro data
colnames(MDA_rf_All)<-"TotalAbundance"

#Run RF
set.seed(3232)
rf_All<-randomForest(scaleenv, as.data.frame(scale(env$Return_total))[,1],
                     ntree=1000, mtry=8, importance=TRUE)
print(rf_All)
MDA_rf_All[,1]<-as.data.frame(importance(rf_All, type = 1))


#Add rownames for rf_All$importance
rownames(MDA_rf_All)<-rownames(rf_All$importance)


#Write table with results
write.table(MDA_rf_All, file="meanMDA_TotalReturns_July6_2023_forcomparison.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

##Open data and make plot 
MDA_data=read.table("meanMDA_TotalReturns_July6_2023_forcomparison.txt", header=T)
MDA_data$Measure<- rownames(MDA_data) 
MDA_data$MDA<- MDA_data$TotalAbundance

#get positive and negative associations
data<- data.frame("variable"= colnames(scaleenv), "cor" = rep(NA, 12))
for(i in 1:ncol(scaleenv)){
  
  lm1 <-cor(as.data.frame(scale(env$Return_total))[,1],scaleenv[,i])
  data[i,2] <- lm1 
}
MDA_dat2<- cbind(MDA_data,data )
MDA_dat2$col <- rep("negative")
MDA_dat2$col[which(MDA_dat2$cor >= 0)] <- "positive"

#Plot results for Total returns RF
ggplot(data=MDA_dat2, aes(x=reorder(Measure, MDA), y=MDA, fill=col))+
  scale_fill_manual(values=c("firebrick","dodgerblue4"), name="Association with\ntotal returns (r)")+
  geom_bar(stat="identity")+theme_bw()+
  coord_flip()+
  geom_hline(yintercept = 0)+
  scale_y_continuous(lim=c(-5,30),expand = c(0, 0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Mean decrease in acccuracy")+ xlab("Variable")
#check significance of variables below - seems results are consistent with just using small salmon abundance - as in same 
# variables rank highly and same variables are significant)

plot(scale(env$Return_total)[,1],scaleenv$CommercialfisherySFA11)

#
set.seed(4389)
test<-rfPermute(scaleenv, as.data.frame(scale(env$Return_total))[,1],
                ntree=1000, mtry=8, importance=TRUE, num.rep = 1000)
predict(test)
summary(test)
plotImportance(test, scale = TRUE)
importance(test, scale = TRUE)
