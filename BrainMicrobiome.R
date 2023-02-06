#set up working directory where the genus and brain data is:
setwd("")

#FORMAT DATA LIKE IN TUTORIAL

library("randomForest")
library("plyr") # for the "arrange" function
library("rfUtilities") # to test model significance
library("caret") # to get leave-one-out cross-validation accuracies and also contains the nearZeroVar function 
library("mosaic")

#reading tutorial data
#OTU
otu_table_brain <- read.table("otu_table_genus_brain.txt", sep="\t", header=T, row.names=1, stringsAsFactors=FALSE, comment.char="")  
metadata_brain <- read.table("metadata_brain.txt", sep="\t", header=T, row.names=1, stringsAsFactors=TRUE, comment.char="")

#explore data
dim(otu_table_brain) #88 samples, 1858 OTUS

#explore metadata
dim(metadata_brain) #3 metadata columns

#Insula model
otu_table_scaled_insula <- data.frame(t(otu_table_brain))  
otu_table_scaled_insula$insula <- metadata_brain[rownames(otu_table_scaled_insula), "insula"]  

set.seed(123)
RF_insula <- randomForest( x=otu_table_scaled_insula[,1:(ncol(otu_table_scaled_insula)-1)] , 
                           y=otu_table_scaled_insula$insula , ntree=501, importance=TRUE, proximities=TRUE )  
RF_insula #%Var 10.21
saveRDS( file = "RF_insula.rda" ,RF_insula )


#Error plot for # of trees
insula_tress_plot<-plot(RF_insula)
insula_tress_plot

#R squared distributions
rsquared <- data.frame(r.squared = RF_insula$rsq)
insula.r.squared.dist<-ggplot(rsquared, aes(r.squared)) + 
  geom_density() +
  theme_bw()
insula.r.squared.dist


#Permutation test (1000) (TAKES A WHILE TO RUN)
RF_insula_sig <- rf.significance( x=RF_insula ,  xdata=otu_table_scaled_insula[,1:(ncol(otu_table_scaled_insula)-1)] , nperm=1000 , ntree=501 )  
RF_insula_sig #Model is signigicant at p=0.004 (r2=0.09)
saveRDS( file = "RF_insula_sig.rda" ,RF_insula_sig )

#Accuracy Estimated by Cross-validation
fit_control <- trainControl( method = "LOOCV" )  
RF_insula_loocv <- train( otu_table_scaled_insula[,1:(ncol(otu_table_scaled_insula)-1)] , y=otu_table_scaled_insula[, ncol(otu_table_scaled_insula)] , method="rf", ntree=501 , tuneGrid=data.frame( mtry=104 ) , trControl=fit_control )
RF_insula_loocv
saveRDS( file = "RF_insula_loocv.rda" ,RF_insula_loocv )


#Identifying important features

varImpPlot_insula<-varImpPlot(RF_insula)
importance(RF_insula)

RF_insula_imp <- as.data.frame( RF_insula$importance )
RF_insula_imp$features <- rownames( RF_insula_imp )
RF_insula_imp_sorted <- arrange( RF_insula_imp  , desc(`%IncMSE`)  )
RF_insula_imp_sorted #top %IncMSE =60758.146060

barplot(RF_insula_imp_sorted$`%IncMSE`, ylab="% Increase in Mean Squared Error (Variable Importance)", main="RF Regression Variable Importance Distribution")
top_10_features_insula<-barplot(RF_insula_imp_sorted[1:10,"%IncMSE"], names.arg=RF_insula_imp_sorted[1:10,"features"] , ylab="% Increase in Mean Squared Error (Variable Importance)", las=2, ylim=c(0,30000), main="Regression RF Insula")  

#Test model without selected features (Veillonella, Enterococcus, Finegoldia, Pseudomonas, Peptacetobacter, Flavobacterium, Acinetobacter, Massilistercora, Lacticaseibacillus, Shewanella)

otu_table_scaled_insula_no_features = subset(otu_table_scaled_insula, select = -c(Veillonella, Enterococcus, Finegoldia, Pseudomonas, Peptacetobacter, Flavobacterium, Acinetobacter, Massilistercora, Lacticaseibacillus, Shewanella))
set.seed(123)
RF_insula_2 <- randomForest( x=otu_table_scaled_insula_no_features[,1:(ncol(otu_table_scaled_insula_no_features)-1)] , 
                           y=otu_table_scaled_insula_no_features$insula , ntree=501, importance=TRUE, proximities=TRUE )  
RF_insula_2 #%Var -0.04 (no longer large R2)

saveRDS( file = "RF_insula_2.rda" ,RF_insula_2 )


#Error plot for # of trees
insula_tress_plot<-plot(RF_insula_2)
insula_tress_plot

#R squared distributions
rsquared <- data.frame(r.squared = RF_insula_2$rsq)
insula.r.squared.dist<-ggplot(rsquared, aes(r.squared)) + 
  geom_density() +
  theme_bw()
insula.r.squared.dist

#Amygdala model
otu_table_scaled_amygdala <- data.frame(t(otu_table_brain))  
otu_table_scaled_amygdala$amygdala <- metadata_brain[rownames(otu_table_scaled_amygdala), "amygdala"]  

set.seed(123)
RF_amygdala <- randomForest( x=otu_table_scaled_amygdala[,1:(ncol(otu_table_scaled_amygdala)-1)] , 
                             y=otu_table_scaled_amygdala$amygdala , ntree=501, importance=TRUE, proximities=TRUE )  
RF_amygdala #%Var -3.14
saveRDS( file = "RF_amygdala.rda" ,RF_amygdala )


#Error plot for # of trees
amygdala_tress_plot<-plot(RF_amygdala)
amygdala_tress_plot

#R squared distributions
rsquared <- data.frame(r.squared = RF_amygdala$rsq)
amygdala.r.squared.dist<-ggplot(rsquared, aes(r.squared)) + 
  geom_density() +
  theme_bw()
amygdala.r.squared.dist


#Hippocampus model
otu_table_scaled_hippocampus <- data.frame(t(otu_table_brain))  
otu_table_scaled_hippocampus$hippocampus <- metadata_brain[rownames(otu_table_scaled_hippocampus), "hippocampus"]  

set.seed(123)
RF_hippocampus <- randomForest( x=otu_table_scaled_hippocampus[,1:(ncol(otu_table_scaled_hippocampus)-1)] , 
                                y=otu_table_scaled_hippocampus$hippocampus , ntree=501, importance=TRUE, proximities=TRUE )  
RF_hippocampus #%V -16.4 
saveRDS( file = "RF_hippocampus.rda" ,RF_hippocampus )


#Error plot for # of trees
hippocampus_tress_plot<-plot(RF_hippocampus)
hippocampus_tress_plot 

#R squared distributions
rsquared <- data.frame(r.squared = RF_hippocampus$rsq)
hippocampus.r.squared.dist<-ggplot(rsquared, aes(r.squared)) + 
  geom_density() +
  theme_bw()
hippocampus.r.squared.dist




