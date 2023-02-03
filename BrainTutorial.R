#set up working directory where the genus and brain data is:
setwd("~/Library/CloudStorage/Dropbox/RandomForest")

#FORMAT DATA LIKE IN TUTORIAL

library("randomForest")
library("plyr") # for the "arrange" function
library("rfUtilities") # to test model significance
library("caret") # to get leave-one-out cross-validation accuracies and also contains the nearZeroVar function 
library("mosaic")

#reading tutorial data
#OTU
otu_table_brain <- read.table("otu_table_brain.txt", sep="\t", header=T, row.names=1, stringsAsFactors=FALSE, comment.char="")  
metadata_brain <- read.table("metadata_brain.txt", sep="\t", header=T, row.names=1, stringsAsFactors=TRUE, comment.char="")

#explore data
dim(otu_table_brain) #86 samples, 1858 OTUS

#explore metadata
dim(metadata_brain) #4 metadata columns

dist.amygdala<-hist(metadata_brain$amygdala)
stats.amygdala<-favstats(metadata_brain$amygdala)
dist.insula<-hist(metadata_brain$insula)
stats.insula<-favstats(metadata_brain$insula)
dist.hippocampus<-hist(metadata_brain$hippocampus)
stats.hippocampus<-favstats(metadata_brain$hippocampus)
dist.acc<-hist(metadata_brain$acc)
stats.acc<-favstats(metadata_brain$acc)


#PRE-PROCESSING

#removing rare features

otu_nonzero_counts <- apply(otu_table_brain, 1, function(y) sum(length(which(y > 0))))
otu_nonzero_counts_plot<-hist(otu_nonzero_counts, breaks=100, col="grey", main="", ylab="Number of OTUs", xlab="Number of Non-Zero Values") #removed OTUS in <=20% samples

remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )  
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}

otu_table_brain_removed <- remove_rare(table=otu_table_brain, cutoff_pro=0.10)

dim(otu_table_brain_removed) #314 <10 0s

#SAVE OBJECTS

#Insula model
otu_table_scaled_insula <- data.frame(t(otu_table_brain_removed))  
otu_table_scaled_insula$insula <- metadata_brain[rownames(otu_table_scaled_insula), "insula"]  

set.seed(123)
RF_insula <- randomForest( x=otu_table_scaled_insula[,1:(ncol(otu_table_scaled_insula)-1)] , 
                           y=otu_table_scaled_insula$insula , ntree=501, importance=TRUE, proximities=TRUE )  
RF_insula #%Var 9.25

saveRDS(file="RF_insula.rda", RF_insula)

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
RF_insula_sig #Model is signigicant at p=0.005

saveRDS(file="RF_insula_sig.rda", RF_insula_sig)

RF_insula_loocv <- train( otu_table_scaled_insula[,1:(ncol(otu_table_scaled_insula)-1)] , y=otu_table_scaled_insula[, ncol(otu_table_scaled_insula)] , method="rf", ntree=501 , tuneGrid=data.frame( mtry=104 ) , trControl=fit_control )
RF_insula_loocv

saveRDS(file="RF_insula_loocv.rda", RF_insula_loocv)

#Identifying important features

varImpPlot_insula<-varImpPlot(RF_insula)
importance(RF_insula)

RF_insula_imp <- as.data.frame( RF_insula$importance )
RF_insula_imp$features <- rownames( RF_insula_imp )
RF_insula_imp_sorted <- arrange( RF_insula_imp  , desc(`%IncMSE`)  )
RF_insula_imp_sorted #top %IncMSE =60758.146060

barplot(RF_insula_imp_sorted$`%IncMSE`, ylab="% Increase in Mean Squared Error (Variable Importance)", main="RF Regression Variable Importance Distribution")
top_10_features_insula<-barplot(RF_insula_imp_sorted[1:10,"%IncMSE"], names.arg=RF_insula_imp_sorted[1:10,"features"] , ylab="% Increase in Mean Squared Error (Variable Importance)", las=2, ylim=c(0,60758.146060), main="Regression RF Insula")  
saveRDS(file="top_10_features_insula.rda",top_10_features_insula)

#Amygdala model
otu_table_scaled_amygdala <- data.frame(t(otu_table_brain_removed))  
otu_table_scaled_amygdala$amygdala <- metadata_brain[rownames(otu_table_scaled_amygdala), "amygdala"]  

set.seed(123)
RF_amygdala <- randomForest( x=otu_table_scaled_amygdala[,1:(ncol(otu_table_scaled_amygdala)-1)] , 
                             y=otu_table_scaled_amygdala$amygdala , ntree=501, importance=TRUE, proximities=TRUE )  
RF_amygdala #%Var 0.52
saveRDS(file="RF_amygdala.rda", RF_amygdala)

#Error plot for # of trees
amygdala_tress_plot<-plot(RF_amygdala)
amygdala_tress_plot

#R squared distributions
rsquared <- data.frame(r.squared = RF_amygdala$rsq)
amygdala.r.squared.dist<-ggplot(rsquared, aes(r.squared)) + 
  geom_density() +
  theme_bw()
amygdala.r.squared.dist

#Permutation test (1000) (TAKES A WHILE TO RUN)
RF_amygdala_sig <- rf.significance( x=RF_amygdala ,  xdata=otu_table_scaled_amygdala[,1:(ncol(otu_table_scaled_amygdala)-1)] , nperm=1000 , ntree=501 )  
RF_amygdala_sig #Model is NOT significant at p=0.074 
saveRDS(file="RF_amygdala_sig.rda", RF_amygdala_sig)

#Accuracy Estimated by Cross-validation
fit_control <- trainControl( method = "LOOCV" )    
RF_amygdala_loocv <- train( otu_table_scaled_amygdala[,1:(ncol(otu_table_scaled_amygdala)-1)] , y=otu_table_scaled_amygdala[, ncol(otu_table_scaled_amygdala)] , method="rf", ntree=501 , tuneGrid=data.frame( mtry=104 ) , trControl=fit_control )
RF_amygdala_loocv
saveRDS(file="RF_amygdala_loocv.rda", RF_amygdala_loocv)

#Identifying important features

varImpPlot_amygdala<-varImpPlot(RF_amygdala)
importance(RF_amygdala)

RF_amygdala_imp <- as.data.frame( RF_amygdala$importance )
RF_amygdala_imp$features <- rownames( RF_amygdala_imp )
RF_amygdala_imp_sorted <- arrange( RF_amygdala_imp  , desc(`%IncMSE`)  )
RF_amygdala_imp_sorted #TOP %IncMSE 177.79103235
barplot(RF_amygdala_imp_sorted$`%IncMSE`, ylab="% Increase in Mean Squared Error (Variable Importance)", main="RF Regression Variable Importance Distribution")
top_10_features_amygdala<-barplot(RF_amygdala_imp_sorted[1:10,"%IncMSE"], names.arg=RF_amygdala_imp_sorted[1:10,"features"] , ylab="% Increase in Mean Squared Error (Variable Importance)", las=2, ylim=c(0,177.79103235), main="Regression RF Amygdala")  

#Hippocampus model
otu_table_scaled_hippocampus <- data.frame(t(otu_table_brain_removed))  
otu_table_scaled_hippocampus$hippocampus <- metadata_brain[rownames(otu_table_scaled_hippocampus), "hippocampus"]  

set.seed(123)
RF_hippocampus <- randomForest( x=otu_table_scaled_hippocampus[,1:(ncol(otu_table_scaled_hippocampus)-1)] , 
                                y=otu_table_scaled_hippocampus$hippocampus , ntree=501, importance=TRUE, proximities=TRUE )  
RF_hippocampus #%V -12.11 
saveRDS(file="RF_hippocampus.rda",RF_hippocampus)

#Error plot for # of trees
hippocampus_tress_plot<-plot(RF_hippocampus)
hippocampus_tress_plot 

#R squared distributions
rsquared <- data.frame(r.squared = RF_hippocampus$rsq)
hippocampus.r.squared.dist<-ggplot(rsquared, aes(r.squared)) + 
  geom_density() +
  theme_bw()
hippocampus.r.squared.dist


#Permutation test (1000) (TAKES A WHILE TO RUN)
RF_hippocampus_sig <- rf.significance( x=RF_hippocampus ,  xdata=otu_table_scaled_hippocampus[,1:(ncol(otu_table_scaled_hippocampus)-1)] , nperm=1000 , ntree=501 )  
RF_hippocampus_sig #not signifiant at p = 0.539 
saveRDS(file="RF_hippocampus_sig.rda", RF_hippocampus_sig)

#ACC model

otu_table_scaled_acc <- data.frame(t(otu_table_brain_removed))  
otu_table_scaled_acc$acc <- metadata_brain[rownames(otu_table_scaled_acc), "acc"]  

set.seed(123)
RF_acc <- randomForest( x=otu_table_scaled_acc[,1:(ncol(otu_table_scaled_acc)-1)] , 
                        y=otu_table_scaled_acc$acc , ntree=501, importance=TRUE, proximities=TRUE )  
RF_acc #%Var -4.05
saveRDS(file="RF_acc.rda", RF_acc)
#Error plot for # of trees
acc_tress_plot<-plot(RF_acc)
acc_tress_plot

#R squared distributions
rsquared <- data.frame(r.squared = RF_acc$rsq)
acc.r.squared.dist<-ggplot(rsquared, aes(r.squared)) + 
  geom_density() +
  theme_bw()
acc.r.squared.dist #R2 dont go >0

#Permutation test (1000) (TAKES A WHILE TO RUN)
RF_acc_sig <- rf.significance( x=RF_acc ,  xdata=otu_table_scaled_acc[,1:(ncol(otu_table_scaled_acc)-1)] , nperm=1000 , ntree=501 )  
RF_acc_sig #Model not signifiant at p = 0.149
saveRDS(file="RF_acc_sig.rda",RF_acc_sig)



