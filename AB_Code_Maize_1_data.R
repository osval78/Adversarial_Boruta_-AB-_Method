rm(list=ls())
library(BGLR) #load BLR
library(dplyr)
library(caret)
library(plyr)
library(tidyr)
library(SKM)
library(Boruta)
library(reshape2)
#####Loading data set#############
#"ADVERSARIAL BORUTA METHOD"
load("Maizev_1.RData")
Markers = Markers[rownames(Markers) %in% unique(Pheno$GID),]
ls()
GIDS_Sort=sort(rownames(Geno))
Geno=Geno[GIDS_Sort,GIDS_Sort]
top_percent <-50
AUC_Threshold=0.8
Traits_to_evaluate=colnames(Pheno)[c(5:ncol(Pheno))]
Traits_to_evaluate
ToP_Selection=20
Threshold_Opt_trn=0.8
# Variables definitions ---------------------------------------------------------
#######Model to be implemented
model_name <- "BGBLUP"
Name_data_set="GDM"
Version=1

########Folds and iterations to use 
cv_folds_number <- 5
tuning_folds_number <-10
iterations_number <- 10000
burn_in <- 2500
Pheno$Line=Pheno$GID
Pheno$Env=Pheno$Loc
# Data preparation -------------------------------------------------------------
Pheno <- Pheno %>% arrange(Env,Line)
Gids_pheno=sort(unique(Pheno$GID))
#cbind(Gids_pheno,GIDS_Sort)
# Select the unique lines both in pheno and geno and sort them

#colnames(Geno)=Gids_pheno
#rownames(Geno)=Gids_pheno
final_geno_lines <- intersect(Pheno$Line, rownames(Geno)) %>% sort()
Geno1=Geno[final_geno_lines, final_geno_lines]
dim(Geno1)
Markers=Markers[final_geno_lines,]
SVD_G=svd(Geno)
U=SVD_G$u
D=diag(sqrt(SVD_G$d))
V=SVD_G$v
PC=Markers
pos_col_No_NA0=which(apply(PC,2,var)>0)
PC=PC[,pos_col_No_NA0]
#colnames(PC)=colnames(Geno)
rownames(PC)=colnames(Geno)
gc()
Pheno_E=Pheno
Predictions_All_traits=data.frame()
AUC_all_traits=data.frame()
Accuracy_all_traits=data.frame()

for (trait in Traits_to_evaluate) {
 # trait=Traits_to_evaluate[1]
  SKM::echo("*** Trait: %s ***", trait)
  Family_test=unique(Pheno_E$Family)
  BLUEs_trait=Pheno_E[, trait]

  Predictions_All_Fam=data.frame()
  AUC_all_Fam=data.frame()
  Accuracy_all_Fam=data.frame()
  n_fam = 0
for (Fam in Family_test) {
 # Fam=Family_test[1]
  n_fam = n_fam + 1
  SKM::echo("*** Fam: %s %i de %i***", Fam, n_fam, length(Family_test))
  Tst_Set=which(Pheno_E$Family==Fam)
  y_Bin=rep(0,nrow(Pheno_E))

   y_Bin[Tst_Set]=1
   Bin_name=paste("Bin",trait,Fam,sep="_")
   Pheno_E=cbind(Pheno_E,Y_Bin=y_Bin)
   #head(Pheno_E)
   colnames(Pheno_E)=c(colnames(Pheno_E)[-ncol(Pheno_E)],Bin_name[1])
   #head(Pheno_E)
   ###Cross validation -----------------------------------------------------------
      PhenoTuning <- Pheno_E
      SKM::echo("*** Iniciando Tunning ***")
      y_tuning_Bin <- PhenoTuning %>% pull(Bin_name)
      tuning_lines <- as.character(unique(PhenoTuning$Line))
      GenoTuning <- PC[tuning_lines,]
      GenoTuning=data.matrix(GenoTuning)
      GenoTuning0 <- Geno[tuning_lines,tuning_lines]
      GenoTuning0=data.matrix(GenoTuning0)
      y_tuning_Bin_Factor=as.factor(y_tuning_Bin)
      ZL=model.matrix(~0+Line,data=PhenoTuning)
      GRMO=ZL%*%GenoTuning0%*%t(ZL)
      GRM=ZL%*%GenoTuning
      Gids_pheno_DM=colnames(ZL)
     # cbind(Gids_pheno,GIDS_Sort,Gids_pheno_DM)
      
      ZE=model.matrix(~0+Env,data=PhenoTuning)[,-1]
      ZE=ZE/sqrt(ncol(ZE))
      K.E=ZE%*%t(ZE)
      K.GEO=GRMO*K.E
     # ZGE=model.matrix(~0+GRM:ZE,data=PhenoTuning)
      ETATuningO=list(Env=list(model='RKHS',K=K.E),Line=list(model='RKHS',K=GRMO),GE=list(model='RKHS',K=K.GEO))
      
      #ZGE=colnames(K.GE)
      ######Boruta Feature selection in R
      yyy <-y_tuning_Bin_Factor
      X=data.frame(cbind(y=yyy,GRM))
      
      boruta_output <- Boruta(y ~ ., data=na.omit(X), doTrace=0)  
      
      boruta_signif <- getSelectedAttributes(boruta_output, withTentative =FALSE)
      
      
      Summary_Features=attStats(boruta_output)
      Pos_IF1=which(Summary_Features$decision=="Confirmed")
      Pos_IF1
      Pos_IF2=which(Summary_Features$decision=="Tentative")
      Pos_IF2
      Pos_IF3=which(Summary_Features$decision=="Rejected")
      Pos_IF3
      Pos_IF=c(Pos_IF1,Pos_IF2,  Pos_IF3)
      length(Pos_IF1)
      length(Pos_IF2)
      Sel_IF=Summary_Features[Pos_IF,1]
      Sel_IF=  Sel_IF+abs(min(Sel_IF))+1
      summary(Sel_IF)

      names(Sel_IF)=rownames(Summary_Features[Pos_IF,])
      Sel_IF
      importance_sorted<-sort(Sel_IF, decreasing =T)
      #hist(importance_sorted)
      P=length(importance_sorted)
      Inv_IF=1/importance_sorted**1
      SKM::echo("Tuning Terminated" )
gc()
#hist(Inv_IF)
weight_Feature=Inv_IF*P/sum(Inv_IF)
sum(weight_Feature)
#hist(weight_Feature)

Matrix_weight_Feature=matrix(rep(weight_Feature,nrow(PC)),nrow=nrow(PC),ncol=P,byrow = T)

top_vars<-names(importance_sorted)
PC_New=PC[,top_vars]
PC_New_Final=PC_New*Matrix_weight_Feature
head(PC_New_Final[,1:8])
#PC_New_Final=scale(PC_New_Final)
Geno_New=PC_New_Final%*%t(PC_New_Final)/ncol(PC_New_Final)
dim(Geno_New)
Geno_New=Geno_New[final_geno_lines, final_geno_lines]
dim(Geno_New)
SKM::echo("Generando GRM2" )
GRM2=ZL%*%Geno_New%*%t(ZL)
K.GE2=GRM2*K.E
gc()
      ETATuning1=list(Env=list(model='RKHS',K=K.E),Line=list(model='RKHS',K=GRM2),GE=list(model='RKHS',K=K.GE2))

      ##########Predictions of testing with full data
      SKM::echo("*** Starting Predictions ***")
      yy=BLUEs_trait
      yy_na=yy
      yy_na[Tst_Set]=NA
      gc()
      SKM::echo("*** Starting model0 ***")
      model0 <- BGLR::BGLR(
        y = yy_na,
        ETA = ETATuningO,
        nIter = iterations_number,
        burnIn = burn_in,
        verbose = FALSE
      )
      SKM::echo("*** Model0 generated ***")
      Predicted_trait0=model0$yHat
      #Predictions$Predicted_trait0=Predicted_trait0
      #Predictions$Observed_trait=yy
      MSE0=mse(yy[Tst_Set],Predicted_trait0[Tst_Set])
      COR0=cor(yy[Tst_Set],Predicted_trait0[Tst_Set])
      NRMSE0=nrmse(yy[Tst_Set],Predicted_trait0[Tst_Set])
      gc()
      SKM::echo("*** Starting Model***")
      model <- BGLR::BGLR(
      y = yy_na,
      ETA = ETATuning1,
        nIter = iterations_number,
        burnIn = burn_in,
        verbose = FALSE
      )
      SKM::echo("*** Model generated ***")
      Predicted_trait=model$yHat
      #Predictions$Predicted_trait=Predicted_trait
      #Predictions$Observed_trait=yy
      MSE=mse(yy[Tst_Set],Predicted_trait[Tst_Set])
      COR=cor(yy[Tst_Set],Predicted_trait[Tst_Set])
      NRMSE=nrmse(yy[Tst_Set],Predicted_trait[Tst_Set])
      
 
      Accuracy_all_Fam=rbind(Accuracy_all_Fam,data.frame(Trait=trait,Family=Fam, MSE= MSE0,COR=COR0, NRMSE= NRMSE0,MSE_Opt1= MSE,COR_Opt1=COR, NRMSE_Opt1= NRMSE))
      gc()
   #   Predictions_All_Fam=rbind(Predictions_All_Fam,Predictions)
      }

#Predictions_All_traits=rbind(Predictions_All_traits,Predictions_All_Fam)
#AUC_all_traits=rbind(AUC_all_traits,AUC_all_Fam)
Accuracy_all_traits=rbind(Accuracy_all_traits,Accuracy_all_Fam)
}

#write.csv(AUC_all_traits,file="AUC_all_traits_Boruta_V17.csv")
write.csv(Accuracy_all_traits,file="Accuracy_all_traits_Boruta_V17.csv")
#write.csv(Predictions_All_traits,file="Predictions_All_traits_Boruta_V17.csv")








