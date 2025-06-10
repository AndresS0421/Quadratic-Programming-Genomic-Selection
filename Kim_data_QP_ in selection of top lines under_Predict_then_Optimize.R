rm(list=ls())
library(SKM)
library(BGLR) #load BLR
library(caret)
library(plyr)
library(tidyr)
library(dplyr)
library(reshape2)
library(CVXR)

# VARIABLES LIST -----------------------------------------------------------------
# Dataset
datasets <- c("Indica_AL", "Japonica_AL", "Groundnut_AL", "Maize_AL", "Rice_Kim_2020", "Pinus_Florida_USA_Resende_2012", "Rice_IRRI_Philippines_Spindel_2015", "Eucalyptus_Australian_Calister_2022")
dataset_selected_index <- 8
dataset_name <- datasets[dataset_selected_index]
# Top percentage
percentages <- c(0.1, 0.2)
# K list
k_list <- c(1)
# Other variables
iterations_number=500
burn_in=300
folds_num=5
testing_proportion=0.2
k=0.5
#name_dataset="Rice_Kim_2020"
Top_percentage_to_select=0.2

#####Loading data set#############
load("Rice_Kim_2020.RData")
ls()
Pheno1=Pheno
head(Pheno1)

Geno=Geno
head(Geno[1:5,1:5])
dim(Geno)
Geno[1:5,1:5]

######Selecting the traits to be evaluated###############
All_Traits=colnames(Pheno1)[c(3)]
All_Traits
All_Env=unique(Pheno1$Env)

###########Directory to save the results#####
results_dir <- file.path(
  "Optimal_Line_Selection_Predict_them_Optimize",
  name_dataset)
#####R code for computing percentage of mathching#######
####Specification of the percentage of top lines to predict######
best_lines_match <- function(Data, proportion = 0.1) {
  best_lines_number <- floor(nrow(Data) * proportion)
  
  best_lines <- Data %>%
    arrange(desc(Observed)) %>%
    slice(seq(best_lines_number)) %>%
    pull(Line) %>%
    as.character()
  
  predicted_lines <- Data %>%
    arrange(desc(Predicted)) %>%
    slice(seq(best_lines_number)) %>%
    pull(Line) %>%
    as.character()
  
  percentage_matching <- sum(predicted_lines %in% best_lines) /
    best_lines_number *
    100
  
  return(percentage_matching)
}

# Data preparation -------------------------------------------------------------
Pheno1 <- Pheno1 %>% arrange(Line)

# Select the unique lines both in pheno and geno and sort them
final_geno_lines <- intersect(Pheno1$Line, rownames(Geno)) %>% sort()
Geno <- Geno[final_geno_lines, final_geno_lines]
dim(Geno)

##########Sorting lines in Geno
geno_lines <- sort(rownames(Geno))
GenoF <- Geno[geno_lines, geno_lines]

Summary_all_Env=data.frame()
Selection_all_Env=data.frame() 
for (e in 1:length(All_Env)){
  #  e=1 
  Env_e=All_Env[e]
  Pheno=Pheno1[Pheno1$Env==Env_e,]
  Pheno=droplevels(Pheno)
  ZL=model.matrix(~0+Line,data=Pheno)
  Geno2=data.matrix(Geno)
  K_G=ZL%*%Geno2%*%t(ZL)
  
  Summary_all_Folds=data.frame()
  Selection_all_Folds=data.frame()
  #########Response variable in Obregon
  
  for(i in 1:folds_num) {
    #        i=1
    ##########Design matrix of line
    set.seed(i)
    folds_e_testing=sample(1:nrow(Pheno),nrow(Pheno)*testing_proportion)
    folds_e_testing  
    Summary_all_traits=data.frame()
    Selection_all_traits=data.frame()
    for (t in 1:length(All_Traits)){
      #        t=1  
      
      Trait=All_Traits[t]
      yy <- Pheno[,Trait]
      y_ff1 <- as.numeric(yy)  
      y_ff=y_ff1
      y_ff[folds_e_testing] <- NA
      
      #######ETA1 sin interaccion===PWI######
      ETA1=list(Line=list(model='RKHS',K=K_G))
      
      ##############Training the regression model with India data#############################
      model_ff<-BGLR::BGLR(
        y = y_ff,
        ETA = ETA1,
        response_type = "gaussian",
        nIter = iterations_number,
        burnIn = burn_in,
        verbose = FALSE
      )
      Predicted=model_ff$yHat[folds_e_testing]
      Observed=yy[folds_e_testing]
      
      
      Observed_trait=Predicted
      BLUE_y <- data.frame(Line=Pheno$Line[folds_e_testing],Trait=Predicted)
      colnames(BLUE_y)=c("Line",Trait)
      head(BLUE_y)
      
      #total_elements_length <- length(BLUE_y[[1]])
      total_elements_length <- nrow(BLUE_y)-0
      
      BLUE_y <- BLUE_y[1:total_elements_length, ]
      Observed_trait=Observed_trait[1:total_elements_length]
      
      # Second solution quadratic problem - quadprog -----------------------------------
      #Lines_testing=Pheno$Line[fold_i$testing]
      Geno_Testing=GenoF[folds_e_testing,folds_e_testing]
      #Geno_Testing2=GenoF[Lines_testing,Lines_testing]
      #cbind(Lines_testing,rownames(Geno_Testing1),rownames(Geno_Testing2))
      D_matrix1 <- Geno_Testing
      D_matrix <- D_matrix1[1:total_elements_length, 1:total_elements_length]
      
      d_vec <- BLUE_y[[2]][1:total_elements_length]
      d_vec_True=Observed[1:total_elements_length]
      
      # Define the variables for optimization
      n <-total_elements_length   # Number of decision variables
      x <- Variable(n,integer=T)  ##Type of decision variables
      # Define the objective function xAx
      
      Q1=as.matrix(D_matrix)
      objective <- Maximize(t(d_vec) %*%x - k *quad_form(x, Q1))
      constraints <- list(
        x >= 0,                    # Ensure non-negativity
        x <= 1,                    # Ensure binary variables
        sum(x)==round(total_elements_length*Top_percentage_to_select)            # At least one line must be selected or ==4
      )
      
      # Solve the problem to obtaing the optimal values
      problem <- Problem(objective, constraints)
      result <- solve(problem)
      
      # Output the solution
      Pred_sol_QP= result$getValue(x)
      Pred_sol_QP=round(Pred_sol_QP)
      Pred_sol_QPO=Observed_trait*c(Pred_sol_QP)
      
      Data=data.frame(Line=BLUE_y$Line,Observed=Observed_trait, Predicted=Pred_sol_QPO)
      Data_True=data.frame(Line=BLUE_y$Line,Observed=Observed, Predicted=Pred_sol_QPO)
      PM=best_lines_match(Data, proportion =Top_percentage_to_select) 
      PM
      PM_True=best_lines_match(Data_True, proportion =Top_percentage_to_select)
      PM_True
      set.seed(i)
      pos_sel_Total_Rand=sample(1:length(Data_True$Observed),length(Data_True$Observed)*Top_percentage_to_select)
      x_sel_Total_Rand=rep(0,length(Data$Observed))
      x_sel_Total_Rand[pos_sel_Total_Rand]=1
      x_sel_Total_Rand
      Data_total_Rand=data.frame(Line=BLUE_y$Line,Observed=Observed, Predicted=x_sel_Total_Rand)
      #plot(Data$Observed,Data$Predicted)
      PM_Total_Rand=best_lines_match(Data_total_Rand, proportion =Top_percentage_to_select)
      PM_Total_Rand
      
      pos_LP_selected=which(Pred_sol_QP==1)
      LP_Selected=Data$Observed[pos_LP_selected]
      
      No_Top_20=length(LP_Selected)
      Sort_Obs=sort(Data$Observed,decreasing = T)
      Trad_Selected=Sort_Obs[1:No_Top_20]
      
      Sort_Obs_True=sort(Data_True$Observed,decreasing = T)
      Trad_Selected_True=Sort_Obs_True[1:No_Top_20]
      
      x_sel=round(result$getValue(x))
      pos_sel_rand=which(Data$Observed>=min(Trad_Selected))
      x_sel_rand=rep(0,length(Data$Observed))
      x_sel_rand[pos_sel_rand]=1
      
      pos_sel_rand_True=which(Data_True$Observed>=min(Trad_Selected_True))
      x_sel_rand_True=rep(0,length(Data_True$Observed))
      x_sel_rand_True[pos_sel_rand_True]=1
      
      Expected_Gain_Opt=t(x_sel)%*%d_vec
      Expected_Gain_Opt
      Expected_Gain_Rand=t(x_sel_rand)%*%d_vec
      Expected_Gain_Rand
      
      Expected_Gain_Rand_True=t(x_sel_rand_True)%*%d_vec_True
      Expected_Gain_Rand_True
      
      Expected_Gain_Total_Rand=t(x_sel_Total_Rand)%*%d_vec_True
      Expected_Gain_Total_Rand
      
      Mean_GBV_LP=Expected_Gain_Opt/No_Top_20
      Mean_GBV_Trad=Expected_Gain_Rand/No_Top_20
      Mean_GBV_Trad_True=Expected_Gain_Rand_True/No_Top_20
      Mean_GBV_Total_Rand=Expected_Gain_Total_Rand/No_Top_20
      
      Var_Opt= k * t(x_sel)%*%Q1%*%x_sel
      Var_Opt
      Var_Rand= k * t(x_sel_rand)%*%Q1%*%x_sel_rand
      Var_Rand
      Var_Rand_True= k * t(x_sel_rand_True)%*%Q1%*%x_sel_rand_True
      Var_Rand_True
      Var_Total_Rand= k * t(x_sel_Total_Rand)%*%Q1%*%x_sel_Total_Rand
      Var_Total_Rand
      
      Expected_Loss_QP=Expected_Gain_Opt-Var_Opt
      Gain_Opt=Expected_Gain_Opt/Var_Opt
      Gain_Opt
      Expected_Loss_Trad=Expected_Gain_Rand-Var_Rand
      Gain_Rand=Expected_Gain_Rand/Var_Rand
      Gain_Rand
      Expected_Loss_Trad_True=Expected_Gain_Rand_True-Var_Rand_True
      Gain_Rand_True=Expected_Gain_Rand_True/Var_Rand_True
      Gain_Rand_True
      Expected_Loss_Total_Rand=Expected_Gain_Total_Rand-Var_Total_Rand
      Gain_Total_Rand=Expected_Gain_Total_Rand/Var_Total_Rand
      Gain_Total_Rand
      
      # Indices of individuals in the sample (e.g., individuals 1 to 5)
      sample_indices_LP=which(round(Pred_sol_QP)==1)
      
      # Subset the GRM
      rownames(D_matrix)=1:nrow(D_matrix)
      colnames(D_matrix)=1:nrow(D_matrix)
      G_sub_LP <- D_matrix[sample_indices_LP,sample_indices_LP]
      
      # Calculate average relatedness
      G_sub_LP=as.matrix(G_sub_LP)
      G_sub_LP
      m_LP <- nrow(G_sub_LP)
      diag(G_sub_LP)=rep(0, nrow(G_sub_LP))
      G_sub_LP
      
      G_rel_LP=c(G_sub_LP)
      off_diag_sum_LP <- sum(G_rel_LP)  # Total off-diagonal sum sum
      average_relatedness_LP <- off_diag_sum_LP / (m_LP * (m_LP - 1))
      
      # Indices of individuals in the sample (e.g., individuals 1 to 5)
      sample_indices_Trad=pos_sel_rand
      sample_indices_Trad
      pos_sel_rand
      # Subset the GRM
      rownames(D_matrix)=1:nrow(D_matrix)
      colnames(D_matrix)=1:nrow(D_matrix)
      G_sub_Trad <- D_matrix[sample_indices_Trad,sample_indices_Trad]
      
      # Calculate average relatedness
      G_sub_Trad=as.matrix(G_sub_Trad)
      m_Trad <- nrow(G_sub_Trad)
      diag(G_sub_Trad)=rep(0, nrow(G_sub_Trad))
      G_rel_Trad=G_sub_Trad
      off_diag_sum_Trad <- sum(G_rel_Trad)  # Total off-diagonal sum
      average_relatedness_Trad <- off_diag_sum_Trad / (m_Trad * (m_Trad - 1))
      
      
      # Indices of individuals in the sample (e.g., individuals 1 to 5)
      sample_indices_Trad_True=pos_sel_rand_True
      sample_indices_Trad_True
      
      G_sub_Trad_True <- D_matrix[sample_indices_Trad_True,sample_indices_Trad_True]
      
      # Calculate average relatedness
      G_sub_Trad_True=as.matrix(G_sub_Trad_True)
      m_Trad_True <- nrow(G_sub_Trad_True)
      diag(G_sub_Trad_True)=rep(0, nrow(G_sub_Trad_True))
      G_rel_Trad_True=G_sub_Trad_True
      off_diag_sum_Trad_True <- sum(G_rel_Trad_True)  # Total off-diagonal sum
      average_relatedness_Trad_True <- off_diag_sum_Trad_True / (m_Trad_True * (m_Trad_True - 1))
      
      
      sample_indices_Total_Rand=pos_sel_Total_Rand
      sample_indices_Total_Rand
      
      G_sub_Total_Rand <- D_matrix[sample_indices_Total_Rand,sample_indices_Total_Rand]
      
      # Calculate average relatedness
      G_sub_Total_Rand=as.matrix(G_sub_Total_Rand)
      m_Total_Rand <- nrow(G_sub_Total_Rand)
      diag(G_sub_Total_Rand)=rep(0, nrow(G_sub_Total_Rand))
      G_rel_Total_Rand=G_sub_Total_Rand
      off_diag_sum_Total_Rand <- sum(G_rel_Total_Rand)  # Total off-diagonal sum
      average_relatedness_Total_Rand <- off_diag_sum_Total_Rand / (m_Total_Rand * (m_Total_Rand - 1))
      
      
      Summary_Trait_t=data.frame(Data=name_dataset,Trait=Trait,testing_proportion=testing_proportion,Fold=i,k=k,Env=Env_e,Top_percentage_to_select=Top_percentage_to_select,PM=PM,PM_True=PM_True,PM_Total_Rand=PM_Total_Rand,TBV_QP=Expected_Gain_Opt,TBV_Trad=Expected_Gain_Rand,TBV_Trad_True=Expected_Gain_Rand_True,TBV_Total_Rand=Expected_Gain_Total_Rand,
                                 Mean_BV_QP=Mean_GBV_LP,Mean_BV_Trad=Mean_GBV_Trad,Mean_BV_Trad_True=Mean_GBV_Trad_True,Mean_BV_Total_Rand=Mean_GBV_Total_Rand,Var_Opt_QP=Var_Opt, Var_Trad=Var_Rand, Var_Trad_True=Var_Rand_True,Var_Total_Rand=Var_Total_Rand,Expected_Loss_QP=Expected_Loss_QP,Expected_Loss_Trad=Expected_Loss_Trad,Expected_Loss_Trad_True=Expected_Loss_Trad_True,Expected_Loss_Total_Rand=Expected_Loss_Total_Rand,
                                 Ratio_QP=Gain_Opt, Ratio_Trad=Gain_Rand,Ratio_Trad_True=Gain_Rand_True,Ratio_Total_Rand=Gain_Total_Rand,Ave_Rel_QP=average_relatedness_LP,Ave_Rel_Trad=average_relatedness_Trad,Ave_Rel_Trad_True=average_relatedness_Trad_True,Ave_Rel_Total_Rand=average_relatedness_Total_Rand)
      
      Summary_all_traits=rbind(Summary_all_traits,Summary_Trait_t)
      Summary_all_traits
      Selection_Trait_t=data.frame(Data=name_dataset,Trait=Trait,Fold=i,k=k,Env=Env_e,Top_percentage_to_select=Top_percentage_to_select,Data,Data_True,Data_total_Rand)
      
      Selection_all_traits=rbind(Selection_all_traits,Selection_Trait_t)
      Selection_all_traits
    }
    Summary_all_Folds=rbind(Summary_all_Folds,Summary_all_traits)
    Summary_all_Folds
    Selection_all_Folds=rbind(Selection_all_Folds,Selection_all_traits)
    Selection_all_Folds
  }
  Summary_all_Env=rbind(Summary_all_Env,Summary_all_Folds)
  Summary_all_Env
  Selection_all_Env=rbind(Selection_all_Env,Selection_all_Folds)
  Selection_all_Env
}
Summary_all_Env
Selection_all_Env
write.csv(Summary_all_Env,file="Summary_Predict_them_Optimize_Rice_Kim_2020.csv")
write.csv(Selection_all_Env,file="Selection_Predict_them_Optimize_Rice_Kim_2020.csv")
