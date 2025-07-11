rm(list=ls())
library(SKM)
library(BGLR) #load BLR
library(caret)
library(plyr)
library(tidyr)
library(dplyr)
library(reshape2)
library(ROI)
library(ROI.plugin.glpk)
library(pracma)  # For optimization
library(quadprog)
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
# Results
results_folder <- "Optimal_Line_Selection"

# GENERAL FILE PATHS -------------------------------------------------------------
general_summary_file_path <- file.path(results_folder, "summary_ALL.csv")
general_selection_file_path <- file.path(results_folder, "selection_ALL.csv")
# GENERAL RESULTS ----------------------------------------------------------------
#general_summary_results <- data.frame()
#general_selection_results <- data.frame()
# If general file data exists ----------------------------------------------------
#if (file.exists(general_summary_file_path) && file.exists(general_selection_file_path)) {
#  general_summary_results <- read.csv(general_summary_file_path) %>% select(-X)
#  general_selection_results <- read.csv(general_selection_file_path) %>% select(-X)
#}

# Directory to save the results --------------------------------------------------
results_dir <- file.path(
  results_folder,
  dataset_name)
mkdir(results_dir)

# LOADING DATA SET ---------------------------------------------------------------
dataset_path <- file.path("Datasets", dataset_name)
# Loading data set ---------------------------------------------------------------
Pheno = data.frame()
Geno = data.frame()
Markers <- data.frame()
load(sprintf("%s.RData", dataset_path), verbose = TRUE)
ls()


columns_list <- c()
Traits_to_evaluate <- c()
if (any(colnames(Pheno) %in% "Env")) {
  columns_list <- c(paste0("T", 1:(length(colnames(Pheno))-2)))
  colnames(Pheno) <- c("Line", "Env", columns_list)
  
  Traits_to_evaluate <- setdiff(colnames(Pheno), c("Line", "Env"))
} else {
  columns_list <- c(paste0("T", 1:(length(colnames(Pheno))-1)))
  colnames(Pheno) <- c("Line", columns_list)
  
  Traits_to_evaluate <- colnames(Pheno)[-1]
}

if (nrow(Pheno) != nrow(Geno)) {
  # If there's no existing line in Geno, from Pheno
  Pheno$Line <- trimws(as.character(Pheno$Line))
  geno_names <- trimws(as.character(colnames(Geno)))
  Pheno <- Pheno[Pheno$Line %in% geno_names, ]
  # If there's a repeated line
  Pheno <- Pheno[!duplicated(Pheno$Line), ]
}

for (i in 1:length(columns_list)) {
  column_name <- columns_list[i]
  Pheno[[column_name]] <- as.numeric(Pheno[[column_name]])
}
#X=wheat.X
#head(X[,1:5])
#XS=scale(X)
#head(XS[,1:5])
#rownames(XS)=rownames(Pheno)
#Geno=XS%*%t(XS)/ncol(XS)
#head(Geno[1:5,1:5])
#Geno=Geno
#head(Pheno)
#dim(Pheno)
#dim(Geno)
######Selecting the traits to be evaluated###############
#Traits_to_evaluate

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
Pheno <- Pheno %>% arrange(Line)

# Select the unique lines both in pheno and geno and sort them
final_geno_lines <- intersect(Pheno$Line, rownames(Geno)) %>% sort()
Geno <- Geno[final_geno_lines, final_geno_lines]
dim(Geno)

##########Files to save the outputs##########
Summary_all_traits=data.frame()
Selection_all_traits=data.frame()

All_envs <- c(NA)
if (any(colnames(Pheno) %in% "Env")) {
  All_envs <- unique(Pheno$Env)
}

All_Pheno <- Pheno

for (k in 1:length(All_envs)) {
  current_env <- All_envs[k]
  print("----------------------------")
  
  if (!is.na(current_env)) {
    print(paste0("ENV = ", current_env))
    Pheno <- All_Pheno[All_Pheno$Env == current_env,]
    Pheno <- droplevels(Pheno)
  }
  
  for (t in 1:length(Traits_to_evaluate)){
    trait=Traits_to_evaluate[t]
    print(paste0("TRAIT = ", trait))
    
    Observed_trait=Pheno[, c(trait)]
    if (dataset_name == "Maize") {
      y <- Pheno[, c("Line", trait)] + 10
    } else {
      y <- Pheno[, c("Line", trait)]
    }
    dim(y)
    
    BLUE_y <- y
    
    total_elements_length <- nrow(y)-0
    
    BLUE_y <- BLUE_y[1:total_elements_length, ]
    Observed_trait=Observed_trait[1:total_elements_length]
    
    # Second solution quadratic problem - quadprog -----------------------------------
    D_matrix <- Geno
    D_matrix <- D_matrix[1:total_elements_length, 1:total_elements_length]
    
    d_vec <- BLUE_y[[2]][1:total_elements_length]
    
    # Define the variables for optimization
    n <-total_elements_length   # Number of decision variables
    x <- Variable(n,integer=T)  ##Type of decision variables
    
    for (i in 1:length(percentages)) {
      Top_percentage_to_select <- percentages[i]
      print(paste0("PERCENTAGE - ", Top_percentage_to_select))
      
      for (j in 1:length(k_list)) {
        # Define the objective function xAx
        k <- k_list[j]
        print(paste0("K VALUE - ", k))
        
        Q1 <- as.matrix(D_matrix)
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
        if (is.na(result$getValue(x)[1])) {
          Summary_Trait_t=data.frame(Data=dataset_name,Env=current_env,Trait=trait,Top_percentage_to_select=Top_percentage_to_select,K_Value=k,PM=NA,TBV_QP=NA,TBV_Trad=NA,
                                     Mean_BV_QP=NA,Mean_BV_Trad=NA,Var_Opt_QP=NA, Var_Trad=NA,Expected_Loss_QP=NA,Expected_Loss_Trad=NA,
                                     Ratio_QP=NA, Ratio_Trad=NA,
                                     Av_Rel_QP=NA,
                                     Ave_Rel_Trad=NA
          )
          Summary_all_traits=rbind(Summary_all_traits,Summary_Trait_t)
          
          Selection_Trait_t=data.frame(Data=dataset_name,Env=current_env,Trait=trait,Top_percentage_to_select=Top_percentage_to_select,K_Value=k,data.frame(Line=NA, Observed=NA, Predicted=NA))
          
          Selection_all_traits=rbind(Selection_all_traits,Selection_Trait_t) 
        } else {
          Pred_sol_QP = result$getValue(x)
          Pred_sol_QP = round(Pred_sol_QP)
          Pred_sol_QPO = Observed_trait*c(Pred_sol_QP)
          
          Data=data.frame(Line=BLUE_y$Line,Observed = Observed_trait, Predicted = Pred_sol_QPO)
          PM=best_lines_match(Data, proportion = Top_percentage_to_select) 
          PM
          
          #plot(Data$Observed, Data$Predicted)
          pos_LP_selected=which(Pred_sol_QP==1)
          LP_Selected=Data$Observed[pos_LP_selected]
          
          
          No_Top_20=length(LP_Selected)
          Sort_Obs=sort(Data$Observed,decreasing = T)
          Trad_Selected=Sort_Obs[1:No_Top_20]
          
          x_sel=round(result$getValue(x))
          pos_sel_rand=which(Data$Observed>=min(Trad_Selected))
          x_sel_rand=rep(0,length(Data$Observed))
          x_sel_rand[pos_sel_rand]=1
          
          Expected_Gain_Opt=t(x_sel)%*%d_vec
          Expected_Gain_Opt
          Expected_Gain_Rand=t(x_sel_rand)%*%d_vec
          Expected_Gain_Rand
          
          Mean_GBV_LP=Expected_Gain_Opt/No_Top_20
          Mean_GBV_Trad=Expected_Gain_Rand/No_Top_20
          
          Var_Opt= k * t(x_sel)%*%Q1%*%x_sel
          Var_Opt
          Var_Rand= k * t(x_sel_rand)%*%Q1%*%x_sel_rand
          Var_Rand
          Expected_Loss_QP=Expected_Gain_Opt-Var_Opt
          Gain_Opt=Expected_Gain_Opt/Var_Opt
          Gain_Opt
          Expected_Loss_Trad=Expected_Gain_Rand-Var_Rand
          Gain_Rand=Expected_Gain_Rand/Var_Rand
          Gain_Rand
          
          # Indices of individuals in the sample (e.g., individuals 1 to 5)
          sample_indices_LP=which(round(Pred_sol_QP)==1)
          
          # Subset the GRM
          rownames(D_matrix)=1:nrow(D_matrix)
          colnames(D_matrix)=1:nrow(D_matrix)
          G_sub_LP <- D_matrix[sample_indices_LP,sample_indices_LP]
          
          # Calculate average relatedness
          G_sub_LP=as.matrix(G_sub_LP)
          m_LP <- nrow(G_sub_LP)
          diag(G_sub_LP)=rep(0, nrow(G_sub_LP))
          
          G_sub_LP <- as.matrix(G_sub_LP)
          diag(G_sub_LP) <- rep(0, nrow(G_sub_LP))
          G_rel_LP=c(G_sub_LP)
          off_diag_sum_LP <- sum(G_rel_LP)  # Total off-diagonal sum sum
          average_relatedness_LP <- off_diag_sum_LP / (m_LP * (m_LP - 1))
          
          # Indices of individuals in the sample (e.g., individuals 1 to 5)
          sample_indices_Trad=1:No_Top_20
          
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

          Summary_Trait_t=data.frame(Data=dataset_name,Env=current_env,Trait=trait,Top_percentage_to_select=Top_percentage_to_select,K_Value=k,PM=PM,TBV_QP=Expected_Gain_Opt,TBV_Trad=Expected_Gain_Rand,
                                     Mean_BV_QP=Mean_GBV_LP,Mean_BV_Trad=Mean_GBV_Trad,Var_Opt_QP=Var_Opt, Var_Trad=Var_Rand,Expected_Loss_QP=Expected_Loss_QP,Expected_Loss_Trad=Expected_Loss_Trad,
                                     Ratio_QP=Gain_Opt, Ratio_Trad=Gain_Rand,
                                     Av_Rel_QP=average_relatedness_LP,
                                     Ave_Rel_Trad=average_relatedness_Trad
          )            
          
          Selection_Trait_t=data.frame(Data=dataset_name,Env=current_env,Trait=trait,Top_percentage_to_select=Top_percentage_to_select,K_Value=k,Data)

          Summary_all_traits=rbind(Summary_all_traits,Summary_Trait_t)
          Selection_all_traits=rbind(Selection_all_traits,Selection_Trait_t) 
        }
      }
    }
  }
}

Summary_all_traits
Selection_all_traits

Summary_all_traits

# Write each dataset results file ------------------------------------------------
write.csv(Summary_all_traits, paste(results_dir, paste0("summary_ALL_", dataset_name, ".csv"), sep = "/"))
write.csv(Selection_all_traits, paste(results_dir, paste0("selection_ALL_", dataset_name, ".csv"), sep = "/"))