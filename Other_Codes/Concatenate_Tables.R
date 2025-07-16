#rm(list = ls(all = true))

library(dplyr)
library(ggplot2)

source("utils.r")
names <- c("Ratio_Mean_Var", "Ratio_Mean", "Ratio_Variance", "Ratio_Relatedness")

# Dataset
datasets <- c("Groundnut", "Indica", "Japonica", "Rice_IRRI_Philippines_Spindel_2015", "Eucalyptus_Australian_Calister_2022", "Across_Dataset")
# Results folder
results_folder_name <- "Optimal_Line_Selection"
# Store info list
all_datasets <- list()

# All info file path
all_datasets_file_path <- file.path(results_folder_name, "summary_ALL_datasets.csv")

for (index in 1:length(datasets)) {
  dataset_name <- datasets[index]
  results_dataset <- read.csv(paste(results_folder_name, dataset_name, paste0("summary_ALL_", dataset_name, ".csv"), sep = "/"))
  
  data_sets <- results_dataset$Data
  data_sets <- unique(data_sets)
  
  # Get repetitions number for each top percentage
  repetitions_0.1 <- sum(results_dataset$Top_percentage_to_select == 0.1)
  repetitions_0.2 <- sum(results_dataset$Top_percentage_to_select == 0.2)
  
  # Results Format To Graph --------------------------------------------------------
  results_dataset_long <- data.frame(
    Dataset = rep(results_dataset$Data, 2),
    Trait = results_dataset$Trait,
    Top_percentage_to_select = factor(rep(results_dataset$Top_percentage_to_select, 2)),
    Method = rep(c("BV_QP", "BV_Trad"), each = nrow(results_dataset)),
    Ratio_Mean = c(results_dataset$Ratio_Mean, rep(1, nrow(results_dataset))),
    Ratio_Variance = c(results_dataset$Ratio_Var, rep(1, nrow(results_dataset))),
    Ratio_Mean_Var = c(results_dataset$Ratio_QP, results_dataset$Ratio_Trad),
    Ratio_Relatedness = c(results_dataset$Ratio_Ave_Rel, rep(1, nrow(results_dataset)))
  )
  
  results_dataset_long$Top_percentage_to_select <- as.factor(results_dataset_long$Top_percentage_to_select)
  
  summary_results_dataset_long <- results_dataset_long %>%
    group_by(Dataset, Top_percentage_to_select, Method) %>%
    summarise(
      across(
        .cols = c(Ratio_Mean, Ratio_Variance, Ratio_Mean_Var, Ratio_Relatedness),
        .fns = list(avg = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE)),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
  
  # Divide sd between the repetitions number
  summary_results_dataset_long <- summary_results_dataset_long %>% 
    mutate(across(ends_with("_sd"), 
                  ~ case_when(
                    Top_percentage_to_select == 0.1 ~ .x / repetitions_0.1,
                    Top_percentage_to_select == 0.2 ~ .x / repetitions_0.2,
                    TRUE ~ .x
                  )))
  
  # Add data to previous data "all_datasets"
  all_datasets[[index]] <- summary_results_dataset_long
}

all_datasets <- bind_rows(all_datasets)

write.csv(all_datasets, all_datasets_file_path)