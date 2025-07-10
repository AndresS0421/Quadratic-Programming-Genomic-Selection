library(dplyr)

# Dataset
datasets <- c("Groundnut", "Indica", "Japonica", "Rice_IRRI_Philippines_Spindel_2015", "Eucalyptus_Australian_Calister_2022")
dataset_selected_index <- 5

dataset_name <- datasets[dataset_selected_index]

results_folder <- "Optimal_Line_Selection"

results_dataset <- read.csv(paste(results_folder, dataset_name, paste0("summary_ALL_", dataset_name, ".csv"), sep = "/"))

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

results_to_graph_types <- c("_By_Dataset")

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

# Write summarised data
write.csv(summary_results_dataset_long, paste(results_folder, dataset_name, paste0("summary_GROUP_", dataset_name, ".csv"), sep = "/"))