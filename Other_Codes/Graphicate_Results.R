#rm(list = ls(all = true))

library(dplyr)
library(ggplot2)

source("utils.r")
names <- c("Ratio_Mean_Var", "Ratio_Mean", "Ratio_Variance", "Ratio_Relatedness")

# Dataset
datasets <- c("Groundnut", "Indica", "Japonica", "Rice_IRRI_Philippines_Spindel_2015", "Eucalyptus_Australian_Calister_2022", "Across_Dataset")
dataset_selected_index <- 6
dataset_name <- datasets[dataset_selected_index]
# Results folder
results_folder_name <- "Optimal_Line_Selection"

results_dataset <- read.csv(paste(results_folder_name, dataset_name, paste0("summary_ALL_", dataset_name, ".csv"), sep = "/"))

head(results_dataset)
dim(results_dataset)

data_sets <- results_dataset$Data
data_sets <- unique(data_sets)


# delete all the previous plots
#unlink(plots_dir, recursive = true)
# create the directory where the plots are goind to be stored in
plots_dir<- paste("Results_Graphics",sep="_")
dir.create(plots_dir)

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

# Divide sd between the repetitions number
summary_results_dataset_long <- summary_results_dataset_long %>% 
  mutate(across(ends_with("_sd"), 
                ~ case_when(
                  Top_percentage_to_select == 0.1 ~ .x / repetitions_0.1,
                  Top_percentage_to_select == 0.2 ~ .x / repetitions_0.2,
                  TRUE ~ .x
                )))

# Renombrar columnas *_avg -> original, mantener *_sd
names(summary_results_dataset_long) <- gsub("_avg$", "", names(summary_results_dataset_long))

# Graphs -------------------------------------------------------------------------
for (data_set in data_sets) {
  data_set = data_set
  results_to_graph <- list(droplevels(summary_results_dataset_long[summary_results_dataset_long$Dataset == data_set, ]))
  
  cat(data_set, "\n")
  
  ### Create Dir -----------------------------------------------------------------
  dir.create(file.path(plots_dir, data_set))  
  ### Graph ALL Data -------------------------------------------------------------
  results_long_colnames <- colnames(summary_results_dataset_long)
  for (col_data_index in seq(4, length(results_long_colnames), by = 2)) {
    col_data <- results_long_colnames[col_data_index]
    col_data_sd <- results_long_colnames[col_data_index + 1]
    cat("Column: ", col_data, "\n")
    
    for (results_index in 1:length(results_to_graph)) {
      results <- results_to_graph[[results_index]]
      results_type <- results_to_graph_types[results_index]
      
      plot <- ggplot(
        results,
        aes(x = Method, y = !!sym(col_data), fill = factor(Method))
      ) +
      geom_bar(stat = "identity", position = "dodge") +  # Cambiado a "dodge" para dividir las barras
      facet_wrap(~Top_percentage_to_select) +
      scale_x_discrete(
        labels = c("0.1" = "0.1", "0.2" = "0.2")
      ) +
        geom_errorbar(
          aes(
            ymin = !!sym(col_data) - !!sym(col_data_sd), 
            ymax = !!sym(col_data) + !!sym(col_data_sd)
          ),
          width = 0.2,
          position = position_dodge(0.9),
          colour = "black"
        )+
      labs(fill = "Proportions")  # Cambia la etiqueta de la leyenda
      plot <- white_theme(plot)
      Plot <- plot + theme(
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = "bold")
      ) +
      theme_bw() +
      theme(text = element_text(size = 20))
      Plot <- vertical_x(Plot, angle = 90)
      
      ### Create Graph File ----------------------------------------------------------
      save_plot(Plot, file = file.path(plots_dir, data_set, paste0(data_set, paste0("_", col_data), results_type, ".png")))
    }
  }
}