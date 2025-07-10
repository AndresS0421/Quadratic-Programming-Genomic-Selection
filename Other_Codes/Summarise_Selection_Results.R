library(dplyr)
source("utils.R")

# Get the results file
results_file_names <- c("Groundnut", "Indica", "Japonica", "Rice_IRRI_Philippines_Spindel_2015", "Eucalyptus_Australian_Calister_2022")
results_folder <- "Optimal_Line_Selection"
summarised_results_folder <- "Across_Dataset"

columns_to_filter <- c("Env", "K_Value")

summarised_results_file_path <- file.path(results_folder, summarised_results_folder, paste("summary_ALL", paste0(summarised_results_folder, ".csv"), sep = "_"))
grouped_results_file_path <- file.path(results_folder, summarised_results_folder, paste("summary_GROUP", paste0(summarised_results_folder, ".csv"), sep = "_"))

# Create an empty list to store datasets
all_grouped_data <- list()
all_summary_data <- list()

mkdir(file.path(results_folder, summarised_results_folder))

for (i in 1:length(results_file_names)) {
  results_file_name <- results_file_names[i]
  
  # summary_GROUP file data
  group_file_name_path <- file.path(results_folder, paste(results_file_name, paste("summary_GROUP", paste0(results_file_name, ".csv"), sep = "_"), sep = "/"))
  
  group_file_data <- read.csv(group_file_name_path)
  group_file_data <- group_file_data[, -c(1:2)]
  
  all_grouped_data[[i]] <- group_file_data
  
  # summary_ALL file data
  all_file_name_path <- file.path(results_folder, paste(results_file_name, paste("summary_ALL", paste0(results_file_name, ".csv"), sep = "_"), sep = "/"))
  
  all_file_data <- read.csv(all_file_name_path)
  all_file_data <- all_file_data[-1]
  all_file_data <- all_file_data[, !colnames(all_file_data) %in% columns_to_filter]
  
  all_summary_data[[i]] <- all_file_data
}

# Bind all datasets by row
grouped_combined_data <- bind_rows(all_grouped_data)
all_combined_data <- bind_rows(all_summary_data)

# Calculate column-wise mean
group_summarised_by_method <- grouped_combined_data %>% 
  group_by(Top_percentage_to_select, Method) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  
# Add dataset name in Data column
all_combined_data$Data <-"Across_Dataset"

# Store the mean of the combined data into a file
write.csv(group_summarised_by_method, grouped_results_file_path)
write.csv(all_combined_data, summarised_results_file_path)