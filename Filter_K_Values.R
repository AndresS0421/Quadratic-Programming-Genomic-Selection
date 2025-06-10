library(dplyr)

# Get the results file
results_file_names <- c("Groundnut", "Indica", "Japonica", "Rice_IRRI_Philippines_Spindel_2015", "Eucalyptus_Australian_Calister_2022")
results_file_name_index <- 5

results_file_name <- results_file_names[results_file_name_index]

results_folder <- "Optimal_Line_Selection"
file_name_path <- file.path(results_folder, paste(results_file_name, paste("summary_ALL", paste0(results_file_name, ".csv"), sep = "_"), sep = "/"))

file_data <- read.csv(file_name_path)
file_data <- file_data[-1]
print(dim(file_data))

# Filter k_value = 1
file_data <- file_data %>% filter(K_Value == 1)
print(dim(file_data))

# Put abs in columns
file_data$Av_Rel_QP <- abs(file_data$Av_Rel_QP)
file_data$Ave_Rel_Trad <- abs(file_data$Ave_Rel_Trad)

# Get ratios
file_data <- file_data %>%
  mutate(
    Ratio_Mean = file_data$Mean_BV_QP / file_data$Mean_BV_Trad,
    Ratio_Var = file_data$Var_Opt_QP / file_data$Var_Trad,
    Ratio_Ave_Rel = file_data$Av_Rel_QP / file_data$Ave_Rel_Trad
  )


write.csv(file_data, file_name_path)