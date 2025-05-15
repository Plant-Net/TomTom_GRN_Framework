library("ggplot2")
library("gridExtra")
library("dplyr")
library("tidyr")
library("ggfortify")
library("tibble")

plot_rnaseq_distributions <- function(long_df) {
  # Check that the required columns exist
  if (!all(c("Gene", "Sample", "Expression") %in% colnames(long_df))) {
    stop("Dataframe must contain 'Gene', 'Sample', and 'Expression' columns.")
  }
  
  # Boxplot
  boxplot <- ggplot(long_df, aes(x = Sample, y = Expression, color=Sample)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate sample labels for readability
    labs(title = "Boxplot of RNAseq Expression per Sample", x = "Sample", y = "Expression")
  
  # Density plot
  density_plot <- ggplot(long_df, aes(x = Expression, color = Sample)) +
    geom_density() +
    labs(title = "Density Plot of RNAseq Expression per Sample", x = "Expression", y = "Density") +
    theme_minimal()
  
  # Arrange both plots side by side
  combined_plot <- grid.arrange(boxplot, density_plot, ncol = 2)
  
  # Return the combined plot
  return(combined_plot)
}


create_split_data <- function(data, sample_to_keep) {
  # Step 1: Extract column names from the data
  col_names <- colnames(data)
  
  # Step 2: Split column names based on the pattern and create metadata
  split_names <- str_split_fixed(col_names, "_", 5)
  metadata <- data.frame(
    species = split_names[, 1],
    type = split_names[, 2],
    exposition = split_names[, 3],
    location = split_names[, 4],
    number = split_names[, 5]
  )
  
  # Step 3: Select columns that match the 'sample_to_keep' pattern
  selected_col <- col_names[str_detect(col_names, sample_to_keep)]
  
  # Step 4: Select the corresponding data and metadata
  selected_data <- data %>%
    select(all_of(selected_col))
  
  selected_metadata <- metadata %>%
    filter(species == sample_to_keep)
  
  # Step 5: Set row names of metadata to match the selected columns
  row.names(selected_metadata) <- selected_col
  
  # Step 6: Convert the selected data to long format for visualization
  selected_data_long <- selected_data %>%
    mutate(Gene = rownames(selected_data)) %>%
    relocate(Gene, .before = 1) %>%
    pivot_longer(cols = -Gene,  # Keep 'Gene' column, pivot all others
                 names_to = "Sample", 
                 values_to = "Expression")
  
  # Step 7: Return the data, metadata, and long-format data in a list
  return(list(selected_data = selected_data,
              selected_metadata = selected_metadata,
              selected_data_long = selected_data_long))
}