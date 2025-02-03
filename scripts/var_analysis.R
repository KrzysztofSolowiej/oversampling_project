library(tidyverse)
library(readxl)
library(tidymodels)
library(ranger)
library(themis)
library(missMDA)
library(pbapply)
library(iml)
library(caret)
library(doParallel)
library(ggplot2)

# Data

path <- "data/original_matrix.xlsx"

gc <- read_excel(path, sheet = 1) %>% 
  mutate(Group = as.factor(Group)) %>% 
  filter(!is.na(Group))

biocrates <- read_excel(path, sheet = 2) %>%
  mutate(Group = as.factor(Group)) %>% 
  filter(!is.na(Group))

common_samples <- biocrates[["Sample Name"]][!grepl(pattern = "^QC ", x = biocrates[["Sample Name"]])]

get_group_id <- function(df, n, group_name) {
  if(n == "full") {
    df[df[["Group"]] == group_name, ]
  } else {
    df[sample(which(df[["Group"]] == group_name), n, replace = FALSE), ]
  }
}

gc_part_dat <- filter(gc, `Sample Name` %in% common_samples) %>%
  filter(!Group %in% c("G2", "G2 (enepdymoma)", "G1", "G3")) %>%
  mutate(Group = if_else(grepl("^G", as.character(Group)), "G4", as.character(Group)))

gc_G4_org <- gc_part_dat %>% filter(Group == "G4")
print("GC G4 dims:")
print(dim(gc_G4_org))

gc_MT_org <- gc_part_dat %>% filter(Group == "MT")
print("GC MT dims:")
print(dim(gc_MT_org))

bio_part_dat <- filter(biocrates, `Sample Name` %in% common_samples) %>%
  filter(!Group %in% c("G2", "G2 (enepdymoma)", "G1", "G3")) %>%
  mutate(Group = if_else(grepl("^G", as.character(Group)), "G4", as.character(Group)))

bio_G4_org <- bio_part_dat %>% filter(Group == "G4")
print("Bio G4 dims:")
print(dim(bio_G4_org))

bio_MT_org <- bio_part_dat %>% filter(Group == "MT")
print("Bio MT dims:")
print(dim(bio_MT_org))

# Functions

perform_smote <- function(data, group_col, original_size) {
  numeric_data <- data %>%
    select(where(is.numeric), !!sym(group_col))
  
  smote_data <- smote(numeric_data, var = group_col, over_ratio = original_size / nrow(data))
  
  smote_data <- as_tibble(smote_data)
  smote_data$`Sample Name` <- paste("Sample", seq(1, nrow(smote_data)))
  smote_data
}

calculate_variances <- function(data, group_col = "Group") {
  data %>%
    select(where(is.numeric)) %>%
    summarise(across(everything(), ~ var(.x, na.rm = TRUE))) %>%
    pivot_longer(cols = everything(), names_to = "Feature", values_to = "Variance")
}

# Loops

set.seed(184502)
n_sampled_values <- c(10, 20, 30)
variance_diffs <- list()
num_iterations <- 100

for (n_sampled in n_sampled_values) {
  for (i in 1:num_iterations) {
    reduced_dat_G4 <- get_group_id(gc_part_dat, n_sampled, "G4") %>%
      droplevels() %>% filter(Group == "G4")
    ov_G4 <- perform_smote(reduced_dat_G4, "Group", nrow(reduced_dat_G4) * (6.5 / (n_sampled / 10)))
    results_full <- calculate_variances(gc_G4_org)
    results_ov <- calculate_variances(ov_G4)
    
    results_full <- results_full %>%
      mutate(Dataset = "Full", n_sampled = n_sampled)
    
    results_ov <- results_ov %>%
      mutate(Dataset = "Oversampled", n_sampled = n_sampled)
    
    variance_comparison <- bind_rows(results_full, results_ov)
    variance_comparison_side_by_side <- variance_comparison %>%
      pivot_wider(names_from = Dataset, values_from = Variance)
    
    variance_diff <- variance_comparison_side_by_side %>%
      mutate(Difference = Oversampled - Full,
             LogDifference = log1p(abs(Difference)))
    
    variance_diffs[[paste(n_sampled, i)]] <- variance_diff
  }
}

variance_all_iterations <- bind_rows(variance_diffs, .id = "Iteration")

variance_summary <- variance_all_iterations %>%
  group_by(Feature, n_sampled) %>%
  summarise(MeanDifference = mean(Difference, na.rm = TRUE),
            MeanLogDifference = mean(LogDifference, na.rm = TRUE))

variance_summary$Feature <- str_sub(variance_summary$Feature, 1, 15)
variance_summary$n_sampled <- factor(variance_summary$n_sampled, levels = c(30, 20, 10))
variance_summary$Feature <- factor(variance_summary$Feature, levels = rev(sort(unique(variance_summary$Feature))))

ggplot(variance_summary, aes(x = Feature, y = MeanLogDifference, fill = factor(n_sampled))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  coord_flip() +
  labs(title = "Log-Transformed Mean Difference in Variance by Sample Size\nGC G4 data - 100 Iterations",
       x = "Metabolite",
       y = "Mean Log-Transformed Variance Difference (Oversampled - Full)") +
  scale_fill_manual(
    values = c("10" = "#ff5a5a", "20" = "#65e133", "30" = "#e798ff"), 
    name = "Sample Size (n)",
    breaks = c(10, 20, 30)
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# GC MT

variance_diffs <- list()

for (n_sampled in n_sampled_values) {
  for (i in 1:num_iterations) {
    reduced_dat_MT <- get_group_id(gc_part_dat, n_sampled, "MT") %>%
      droplevels() %>% filter(Group == "MT")

    ov_MT <- perform_smote(reduced_dat_MT, "Group", nrow(reduced_dat_MT) * (6.8 / (n_sampled / 10)))

    results_full <- calculate_variances(gc_MT_org)
    results_ov <- calculate_variances(ov_MT)
    
    results_full <- results_full %>%
      mutate(Dataset = "Full", n_sampled = n_sampled)
    
    results_ov <- results_ov %>%
      mutate(Dataset = "Oversampled", n_sampled = n_sampled)
    
    variance_comparison <- bind_rows(results_full, results_ov)
    variance_comparison_side_by_side <- variance_comparison %>%
      pivot_wider(names_from = Dataset, values_from = Variance)
    
    variance_diff <- variance_comparison_side_by_side %>%
      mutate(Difference = `Oversampled` - `Full`,
             LogDifference = log1p(abs(Difference)))
    
    variance_diffs[[paste(n_sampled, i)]] <- variance_diff
  }
}

variance_all_iterations <- bind_rows(variance_diffs, .id = "Iteration")

variance_summary <- variance_all_iterations %>%
  group_by(Feature, n_sampled) %>%
  summarise(MeanDifference = mean(Difference, na.rm = TRUE),
            MeanLogDifference = mean(LogDifference, na.rm = TRUE))

variance_summary$Feature <- str_sub(variance_summary$Feature, 1, 15)
variance_summary$n_sampled <- factor(variance_summary$n_sampled, levels = c(30, 20, 10))
variance_summary$Feature <- factor(variance_summary$Feature, levels = rev(sort(unique(variance_summary$Feature))))

ggplot(variance_summary, aes(x = Feature, y = MeanLogDifference, fill = factor(n_sampled))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  coord_flip() +
  labs(title = "Log-Transformed Mean Difference in Variance by Sample Size\nGC MT data - 100 Iterations",
       x = "Metabolite",
       y = "Mean Log-Transformed Variance Difference (Oversampled - Full)") +
  scale_fill_manual(
    values = c("10" = "#ff5a5a", "20" = "#65e133", "30" = "#e798ff"), 
    name = "Sample Size (n)",
    breaks = c(10, 20, 30)
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Biocrates G4

variance_diffs <- list()

for (n_sampled in n_sampled_values) {
  for (i in 1:num_iterations) {
    reduced_dat_G4 <- get_group_id(bio_part_dat, n_sampled, "G4") %>%
      droplevels() %>% filter(Group == "G4")
    
    ov_G4 <- perform_smote(reduced_dat_G4, "Group", nrow(reduced_dat_G4) * (6.5 / (n_sampled / 10)))
    
    results_full <- calculate_variances(bio_G4_org)
    results_ov <- calculate_variances(ov_G4)
    
    results_full <- results_full %>%
      mutate(Dataset = "Full", n_sampled = n_sampled)
    
    results_ov <- results_ov %>%
      mutate(Dataset = "Oversampled", n_sampled = n_sampled)
    
    variance_comparison <- bind_rows(results_full, results_ov)
    variance_comparison_side_by_side <- variance_comparison %>%
      pivot_wider(names_from = Dataset, values_from = Variance)
    
    variance_diff <- variance_comparison_side_by_side %>%
      mutate(Difference = `Oversampled` - `Full`,
             LogDifference = log1p(abs(Difference)))
    
    variance_diffs[[paste(n_sampled, i)]] <- variance_diff
  }
}

variance_all_iterations <- bind_rows(variance_diffs, .id = "Iteration")

variance_summary <- variance_all_iterations %>%
  group_by(Feature, n_sampled) %>%
  summarise(MeanDifference = mean(Difference, na.rm = TRUE),
            MeanLogDifference = mean(LogDifference, na.rm = TRUE))


variance_summary$Feature <- str_sub(variance_summary$Feature, 1, 15)
variance_summary$n_sampled <- factor(variance_summary$n_sampled, levels = c(30, 20, 10))
variance_summary$Feature <- factor(variance_summary$Feature, levels = rev(sort(unique(variance_summary$Feature))))

ggplot(variance_summary, aes(x = Feature, y = MeanLogDifference, fill = factor(n_sampled))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  coord_flip() +
  labs(title = "Log-Transformed Mean Difference in Variance by Sample Size\nBiocrates G4 data - 100 Iterations",
       x = "Metabolite",
       y = "Mean Log-Transformed Variance Difference (Oversampled - Full)") +
  scale_fill_manual(
    values = c("10" = "#ff5a5a", "20" = "#65e133", "30" = "#e798ff"), 
    name = "Sample Size (n)",
    breaks = c(10, 20, 30)
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Biocrates MT


variance_diffs <- list()

for (n_sampled in n_sampled_values) {
  for (i in 1:num_iterations) {
    reduced_dat_MT <- get_group_id(bio_part_dat, n_sampled, "MT") %>%
      droplevels() %>% filter(Group == "MT")

    ov_MT <- perform_smote(reduced_dat_MT, "Group", nrow(reduced_dat_MT) * (6.8 / (n_sampled / 10)))
    
    results_full <- calculate_variances(bio_MT_org)
    results_ov <- calculate_variances(ov_MT)
    
    results_full <- results_full %>%
      mutate(Dataset = "Full", n_sampled = n_sampled)
    
    results_ov <- results_ov %>%
      mutate(Dataset = "Oversampled", n_sampled = n_sampled)
    
    variance_comparison <- bind_rows(results_full, results_ov)
    variance_comparison_side_by_side <- variance_comparison %>%
      pivot_wider(names_from = Dataset, values_from = Variance)
    
    variance_diff <- variance_comparison_side_by_side %>%
      mutate(Difference = `Oversampled` - `Full`,
             LogDifference = log1p(abs(Difference)))
    
    variance_diffs[[paste(n_sampled, i)]] <- variance_diff
  }
}

variance_all_iterations <- bind_rows(variance_diffs, .id = "Iteration")

variance_summary <- variance_all_iterations %>%
  group_by(Feature, n_sampled) %>%
  summarise(MeanDifference = mean(Difference, na.rm = TRUE),
            MeanLogDifference = mean(LogDifference, na.rm = TRUE))


variance_summary$Feature <- str_sub(variance_summary$Feature, 1, 15)
variance_summary$n_sampled <- factor(variance_summary$n_sampled, levels = c(30, 20, 10))
variance_summary$Feature <- factor(variance_summary$Feature, levels = rev(sort(unique(variance_summary$Feature))))

ggplot(variance_summary, aes(x = Feature, y = MeanLogDifference, fill = factor(n_sampled))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  coord_flip() +
  labs(title = "Log-Transformed Mean Difference in Variance by Sample Size\nBiocrates MT data - 100 Iterations",
       x = "Metabolite",
       y = "Mean Log-Transformed Variance Difference (Oversampled - Full)") +
  scale_fill_manual(
    values = c("10" = "#ff5a5a", "20" = "#65e133", "30" = "#e798ff"), 
    name = "Sample Size (n)",
    breaks = c(10, 20, 30)
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
