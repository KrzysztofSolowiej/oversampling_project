if (!requireNamespace("recipes", quietly = TRUE)) {
  install.packages("recipes")
}

library(readxl)
library(tidyverse)
library(recipes)
library(themis)


path <- "data/original_matrix.xlsx"

gc <- read_excel(path, sheet = 1) %>% 
  mutate(Group = as.factor(Group)) %>% 
  filter(!is.na(Group))

biocrates <- read_excel(path, sheet = 2) %>% 
  mutate(Group = as.factor(Group)) %>% 
  filter(!is.na(Group))

tof <- read_excel(path, sheet = 3, col_types = c("text", "text", "skip", replicate(809, "numeric"))) %>% 
  mutate(Group = as.factor(Group)) %>% 
  filter(!is.na(Group)) %>% 
  na.omit()

common_samples <- biocrates[["Sample Name"]][!grepl(pattern = "^QC ", x = biocrates[["Sample Name"]])]

# Functions

get_group_id <- function(df, n, group_name) {
  if (n == "full") {
    df[df[["Group"]] == group_name, ]
  } else {
    df[sample(which(df[["Group"]] == group_name), n, replace = FALSE), ]
  }
}

perform_smote <- function(data, group_col, original_size) {
  numeric_data <- data %>%
    select(where(is.numeric), !!sym(group_col))
  
  smote_data <- smote(numeric_data, var = group_col, over_ratio = original_size / nrow(data))
  
  smote_data <- as_tibble(smote_data)
  smote_data$`Sample Name` <- paste("Sample", seq(1, nrow(smote_data)))
  smote_data
}

do_ks_test <- function(original_dat, oversampled_dat) {
  sapply(setdiff(colnames(original_dat), c("Group", "Sample Name")), function(ith_metabolite) {
    ks_result <- ks.test(original_dat[[ith_metabolite]], oversampled_dat[[ith_metabolite]])
    c(D = ks_result$statistic, p_value = ks_result$p.value)
  }) %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("metabolite")
}

process_big_data <- function(data) {
  data %>%
    select(metabolite, starts_with("rep")) %>%
    pivot_longer(
      cols = starts_with("rep"),
      names_to = "iteration",
      values_to = "pval",
      names_prefix = "rep"
    )
}


# GC ------------------

part_dat <- gc %>%
  filter(`Sample Name` %in% common_samples) %>%
  filter(!Group %in% c("G2", "G2 (enepdymoma)", "G1", "G3")) %>%
  mutate(Group = if_else(grepl("^G", as.character(Group)), "G4", as.character(Group)))

# Main loop for n_sampled and 100 repetitions
n_sampled_values <- c(10, 20, 30)
n_repeats <- 100
results_list <- vector("list", length(n_sampled_values))

for (n_index in seq_along(n_sampled_values)) {
  n_sampled <- n_sampled_values[n_index]
  rep_results <- vector("list", n_repeats)
  
  for (i in 1:n_repeats) {
    reduced_dat_MT <- get_group_id(part_dat, n_sampled, "MT") %>%
      droplevels() %>% filter(Group == "MT")
    reduced_dat_G4 <- get_group_id(part_dat, n_sampled, "G4") %>%
      droplevels() %>% filter(Group == "G4")
    
    oversampled_data <- list(
      MT = perform_smote(reduced_dat_MT, "Group", nrow(reduced_dat_MT) * (6.8 / (n_sampled / 10))),
      G4 = perform_smote(reduced_dat_G4, "Group", nrow(reduced_dat_G4) * (6.5 / (n_sampled / 10)))
    )
    
    # Original data for KS test
    original_dat_MT <- part_dat %>% filter(Group == "MT")
    original_dat_G4 <- part_dat %>% filter(Group == "G4")
    
    # Perform KS tests
    ks_results_MT <- do_ks_test(original_dat_MT, oversampled_data$MT) %>%
      rename(!!paste0("rep", i, "_D") := D.D, !!paste0("rep", i, "_pval") := p_value)
    ks_results_G4 <- do_ks_test(original_dat_G4, oversampled_data$G4) %>%
      rename(!!paste0("rep", i, "_D") := D.D, !!paste0("rep", i, "_pval") := p_value)
    
    # Store results
    rep_results[[i]] <- list(
      MT = ks_results_MT,
      G4 = ks_results_G4
    )
  }
  
  # Combine all repetitions for MT and G4
  final_results_MT <- reduce(map(rep_results, "MT"), full_join, by = "metabolite")
  final_results_G4 <- reduce(map(rep_results, "G4"), full_join, by = "metabolite")
  
  # Store results for this n_sampled
  results_list[[n_index]] <- list(
    n_sampled = n_sampled,
    MT = final_results_MT,
    G4 = final_results_G4
  )
}

# View results for each n_sampled
names(results_list) <- paste0("n_sampled_", n_sampled_values)

# GC G4 data

# Clean and reshape data
g4_long <- list(
  `n_sampled_10` = process_big_data(results_list[["n_sampled_10"]]$G4) %>% mutate(n_sampled = "10"),
  `n_sampled_20` = process_big_data(results_list[["n_sampled_20"]]$G4) %>% mutate(n_sampled = "20"),
  `n_sampled_30` = process_big_data(results_list[["n_sampled_30"]]$G4) %>% mutate(n_sampled = "30")
) %>%
  bind_rows()

# Summarize data for CI
g4_summary <- g4_long %>%
  group_by(metabolite, n_sampled) %>%
  summarise(
    mean_pval = mean(pval, na.rm = TRUE),
    ci_lower = mean_pval - 1.96 * sd(pval, na.rm = TRUE) / sqrt(n()),
    ci_upper = mean_pval + 1.96 * sd(pval, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# GC MT data

results_list[["n_sampled_10"]]$G4_clean <- results_list[["n_sampled_10"]]$G4 %>%
  select(-matches("_D$"))
results_list[["n_sampled_20"]]$G4_clean <- results_list[["n_sampled_20"]]$G4 %>%
  select(-matches("_D$"))
results_list[["n_sampled_30"]]$G4_clean <- results_list[["n_sampled_30"]]$G4 %>%
  select(-matches("_D$"))


g4_long <- list(
  `n_sampled_10` = process_big_data(results_list[["n_sampled_10"]]$G4_clean) %>% mutate(n_sampled = "10", group = "G4"),
  `n_sampled_20` = process_big_data(results_list[["n_sampled_20"]]$G4_clean) %>% mutate(n_sampled = "20", group = "G4"),
  `n_sampled_30` = process_big_data(results_list[["n_sampled_30"]]$G4_clean) %>% mutate(n_sampled = "30", group = "G4")
) %>%
  bind_rows()

results_list[["n_sampled_10"]]$MT_clean <- results_list[["n_sampled_10"]]$MT %>%
  select(-matches("_D$"))
results_list[["n_sampled_20"]]$MT_clean <- results_list[["n_sampled_20"]]$MT %>%
  select(-matches("_D$"))
results_list[["n_sampled_30"]]$MT_clean <- results_list[["n_sampled_30"]]$MT %>%
  select(-matches("_D$"))

mt_long <- list(
  `n_sampled_10` = process_big_data(results_list[["n_sampled_10"]]$MT_clean) %>%
    select(-matches("_D$")) %>% mutate(n_sampled = "10", group = "MT"),
  `n_sampled_20` = process_big_data(results_list[["n_sampled_20"]]$MT_clean) %>%
    select(-matches("_D$")) %>% mutate(n_sampled = "20", group = "MT"),
  `n_sampled_30` = process_big_data(results_list[["n_sampled_30"]]$MT_clean) %>%
    select(-matches("_D$")) %>% mutate(n_sampled = "30", group = "MT")
) %>%
  bind_rows()

data_summary <- bind_rows(
  g4_long,
  mt_long
) %>%
  group_by(metabolite, n_sampled, group) %>%
  summarise(
    mean_pval = mean(pval, na.rm = TRUE),
    median_pval = median(pval, na.rm = TRUE),
    sd_pval = sd(pval, na.rm = TRUE),
    iqr_pval = IQR(pval, na.rm = TRUE),
    range_pval = paste(range(pval, na.rm = TRUE), collapse = " - "),
    prop_significant = mean(pval < 0.05, na.rm = TRUE),
    mad_pval = mad(pval, na.rm = TRUE),
    .groups = "drop"
  )

# Visualize
ggplot(data_summary, aes(x = n_sampled, y = mean_pval, color = group)) +
  geom_point(alpha = 0.6, position = position_jitter(width = 0.1)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(
    title = "Distribution of Mean P-values of KS test Across n_sampled Groups for G4 and MT data",
    x = "n sampled",
    y = "Mean p-value calculated from 100 results",
    color = "Group"
  ) +
  theme_minimal() +
  scale_y_continuous() +
  expand_limits(y = c(0, 0.8)) +
  scale_color_manual(values = c("G4" = "#ee55f8", "MT" = "#42ff59"))


# Biocrates ------------------

part_dat <- biocrates %>%
  filter(`Sample Name` %in% common_samples) %>%
  filter(!Group %in% c("G2", "G2 (enepdymoma)", "G1", "G3")) %>%
  mutate(Group = if_else(grepl("^G", as.character(Group)), "G4", as.character(Group)))

# Main loop for n_sampled and 100 repetitions
n_sampled_values <- c(10, 20, 30)
n_repeats <- 100
results_list <- vector("list", length(n_sampled_values))

for (n_index in seq_along(n_sampled_values)) {
  n_sampled <- n_sampled_values[n_index]
  rep_results <- vector("list", n_repeats)
  
  for (i in 1:n_repeats) {
    reduced_dat_MT <- get_group_id(part_dat, n_sampled, "MT") %>%
      droplevels() %>% filter(Group == "MT")
    reduced_dat_G4 <- get_group_id(part_dat, n_sampled, "G4") %>%
      droplevels() %>% filter(Group == "G4")
    
    oversampled_data <- list(
      MT = perform_smote(reduced_dat_MT, "Group", nrow(reduced_dat_MT) * (6.8 / (n_sampled / 10))),
      G4 = perform_smote(reduced_dat_G4, "Group", nrow(reduced_dat_G4) * (6.5 / (n_sampled / 10)))
    )
    
    original_dat_MT <- part_dat %>% filter(Group == "MT")
    original_dat_G4 <- part_dat %>% filter(Group == "G4")
    
    ks_results_MT <- do_ks_test(original_dat_MT, oversampled_data$MT) %>%
      rename(!!paste0("rep", i, "_D") := D.D, !!paste0("rep", i, "_pval") := p_value)
    ks_results_G4 <- do_ks_test(original_dat_G4, oversampled_data$G4) %>%
      rename(!!paste0("rep", i, "_D") := D.D, !!paste0("rep", i, "_pval") := p_value)
    
    rep_results[[i]] <- list(
      MT = ks_results_MT,
      G4 = ks_results_G4
    )
  }
  
  final_results_MT <- reduce(map(rep_results, "MT"), full_join, by = "metabolite")
  final_results_G4 <- reduce(map(rep_results, "G4"), full_join, by = "metabolite")
  
  results_list[[n_index]] <- list(
    n_sampled = n_sampled,
    MT = final_results_MT,
    G4 = final_results_G4
  )
}

names(results_list) <- paste0("n_sampled_", n_sampled_values)
bio_results_list <- results_list

bio_results_list[["n_sampled_10"]]$G4_clean <- bio_results_list[["n_sampled_10"]]$G4 %>%
  select(-matches("_D$"))
bio_results_list[["n_sampled_20"]]$G4_clean <- bio_results_list[["n_sampled_20"]]$G4 %>%
  select(-matches("_D$"))
bio_results_list[["n_sampled_30"]]$G4_clean <- bio_results_list[["n_sampled_30"]]$G4 %>%
  select(-matches("_D$"))

g4_long <- list(
  `n_sampled_10` = process_big_data(bio_results_list[["n_sampled_10"]]$G4_clean) %>% mutate(n_sampled = "10"),
  `n_sampled_20` = process_big_data(bio_results_list[["n_sampled_20"]]$G4_clean) %>% mutate(n_sampled = "20"),
  `n_sampled_30` = process_big_data(bio_results_list[["n_sampled_30"]]$G4_clean) %>% mutate(n_sampled = "30")
) %>%
  bind_rows()

g4_summary <- g4_long %>%
  group_by(metabolite, n_sampled) %>%
  summarise(
    mean_pval = mean(pval, na.rm = TRUE),
    ci_lower = mean_pval - 1.96 * sd(pval, na.rm = TRUE) / sqrt(n()),
    ci_upper = mean_pval + 1.96 * sd(pval, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

bio_results_list[["n_sampled_10"]]$MT_clean <- bio_results_list[["n_sampled_10"]]$MT %>%
  select(-matches("_D$"))
bio_results_list[["n_sampled_20"]]$MT_clean <- bio_results_list[["n_sampled_20"]]$MT %>%
  select(-matches("_D$"))
bio_results_list[["n_sampled_30"]]$MT_clean <- bio_results_list[["n_sampled_30"]]$MT %>%
  select(-matches("_D$"))

mt_long <- list(
  `n_sampled_10` = process_big_data(bio_results_list[["n_sampled_10"]]$MT_clean) %>% mutate(n_sampled = "10"),
  `n_sampled_20` = process_big_data(bio_results_list[["n_sampled_20"]]$MT_clean) %>% mutate(n_sampled = "20"),
  `n_sampled_30` = process_big_data(bio_results_list[["n_sampled_30"]]$MT_clean) %>% mutate(n_sampled = "30")
) %>%
  bind_rows()

mt_summary <- mt_long %>%
  group_by(metabolite, n_sampled) %>%
  summarise(
    mean_pval = mean(pval, na.rm = TRUE),
    ci_lower = mean_pval - 1.96 * sd(pval, na.rm = TRUE) / sqrt(n()),
    ci_upper = mean_pval + 1.96 * sd(pval, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

g4_long <- list(
  `n_sampled_10` = process_big_data(bio_results_list[["n_sampled_10"]]$G4_clean) %>% mutate(n_sampled = "10", group = "G4"),
  `n_sampled_20` = process_big_data(bio_results_list[["n_sampled_20"]]$G4_clean) %>% mutate(n_sampled = "20", group = "G4"),
  `n_sampled_30` = process_big_data(bio_results_list[["n_sampled_30"]]$G4_clean) %>% mutate(n_sampled = "30", group = "G4")
) %>%
  bind_rows()

mt_long <- list(
  `n_sampled_10` = process_big_data(bio_results_list[["n_sampled_10"]]$MT_clean) %>% mutate(n_sampled = "10", group = "MT"),
  `n_sampled_20` = process_big_data(bio_results_list[["n_sampled_20"]]$MT_clean) %>% mutate(n_sampled = "20", group = "MT"),
  `n_sampled_30` = process_big_data(bio_results_list[["n_sampled_30"]]$MT_clean) %>% mutate(n_sampled = "30", group = "MT")
) %>%
  bind_rows()

bio_ks_data_summary <- bind_rows(
  g4_long,
  mt_long
) %>%
  group_by(metabolite, n_sampled, group) %>%
  summarise(
    mean_pval = mean(pval, na.rm = TRUE),
    ci_lower = mean_pval - 1.96 * sd(pval, na.rm = TRUE) / sqrt(n()),
    ci_upper = mean_pval + 1.96 * sd(pval, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Visualize
ggplot(bio_ks_data_summary, aes(x = n_sampled, y = mean_pval, color = group)) +
  geom_point(alpha = 0.6, position = position_jitter(width = 0.1)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(
    title = "Distribution of Mean P-values of KS test Across n_sampled Groups for Biocrates data",
    x = "n sampled",
    y = "Mean p-value calculated from 100 results",
    color = "Group"
  ) +
  theme_minimal() +
  scale_y_continuous() +
  expand_limits(y = c(0, 0.8)) +
  scale_color_manual(values = c("G4" = "#ee55f8", "MT" = "#42ff59"))

# QTOF ----------------------

part_dat <- tof %>%
  filter(`Sample Name` %in% common_samples) %>%
  filter(!Group %in% c("G2", "G2 (enepdymoma)", "G1", "G3")) %>%
  mutate(Group = if_else(grepl("^G", as.character(Group)), "G4", as.character(Group)))

# Main loop for n_sampled and 100 repetitions
n_sampled_values <- c(10, 20, 30)
n_repeats <- 100
results_list <- vector("list", length(n_sampled_values))

for (n_index in seq_along(n_sampled_values)) {
  n_sampled <- n_sampled_values[n_index]
  rep_results <- vector("list", n_repeats)
  
  for (i in 1:n_repeats) {
    reduced_dat_MT <- get_group_id(part_dat, n_sampled, "MT") %>%
      droplevels() %>% filter(Group == "MT")
    reduced_dat_G4 <- get_group_id(part_dat, n_sampled, "G4") %>%
      droplevels() %>% filter(Group == "G4")
    
    oversampled_data <- list(
      MT = perform_smote(reduced_dat_MT, "Group", nrow(reduced_dat_MT) * (6.8 / (n_sampled / 10))),
      G4 = perform_smote(reduced_dat_G4, "Group", nrow(reduced_dat_G4) * (6.5 / (n_sampled / 10)))
    )
    
    original_dat_MT <- part_dat %>% filter(Group == "MT")
    original_dat_G4 <- part_dat %>% filter(Group == "G4")
    
    ks_results_MT <- do_ks_test(original_dat_MT, oversampled_data$MT) %>%
      rename(!!paste0("rep", i, "_D") := D.D, !!paste0("rep", i, "_pval") := p_value)
    ks_results_G4 <- do_ks_test(original_dat_G4, oversampled_data$G4) %>%
      rename(!!paste0("rep", i, "_D") := D.D, !!paste0("rep", i, "_pval") := p_value)
    
    rep_results[[i]] <- list(
      MT = ks_results_MT,
      G4 = ks_results_G4
    )
  }
  
  final_results_MT <- reduce(map(rep_results, "MT"), full_join, by = "metabolite")
  final_results_G4 <- reduce(map(rep_results, "G4"), full_join, by = "metabolite")
  
  results_list[[n_index]] <- list(
    n_sampled = n_sampled,
    MT = final_results_MT,
    G4 = final_results_G4
  )
}

# View results for each n_sampled
names(results_list) <- paste0("n_sampled_", n_sampled_values)

tof_results_list <- results_list

tof_results_list[["n_sampled_10"]]$G4_clean <- tof_results_list[["n_sampled_10"]]$G4 %>%
  select(-matches("_D$"))
tof_results_list[["n_sampled_20"]]$G4_clean <- tof_results_list[["n_sampled_20"]]$G4 %>%
  select(-matches("_D$"))
tof_results_list[["n_sampled_30"]]$G4_clean <- tof_results_list[["n_sampled_30"]]$G4 %>%
  select(-matches("_D$"))

# Clean and reshape data
g4_long <- list(
  `n_sampled_10` = process_big_data(tof_results_list[["n_sampled_10"]]$G4_clean) %>% mutate(n_sampled = "10"),
  `n_sampled_20` = process_big_data(tof_results_list[["n_sampled_20"]]$G4_clean) %>% mutate(n_sampled = "20"),
  `n_sampled_30` = process_big_data(tof_results_list[["n_sampled_30"]]$G4_clean) %>% mutate(n_sampled = "30")
) %>%
  bind_rows()

tof_results_list[["n_sampled_10"]]$MT_clean <- tof_results_list[["n_sampled_10"]]$MT %>%
  select(-matches("_D$"))
tof_results_list[["n_sampled_20"]]$MT_clean <- tof_results_list[["n_sampled_20"]]$MT %>%
  select(-matches("_D$"))
tof_results_list[["n_sampled_30"]]$MT_clean <- tof_results_list[["n_sampled_30"]]$MT %>%
  select(-matches("_D$"))

# Clean and reshape data
mt_long <- list(
  `n_sampled_10` = process_big_data(tof_results_list[["n_sampled_10"]]$MT_clean) %>% mutate(n_sampled = "10"),
  `n_sampled_20` = process_big_data(tof_results_list[["n_sampled_20"]]$MT_clean) %>% mutate(n_sampled = "20"),
  `n_sampled_30` = process_big_data(tof_results_list[["n_sampled_30"]]$MT_clean) %>% mutate(n_sampled = "30")
) %>%
  bind_rows()


g4_long <- list(
  `n_sampled_10` = process_big_data(tof_results_list[["n_sampled_10"]]$G4_clean) %>% mutate(n_sampled = "10", group = "G4"),
  `n_sampled_20` = process_big_data(tof_results_list[["n_sampled_20"]]$G4_clean) %>% mutate(n_sampled = "20", group = "G4"),
  `n_sampled_30` = process_big_data(tof_results_list[["n_sampled_30"]]$G4_clean) %>% mutate(n_sampled = "30", group = "G4")
) %>%
  bind_rows()

mt_long <- list(
  `n_sampled_10` = process_big_data(tof_results_list[["n_sampled_10"]]$MT_clean) %>% mutate(n_sampled = "10", group = "MT"),
  `n_sampled_20` = process_big_data(tof_results_list[["n_sampled_20"]]$MT_clean) %>% mutate(n_sampled = "20", group = "MT"),
  `n_sampled_30` = process_big_data(tof_results_list[["n_sampled_30"]]$MT_clean) %>% mutate(n_sampled = "30", group = "MT")
) %>%
  bind_rows()

# Combine the data and summarize for CI
tof_ks_data_summary <- bind_rows(
  g4_long,
  mt_long
) %>%
  group_by(metabolite, n_sampled, group) %>%
  summarise(
    mean_pval = mean(pval, na.rm = TRUE),
    ci_lower = mean_pval - 1.96 * sd(pval, na.rm = TRUE) / sqrt(n()),
    ci_upper = mean_pval + 1.96 * sd(pval, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


# Visualize the data
ggplot(tof_ks_data_summary, aes(x = n_sampled, y = mean_pval, color = group)) +
  geom_point(alpha = 0.6, position = position_jitter(width = 0.1)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(
    title = "Distribution of Mean P-values of KS test Across n_sampled Groups for QTOF data",
    x = "n sampled",
    y = "Mean p-value calculated from 100 results",
    color = "Group"
  ) +
  theme_minimal() +
  scale_y_continuous() +
  expand_limits(y = c(0, 0.8)) +
  scale_color_manual(values = c("G4" = "#ee55f8", "MT" = "#42ff59"))