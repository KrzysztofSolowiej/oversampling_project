if (!requireNamespace("recipes", quietly = TRUE)) {
  install.packages("recipes")
}

library(tidyverse)
library(readxl)
library(dplyr)
library(ggplot2)
library(recipes)
library(themis)
library(Rtsne)

path <- "data/original_matrix.xlsx"

# Read datasets
gc <- read_excel(path, sheet = 1) %>%
  mutate(Group = as.factor(Group)) %>%
  filter(!is.na(Group))

biocrates <- read_excel(path, sheet = 2) %>%
  mutate(Group = as.factor(Group)) %>%
  filter(!is.na(Group))

common_samples <- biocrates[["Sample Name"]][!grepl(pattern = "^QC ", x = biocrates[["Sample Name"]])]

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

# Filter and preprocess data
gc_data <- gc %>%
  filter(`Sample Name` %in% common_samples) %>%
  filter(!Group %in% c("G2", "G2 (enepdymoma)", "G1", "G3")) %>%
  mutate(Group = if_else(grepl("^G", as.character(Group)), "G4", as.character(Group)))

bio_data <- biocrates %>%
  filter(`Sample Name` %in% common_samples) %>%
  filter(!Group %in% c("G2", "G2 (enepdymoma)", "G1", "G3")) %>%
  mutate(Group = if_else(grepl("^G", as.character(Group)), "G4", as.character(Group)))

gc_features <- c("Alanine", "Creatinine", "Glutamic acid", "Glycine", "Isoleucine", "Leucine", "Methionine", "Phenylalanine", "Proline", "Threonine", "trans-4-hydroxy-L-proline", "Tyrosine", "Valine")
bio_features <- c("Ala", "Creatinine", "Glu", "Gly", "Ile", "Leu", "Met", "Phe", "Pro", "Thr", "t4-OH-Pro", "Tyr", "Val")

gc_data_filt <- gc_data %>%
  select(`Sample Name`, Group, all_of(gc_features)) %>%
  mutate(Group = factor(Group, levels = c("G4", "MT", "Con")))

bio_data_filt <- bio_data %>%
  select(`Sample Name`, Group, all_of(bio_features)) %>%
  mutate(Group = factor(Group, levels = c("G4", "MT", "Con")))

feature_mapping <- setNames(bio_features, gc_features)

gc_dat <- gc_data_filt %>%
  rename_with(~ feature_mapping[.x], .cols = all_of(gc_features))
bio_dat <- bio_data_filt

# Correlations ------------------
calculate_correlations <- function(group_name, gc_data, bio_data) {
  gc_group <- gc_data %>% filter(Group == group_name)
  bio_group <- bio_data %>% filter(Group == group_name)
  
  gc_group <- gc_group %>% arrange(`Sample Name`)
  bio_group <- bio_group %>% arrange(`Sample Name`)
  
  metabolites <- setdiff(colnames(gc_group), c("Sample Name", "Group"))
  
  results <- lapply(metabolites, function(metabolite) {
    x <- gc_group[[metabolite]]
    y <- bio_group[[metabolite]]
    if (length(x) != length(y)) {
      stop("Mismatch in lengths of data for metabolite: ", metabolite, " in group: ", group_name)
    }
    test_result <- cor.test(x, y, use = "complete.obs", method = "pearson")
    list(
      Correlation = test_result$estimate,
      P_Value = test_result$p.value
    )
  })
  
  tibble(
    Metabolite = metabolites,
    Correlation = sapply(results, function(res) res$Correlation),
    P_Value = sapply(results, function(res) res$P_Value),
    Group = group_name
  )
}


# Groups to analyze
groups <- unique(gc_data$Group)

all_correlations <- lapply(groups, calculate_correlations, gc_data = gc_dat, bio_data = bio_dat) %>%
  bind_rows()

print(all_correlations, n = 40)


# Visualize ----------------

fixed_limits <- c(-1, 1)

heatmap_data <- all_correlations %>%
  pivot_wider(names_from = Group, values_from = c(Correlation, P_Value)) 

long_heatmap_data <- heatmap_data %>%
  pivot_longer(
    cols = starts_with("Correlation_"), 
    names_to = "Group", 
    names_prefix = "Correlation_", 
    values_to = "Correlation"
  ) %>%
  left_join(
    heatmap_data %>%
      pivot_longer(
        cols = starts_with("P_Value_"), 
        names_to = "Group", 
        names_prefix = "P_Value_", 
        values_to = "P_Value"
      ),
    by = c("Metabolite", "Group")
  )

long_heatmap_data$Group <- factor(long_heatmap_data$Group, levels = c("G4", "MT", "Con"))
long_heatmap_data$Metabolite <- factor(long_heatmap_data$Metabolite, levels = rev(unique(long_heatmap_data$Metabolite)))

ggplot(long_heatmap_data, aes(x = Group, y = Metabolite, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = sprintf("corr=%.2f", Correlation)), color = "black", size = 3, vjust = -0.5) +
  geom_text(aes(label = ifelse(P_Value < 0.05, sprintf("p-val=%.2e", P_Value), "")), 
            color = "black", size = 3, vjust = 1.5) +
  
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limits = fixed_limits
  ) +
  
  theme_minimal() +
  labs(
    title = "Original Metabolite Correlations Between GC and Biocrates",
    x = "Group",
    y = "Metabolite",
    fill = "Correlation"
  )

# Oversampled data - 100 iteration --------------

n_sampled_values <- c(10, 20, 30)
n_repeats <- 100
results_list <- vector("list", length(n_sampled_values))

for (n_index in seq_along(n_sampled_values)) {
  n_sampled <- n_sampled_values[n_index]
  rep_results <- vector("list", n_repeats)
  
  for (i in 1:n_repeats) {
    reduced_gc_g4_dat <- get_group_id(gc_dat, n_sampled, "G4") %>%
      droplevels() %>% filter(Group == "G4")
    reduced_gc_mt_dat <- get_group_id(gc_dat, n_sampled, "MT") %>%
      droplevels() %>% filter(Group == "MT")
    reduced_bio_g4_dat <- bio_dat %>%
      filter(`Sample Name` %in% reduced_gc_g4_dat$'Sample Name') %>% droplevels()
    reduced_bio_mt_dat <- bio_dat %>%
      filter(`Sample Name` %in% reduced_gc_mt_dat$'Sample Name') %>% droplevels()
    
    reduced_gc_g4_dat <- reduced_gc_g4_dat %>%
      arrange(`Sample Name`) %>%
      mutate(Sample = row_number())
    
    reduced_bio_g4_dat <- reduced_bio_g4_dat %>%
      arrange(`Sample Name`) %>%
      mutate(Sample = row_number())
    
    reduced_gc_mt_dat <- reduced_gc_mt_dat %>%
      arrange(`Sample Name`) %>%
      mutate(Sample = row_number())
    
    reduced_bio_mt_dat <- reduced_bio_mt_dat %>%
      arrange(`Sample Name`) %>%
      mutate(Sample = row_number())
    
    gc_g4_ov <- perform_smote(reduced_gc_g4_dat, "Group", nrow(reduced_gc_g4_dat) * (6.5 / (n_sampled / 10)))
    gc_mt_ov <- perform_smote(reduced_gc_mt_dat, "Group", nrow(reduced_gc_mt_dat) * (6.8 / (n_sampled / 10)))
    
    ov_gc_dat <- bind_rows(gc_g4_ov, gc_mt_ov)
    
    bio_g4_ov <- perform_smote(reduced_bio_g4_dat, "Group", nrow(reduced_bio_g4_dat) * (6.5 / (n_sampled / 10)))
    bio_mt_ov <- perform_smote(reduced_bio_mt_dat, "Group", nrow(reduced_bio_mt_dat) * (6.8 / (n_sampled / 10)))
    
    ov_bio_dat <- bind_rows(bio_g4_ov, bio_mt_ov)
    
    groups <- unique(ov_gc_dat$Group)
    
    ov_gc_dat <- ov_gc_dat %>%
      select(!Sample)
    
    ov_bio_dat <- ov_bio_dat %>%
      select(!Sample)
    
    ov_correlations <- lapply(groups, calculate_correlations, gc_data = ov_gc_dat, bio_data = ov_bio_dat) %>%
      bind_rows()
    
    rep_results[[i]] <- list(
      corr_results = ov_correlations
    )
  }
  
  results_list[[n_index]] <- list(
    n_sampled = n_sampled,
    my_results = rep_results
  )
}

process_results <- function(results_element) {
  correlation_data <- lapply(1:length(results_element$my_results), function(iter) {
    data <- results_element$my_results[[iter]]$corr_results
    
    colnames(data)[2] <- paste0("rep", iter, "_Correlation")
    colnames(data)[3] <- paste0("rep", iter, "_PValue")
    
    data[, c("Metabolite", paste0("rep", iter, "_Correlation"), 
             paste0("rep", iter, "_PValue"), "Group")]
  })
  
  final_df <- Reduce(function(df1, df2) {
    merge(df1, df2, by = c("Metabolite", "Group"))
  }, correlation_data)
  
  return(final_df)
}

# Process results for each n_sampled value
results_df1 <- process_results(results_list[[1]])
results_df2 <- process_results(results_list[[2]])
results_df3 <- process_results(results_list[[3]])

# Add mean and median calculations for both correlation and p-values
results_df1 <- results_df1 %>%
  mutate(
    corr_mean = rowMeans(across(starts_with("rep") & ends_with("_Correlation")), na.rm = TRUE),
    corr_median = apply(across(starts_with("rep") & ends_with("_Correlation")), 1, median, na.rm = TRUE),
    pval_mean = rowMeans(across(starts_with("rep") & ends_with("_PValue")), na.rm = TRUE),
    pval_median = apply(across(starts_with("rep") & ends_with("_PValue")), 1, median, na.rm = TRUE)
  )

results_df2 <- results_df2 %>%
  mutate(
    corr_mean = rowMeans(across(starts_with("rep") & ends_with("_Correlation")), na.rm = TRUE),
    corr_median = apply(across(starts_with("rep") & ends_with("_Correlation")), 1, median, na.rm = TRUE),
    pval_mean = rowMeans(across(starts_with("rep") & ends_with("_PValue")), na.rm = TRUE),
    pval_median = apply(across(starts_with("rep") & ends_with("_PValue")), 1, median, na.rm = TRUE)
  )

results_df3 <- results_df3 %>%
  mutate(
    corr_mean = rowMeans(across(starts_with("rep") & ends_with("_Correlation")), na.rm = TRUE),
    corr_median = apply(across(starts_with("rep") & ends_with("_Correlation")), 1, median, na.rm = TRUE),
    pval_mean = rowMeans(across(starts_with("rep") & ends_with("_PValue")), na.rm = TRUE),
    pval_median = apply(across(starts_with("rep") & ends_with("_PValue")), 1, median, na.rm = TRUE)
  )

fixed_limits <- c(-1, 1)

results_df1_heat_dat <- results_df1 %>%
  select(Metabolite, Group, corr_mean, pval_mean)

results_df1_heat_dat$Group <- factor(results_df1_heat_dat$Group, levels = c("G4", "MT"))
results_df1_heat_dat$Metabolite <- factor(results_df1_heat_dat$Metabolite, levels = rev(unique(results_df1_heat_dat$Metabolite)))

ggplot(results_df1_heat_dat, aes(x = Group, y = Metabolite, fill = corr_mean)) +
  geom_tile() +
  geom_text(aes(label = sprintf("corr=%.2f", corr_mean)), color = "black", size = 3, vjust = -0.5) +
  geom_text(aes(label = ifelse(pval_mean < 0.05, sprintf("p-val=%.2e", pval_mean), "")), 
            color = "black", size = 3, vjust = 1.5) +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limits = fixed_limits
  ) +
  theme_minimal() +
  labs(
    title = "Mean Metabolite Correlations between GC and Biocrates oversampled data\nafter oversampling done 100 times (n = 10)",
    x = "Group",
    y = "Metabolite",
    fill = "Correlation"
  )



results_df2_heat_dat <- results_df2 %>%
  select(Metabolite, Group, corr_mean, pval_mean)

results_df2_heat_dat$Group <- factor(results_df2_heat_dat$Group, levels = c("G4", "MT"))
results_df2_heat_dat$Metabolite <- factor(results_df2_heat_dat$Metabolite, levels = rev(unique(results_df2_heat_dat$Metabolite)))


ggplot(results_df2_heat_dat, aes(x = Group, y = Metabolite, fill = corr_mean)) +
  geom_tile() +
  geom_text(aes(label = sprintf("corr=%.2f", corr_mean)), color = "black", size = 3, vjust = -0.5) +
  geom_text(aes(label = ifelse(pval_mean < 0.05, sprintf("p-val=%.2e", pval_mean), "")), 
            color = "black", size = 3, vjust = 1.5) +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limits = fixed_limits
  ) +
  theme_minimal() +
  labs(
    title = "Mean Metabolite Correlations between GC and Biocrates oversampled data\nafter oversampling done 100 times (n = 20)",
    x = "Group",
    y = "Metabolite",
    fill = "Correlation"
  )

results_df3_heat_dat <- results_df3 %>%
  select(Metabolite, Group, corr_mean, pval_mean)

results_df3_heat_dat$Group <- factor(results_df3_heat_dat$Group, levels = c("G4", "MT"))
results_df3_heat_dat$Metabolite <- factor(results_df3_heat_dat$Metabolite, levels = rev(unique(results_df3_heat_dat$Metabolite)))


ggplot(results_df3_heat_dat, aes(x = Group, y = Metabolite, fill = corr_mean)) +
  geom_tile() +
  geom_text(aes(label = sprintf("corr=%.2f", corr_mean)), color = "black", size = 3, vjust = -0.5) +
  geom_text(aes(label = ifelse(pval_mean < 0.05, sprintf("p-val=%.2e", pval_mean), "")), 
            color = "black", size = 3, vjust = 1.5) +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limits = fixed_limits
  ) +
  theme_minimal() +
  labs(
    title = "Mean Metabolite Correlations between GC and Biocrates oversampled data\nafter oversampling done 100 times (n = 30)",
    x = "Group",
    y = "Metabolite",
    fill = "Correlation"
  )
