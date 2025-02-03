if (!requireNamespace("caret", quietly = TRUE)) {
  install.packages("caret")
}

library(tidyverse)
library(readxl)
library(tidymodels)
library(ranger)
library(missMDA)
library(pbapply)
library(iml)
library(caret)
library(doParallel)
library(ggplot2)
library(cowplot)

cores <- parallel::detectCores() - 1
cl <- makeCluster(cores)
registerDoParallel(cl)

path <- "data/original_matrix.xlsx"

gc <- read_excel(path, sheet = 1) %>%
  filter(!Group %in% c("G2", "G2 (enepdymoma)", "G1", "G3")) %>%
  mutate(Group = if_else(grepl("^G", as.character(Group)), "G4", as.character(Group))) %>%
  mutate(Group = as.factor(Group)) %>%
  filter(!is.na(Group))

biocrates <- read_excel(path, sheet = 2) %>%
  filter(!Group %in% c("G2", "G2 (enepdymoma)", "G1", "G3")) %>%
  mutate(Group = if_else(grepl("^G", as.character(Group)), "G4", as.character(Group))) %>%
  mutate(Group = as.factor(Group)) %>%
  filter(!is.na(Group))

tof <- read_excel(path, sheet = 3, col_types = c("text", "text", "skip", replicate(809, "numeric"))) %>%
  filter(!Group %in% c("G2", "G2 (enepdymoma)", "G1", "G3")) %>%
  mutate(Group = if_else(grepl("^G", as.character(Group)), "G4", as.character(Group))) %>%
  mutate(Group = as.factor(Group)) %>%
  filter(!is.na(Group)) %>%
  na.omit()

common_samples <- biocrates[["Sample Name"]][!grepl(pattern = "^QC ", x = biocrates[["Sample Name"]])]

get_group_id <- function(df, n, group_name) {
  if(n == "full") {
    df[df[["Group"]] == group_name, ]
  } else {
    df[sample(which(df[["Group"]] == group_name), n, replace = FALSE), ]
  }
}

set.seed(184502)

do_wilcoxon <- function(dat, cnt) {
  sapply(setdiff(colnames(dat), "Group"), function(ith_metabolite) {
    wilcox.test(dat[[ith_metabolite]], cnt[[ith_metabolite]])[["p.value"]]
  })
}

do_forest <- function(dat, cnt) {
  full_dat <- droplevels(rbind(dat, cnt[setdiff(colnames(cnt), "Sample Name")]))
  ranger(x = full_dat[setdiff(colnames(full_dat), "Group")],
         y = full_dat[["Group"]], oob.error = TRUE)[["prediction.error"]]
}

do_svm <- function(dat, cnt) {
  full_dat <- droplevels(bind_rows(dat, select(cnt, -`Sample Name`)))
  x <- as.data.frame(full_dat[setdiff(colnames(full_dat), "Group")])
  y <- full_dat[["Group"]]
  train_control <- trainControl(method="cv",
                                number=5,
                                allowParallel = TRUE)
  
  svm_model <- caret::train(
    x = x, y = y,
    method = "svmRadial",
    trControl = train_control,
    tuneLength = 5
  )
  
  1 - max(svm_model$results$Accuracy)  # Out-of-sample error
}

do_testing <- function(oversampled_data, cnt, n_sampled, rep_id, ith_dataset_name, ith_oa_name) {
  if (inherits(oversampled_data, "try-error")) {
    list(
      svm = data.frame(
        oos = NA, n_sampled = as.character(n_sampled),
        rep_id = rep_id, dataset = ith_dataset_name,
        group_name = c("MT", "G4"), oa = ith_oa_name
      )
    )
  } else {
    svm_res <- lapply(names(oversampled_data), function(ith_dat_name) {
      data.frame(
        oos = do_svm(oversampled_data[[ith_dat_name]], cnt),
        n_sampled = as.character(n_sampled), rep_id = rep_id,
        dataset = ith_dataset_name, group_name = ith_dat_name,
        oa = ith_oa_name
      )
    }) %>%
      bind_rows()
    
    list(svm = svm_res)
  }
}


dataset_list <- list(gc = gc, tof = tof, biocrates = biocrates)
alg_list <- list(smote = themis::smote, adasyn = themis::adasyn, 
                 bsmote = themis::bsmote, smotenc = themis::smotenc)

os_df <- pblapply(paste0("rep", 1L:500), function(rep_id) {
  lapply(names(dataset_list), function(ith_dataset_name) {
    part_dat <- filter(dataset_list[[ith_dataset_name]], 
                       `Sample Name` %in% common_samples) 
    
    cnt <- get_group_id(part_dat, "full", "Con")
    
    lapply(list(n10 = 10, n20 = 20, n30 = 30, n_full = "full"), function(n_sampled) {
      
      sampled_dat_MT <- rbind(cnt,
                              get_group_id(part_dat, n_sampled, "MT")) %>% 
        droplevels()

      sampled_dat_G4 <- rbind(cnt,
                              get_group_id(part_dat, n_sampled, "G4")) %>% 
        droplevels()

      if(n_sampled == "full") {
        oversampled_dfs <- list(list(MT = select(sampled_dat_MT, -`Sample Name`) %>% 
                                  filter(Group == "MT"),
                                G4 = select(sampled_dat_G4, -`Sample Name`) %>% 
                                  filter(Group == "G4"))) %>% 
          setNames("no oversampling")
      } else {
        oversampled_dfs <- lapply(names(alg_list), function(ith_oa_name) {
          ith_oa <- alg_list[[ith_oa_name]]
          
          try({
            list(MT = ith_oa(select(sampled_dat_MT, -`Sample Name`), var = "Group", over_ratio = 1.09) %>% 
                   filter(Group == "MT"),
                 G4 = ith_oa(select(sampled_dat_G4, -`Sample Name`), var = "Group", over_ratio = 1.14) %>% 
                   filter(Group == "G4")) 
          }, silent = TRUE)
        }) %>% 
          setNames(names(alg_list))
      }
      
      lapply(names(oversampled_dfs), function(ith_oa_name)
        do_testing(oversampled_dfs[[ith_oa_name]], cnt, n_sampled, rep_id, ith_dataset_name, ith_oa_name)
      )
    })
  })
})

unlisted_data <- unlist(os_df, recursive = FALSE) %>% 
  unlist(recursive = FALSE) %>% 
  unlist(recursive = FALSE) %>% 
  unlist(recursive = FALSE)


all_svm <- bind_rows(unlisted_data[grepl(pattern = "svm", names(unlisted_data))]) %>%
  select(group_name, n_sampled, dataset, rep_id, oa, oos)

stopCluster(cl)
registerDoSEQ()

mean_oos_data <- all_svm %>%
  group_by(group_name, n_sampled, dataset, oa) %>%
  summarize(mean_oos = mean(oos, na.rm = TRUE)) %>%
  ungroup()


mean_oos_data <- mean_oos_data %>%
  mutate(dataset = recode(dataset, "tof" = "QTOF")) %>%
  mutate(dataset = recode(dataset, "gc" = "GC")) %>%
  mutate(dataset = recode(dataset, "biocrates" = "Biocrates"))

mean_oos_data <- mean_oos_data %>%
  mutate(oa = recode(oa, "smote" = "SMOTE")) %>%
  mutate(oa = recode(oa, "bsmote" = "BSMOTE")) %>%
  mutate(oa = recode(oa, "smotenc" = "SMOTENC")) %>%
  mutate(oa = recode(oa, "adasyn" = "ADASYN"))

custom_colors <- c("Biocrates" = "#b5fdd5", "GC" = "#ffd893", "QTOF" = "#e798ff")

plot_1_data <- mean_oos_data %>%
  filter(n_sampled != "full")

plot_2_data <- mean_oos_data %>%
  filter(n_sampled == "full")

plot_1_data$mean_oos <- as.numeric(plot_1_data$mean_oos)
plot_1_data$group_name <- factor(plot_1_data$group_name, levels = c("G4", "MT"))
plot_1_data$oa <- factor(plot_1_data$oa, levels = c("SMOTE", "BSMOTE", "SMOTENC", "ADASYN"))
plot_1_data$dataset <- factor(plot_1_data$dataset, levels = c("Biocrates", "GC", "QTOF"))

plot_2_data$mean_oos <- as.numeric(plot_2_data$mean_oos)
plot_2_data$group_name <- factor(plot_2_data$group_name, levels = c("G4", "MT"))
plot_2_data$dataset <- factor(plot_2_data$dataset, levels = c("Biocrates", "GC", "QTOF"))

cow_1 <- ggplot(plot_1_data, aes(x = oa, y = mean_oos, fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = round(mean_oos, 2)),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 2) +
  facet_grid(group_name ~ n_sampled, scales = "fixed") +
  labs(title = "OOS (non-linear kernel) of Algorithms by Tumor Type, Matrix Number and Source \ncv = 5",
       x = "Oversampled sets",
       y = "Mean OOS Error",
       fill = "Source") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 0.25), expand = expansion(mult = c(0, 0.2))) +
  geom_vline(xintercept = c(4.5, 7.5), linetype = "dashed", color = "black", size = 0.5) +
  scale_fill_manual(values = custom_colors)

cow_2 <- ggplot(plot_2_data, aes(x = dataset, y = mean_oos, fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = round(mean_oos, 2)),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 2) +
  facet_wrap(~ group_name, ncol = 1, scales = "fixed") +
  labs(title = " \n ",
       x = "Full set",
       y = NULL,
       fill = "Source") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text(margin = margin(t = 36)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        strip.placement = "outside",
        plot.title = element_text(size = 16)) +
  scale_y_continuous(limits = c(0, 0.25), expand = expansion(mult = c(0, 0.2))) +
  scale_fill_manual(values = custom_colors)

cow <- plot_grid(cow_1, cow_2, ncol = 2, rel_widths = c(3.4, 1), axis = "tb")

print(cow)
