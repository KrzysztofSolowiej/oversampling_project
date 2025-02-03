library(tidyverse)
library(readxl)
library(tidymodels)
library(ranger)
library(missMDA)
library(pbapply)

n_rep <- 500

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

do_testing <- function(oversampled_data, cnt, n_sampled, rep_id, ith_dataset_name, ith_oa_name) {
  if(inherits(oversampled_data, "try-error")) {
    list(wilcoxon = data.frame(metabolite = NA,
                               apval = NA,
                               group_name = c("MT", "G4"),
                               n_sampled = as.character(n_sampled),
                               rep_id = rep_id, 
                               dataset = ith_dataset_name,
                               oa = ith_oa_name), 
         rf = data.frame(oob = NA, 
                         n_sampled = as.character(n_sampled),
                         rep_id = rep_id, 
                         dataset = ith_dataset_name,
                         group_name = c("MT", "G4"),
                         oa = ith_oa_name))
  } else {
    wilcoxon_res <- lapply(names(oversampled_data), function(ith_dat_name) {
      
      pvals <- do_wilcoxon(oversampled_data[[ith_dat_name]], cnt)
      
      data.frame(name = ith_dat_name, 
                 metabolite = names(pvals),
                 pval = unname(pvals),
                 group_name = ith_dat_name)
    }) %>% 
      bind_rows() %>% 
      mutate(apval = p.adjust(pval, method = "BH"),
             n_sampled = as.character(n_sampled),
             rep_id = rep_id, 
             dataset = ith_dataset_name,
             oa = ith_oa_name) %>% 
      select(-pval)
    
    rf_res <- lapply(names(oversampled_data), function(ith_dat_name) {
      data.frame(oob = do_forest(oversampled_data[[ith_dat_name]], cnt), 
                 n_sampled = as.character(n_sampled),
                 rep_id = rep_id, 
                 dataset = ith_dataset_name,
                 group_name = ith_dat_name,
                 oa = ith_oa_name)
    }) %>% 
      bind_rows()
    
    
    list(wilcoxon = wilcoxon_res, rf = rf_res)
  }
}


dataset_list <- list(gc = gc, tof = tof, biocrates = biocrates)
alg_list <- list(smote = themis::smote, adasyn = themis::adasyn, 
                 bsmote = themis::bsmote, smotenc = themis::smotenc)

os_df <- pblapply(paste0("rep", 1L:n_rep), function(rep_id) {
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

all_rf <- bind_rows(unlisted_data[grepl(pattern = "rf", names(unlisted_data))]) %>% 
  select(group_name, n_sampled, dataset, rep_id, oa, oob)

all_wilcoxon <- bind_rows(unlisted_data[grepl(pattern = "wilcoxon", names(unlisted_data))]) %>% 
  select(group_name, n_sampled, dataset, rep_id, oa, metabolite, apval)

all_wilcoxon_select <- all_wilcoxon %>%
  select(metabolite, group_name, oa, apval, n_sampled, rep_id, dataset)

all_wilcoxon_select_wide <- all_wilcoxon_select %>%
  pivot_wider(
    names_from = rep_id,
    values_from = apval
  )

alpha <- 0.05
final_rep <- paste("rep", n_rep, sep = "")

wilcoxon_summary_stats_better <- all_wilcoxon_select_wide %>%
  mutate(
    mean_p = rowMeans(select(., rep1:final_rep)), 
    median_p = apply(select(., rep1:final_rep), 1, median),
    prop_significant = rowMeans(select(., rep1:final_rep) < alpha),
    sd = apply(select(., rep1:final_rep), 1, sd)
  )