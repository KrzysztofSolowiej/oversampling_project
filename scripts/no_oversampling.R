library(tidyverse)
library(readxl)
library(tidymodels)
library(ranger)
library(missMDA)
library(pbapply)
library(iml)
library(purrr)

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


bare_sample_workflow <- function(dataset_list, common_samples, n_sampled_list, rep_id) {
  pblapply(names(dataset_list), function(ith_dataset_name) {
    part_dat <- filter(dataset_list[[ith_dataset_name]], 
                       `Sample Name` %in% common_samples) 
    
    cnt <- get_group_id(part_dat, "full", "Con")
    
    lapply(n_sampled_list, function(n_sampled) {
      
      sampled_dat_MT <- rbind(cnt,
                              get_group_id(part_dat, n_sampled, "MT")) %>% 
        droplevels()
      
      sampled_dat_G4 <- rbind(cnt,
                              get_group_id(part_dat, n_sampled, "G4")) %>% 
        droplevels()
      
      bare_samples <- list(MT = select(sampled_dat_MT, -`Sample Name`) %>% filter(Group == "MT"),
                           G4 = select(sampled_dat_G4, -`Sample Name`) %>% filter(Group == "G4"))
      
      do_testing(bare_samples, cnt, n_sampled, rep_id, ith_dataset_name, "no oversampling")
    })
  })
}

dataset_list <- list(gc = gc, tof = tof, biocrates = biocrates)


set.seed(184502)

n_sampled_list <- list(n10 = 10, n20 = 20, n30 = 30, n_full = "full")
reps <- paste0("rep", 1L:n_rep)

bare_sample_results <- pblapply(reps, function(rep_id) {
  bare_sample_workflow(dataset_list, common_samples, n_sampled_list, rep_id)
})

bare_unlisted <- unlist(bare_sample_results, recursive = FALSE) %>% 
  unlist(recursive = FALSE)

names(bare_unlisted) <- paste0(rep(names(n_sampled_list), times = length(reps)), "_", rep(reps, each = length(n_sampled_list)))

all_bare_rf <- bare_unlisted %>%
  map("rf") %>%
  bind_rows() %>%
  mutate(oa = "no oversampling")

all_bare_wilcoxon <- bare_unlisted %>%
  map("wilcoxon") %>%
  bind_rows() %>%
  mutate(oa = "no oversampling")

noova_raw <- all_bare_wilcoxon %>%
  select(metabolite, group_name, apval, n_sampled, rep_id, dataset)

noova_wide <- noova_raw %>%
  pivot_wider(
    names_from = rep_id,
    values_from = apval
  )

alpha <- 0.05
final_rep <- paste("rep", n_rep, sep = "")

summary_stats <- noova_wide %>%
  mutate(
    mean_p = rowMeans(select(., rep1:final_rep)), 
    median_p = apply(select(., rep1:final_rep), 1, median),
    prop_significant = rowMeans(select(., rep1:final_rep) < alpha),
    sd = apply(select(., rep1:final_rep), 1, sd)
  )