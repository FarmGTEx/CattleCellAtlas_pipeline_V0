#! /work/home/sdxgroup01/00envs/anaconda3/envs/R4.2.0/bin/Rscript

#load packages
library(dplyr)
library(reshape2)
library(tidyr)

files <- list.files("~/dat/pyscenic/sex_no1day_tf",pattern="\\.csv$", full.names= TRUE)
files <- files[-c(17, 18, 21)]
ars_df <- data.frame()
all_res <- data.frame()
for (file in files) {
  data <- read.csv(file)
  long_data <- pivot_longer(data, cols = -X, names_to = "Gender_Celltype", values_to = "ARS")
  long_data <- separate(long_data, Gender_Celltype, into = c("Gender", "Celltype"), sep = "_")
  long_data <- long_data %>%
    mutate(Gender = factor(Gender, levels = c("Male", "Female")))
  t_df <- data.frame()
  for (cell_type in unique(long_data$Celltype)) {
    data_tmp <- long_data[long_data$Celltype == cell_type, ]
    tissue_tmp <- strsplit(as.character(data_tmp[1, 1]), "_")[[1]][2]
    if (length(unique(data_tmp$Gender)) > 1) {
      nor <- shapiro.test(data_tmp$ARS)
      data_tmp1 <- pivot_wider(data_tmp, names_from = Gender, values_from = ARS)
      # res_tmp <- t.test(ARS ~ Gender, data = data_tmp)
      res_tmp <- wilcox.test(data_tmp1$Female, data_tmp1$Male, paired = TRUE)
      t_df <- rbind(t_df, data.frame(Tissue = tissue_tmp, Celltype = cell_type, pvalue = res_tmp$p.value))
    } else {
      t_df <- rbind(t_df, data.frame(Tissue = tissue_tmp, Celltype = cell_type, pvalue = 0.0001))
    }
  }
  # t_df$adjust_pvalue <- p.adjust(t_df$pvalue, method = "bonferroni")
  t_df$adjust_pvalue <- p.adjust(t_df$pvalue, method = "bonferroni")
  # t_df$adjust_pvalue <- p.adjust(t_df$pvalue, method = "BH")
  per <- paste0(sum(t_df$pvalue < 0.05) / length(t_df$pvalue) * 100, "%")
  per1 <- paste0(sum(t_df$adjust_pvalue < 0.05) / length(t_df$pvalue) * 100, "%")
  # per <- paste0(sum(t_df$adjust_pvalue < 0.05) / length(t_df$pvalue) * 100, "%")
  # print(paste0(tissue_tmp, " : ", per, " have sex bias on tf reulgons!"))
  print(paste0(tissue_tmp, " : ", per1, " have sex bias on tf reulgons(adjusted)!"))
  # all_res <- rbind(all_res, t_df)
  # ars_df <- rbind(ars_df, long_data)
}
write.csv(ars_df,"tf_info.csv")


