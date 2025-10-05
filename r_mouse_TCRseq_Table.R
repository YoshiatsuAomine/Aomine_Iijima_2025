install.packages("ggplot2")  
install.packages("dplyr")    
install.packages("crayon")
install.packages("gridExtra")
install.packages("scales")



library(ggplot2)
library(dplyr)
library(crayon)
library(gridExtra)
library(scales)

base_dir <- "/Volumes/AOMINE-E1TB//TCRab-seq_All_Aomine/04R"

setwd(base_dir)

tcr_options <- c("A", "B")
namelists <- list(
  CD4 = c("CD4_Br_S6", "CD4_SK_S8", "CD4_SP_S1", "CD4_Stem_S7","CD4_TG_S10", "Na_SP_CD4_S32", "Na_TG_CD4_S30", "CD4_Br_TG_S6_L001"),
  CD8 = c("CD8_Br_S1", "CD8_SK_S4", "CD8_SP_S5", "CD8_Stem_S3","CD8_TG_S2", "Na_SP_CD8_S33", "Na_TG_CD8_S31", "Tet_TGandBr_CD8_S1_L001")
)

excluded_patterns <- list(
  CD4 = list(
    c(),                                        
    c("Na_SP_CD4_S32"),                         
    c("Na_TG_CD4_S30"),                         
    c("Na_SP_CD4_S32", "Na_TG_CD4_S30")        
  ),
  CD8 = list(
    c(),                                        
    c("Na_SP_CD8_S33"),                         
    c("Na_TG_CD8_S31"),                         
    c("Na_SP_CD8_S33", "Na_TG_CD8_S31")         
  )
)



for (TCR in tcr_options) {
  for (CD in names(namelists)) {
    namelist <- namelists[[CD]]
    
    for (Excluded_sample in excluded_patterns[[CD]]) {
      
df <- read.csv(paste0("output_",CD,"_TCR-",TCR,".csv"), stringsAsFactors = FALSE)

df_filtered <- df[ rowSums(df[Excluded_sample], na.rm = TRUE) == 0, ]

sample_cols <- setdiff(colnames(df_filtered), c("V_J_CDR3aa", "allVHitsWithScore", "allJHitsWithScore", "aaSeqCDR3", "nSeqCDR3", Excluded_sample))

results <- data.frame(Sample1 = character(),
                      Sample2 = character(),
                      Shared_Count = numeric(),
                      Unique_Count = numeric(),
                      stringsAsFactors = FALSE)


for (i in sample_cols) {
  for (j in sample_cols) {
    if (i != j) {
      shared <- subset(df_filtered,
                       !is.na(df_filtered[[i]]) & df_filtered[[i]] > 0 &
                         !is.na(df_filtered[[j]]) & df_filtered[[j]] > 0)
      
      unique <- subset(df_filtered,
                       !is.na(df_filtered[[i]]) & df_filtered[[i]] > 0 &
                         (is.na(df_filtered[[j]]) | df_filtered[[j]] == 0))
      
      results <- rbind(results, data.frame(Sample1 = i,
                                           Sample2 = j,
                                           Sum_of_Sample1_Counts_Common_to_Sample2 = sum(shared[[i]], na.rm = TRUE),
                                           Sum_of_Sample1_Counts_Unique = sum(unique[[i]], na.rm = TRUE)))
    }
  }
}


results <- results %>%
  mutate(Common_ratio = percent(Sum_of_Sample1_Counts_Common_to_Sample2 / (Sum_of_Sample1_Counts_Common_to_Sample2 + Sum_of_Sample1_Counts_Unique), accuracy = 0.1))
results <- results %>%
  mutate(Unique_ratio = percent(Sum_of_Sample1_Counts_Unique / (Sum_of_Sample1_Counts_Common_to_Sample2 + Sum_of_Sample1_Counts_Unique), accuracy = 0.1))

print(results)

excluded_str <- paste(Excluded_sample, collapse = "_")

if (nzchar(excluded_str)) {
  excluded_str <- paste0("_", excluded_str)
}

output_filename <- paste0("Common_Unique_Table_", CD, "_TCR-", TCR, excluded_str, ".csv")
write.csv(results, file = output_filename, row.names = FALSE)

    }
  }
}
