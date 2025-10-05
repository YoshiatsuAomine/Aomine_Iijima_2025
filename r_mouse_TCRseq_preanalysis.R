# 必要なパッケージをインストール
install.packages("tidyverse")
install.packages("readr")


# 必要なパッケージを読み込み
library(tidyverse)
library(readr)

setwd("/Volumes/AOMINE-E1TB/基盤＿飯島先生/20250520_TCRab-seq_All_Aomine/03mixcr_SMARTer")

#　ここを編集する
CD4orCD8 <- "CD4"
AorB <- "B"

# `CD4orCD8`と同じ方のリストを使用する
namelist <- c("CD4_Br_S6", "CD4_SK_S8", "CD4_SP_S1", "CD4_Stem_S7", "CD4_TG_S10", "Na_SP_CD4_S32", "Na_TG_CD4_S30","CD4_Br_TG_S6_L001")
#namelist <- c("CD8_Br_S1", "CD8_SK_S4", "CD8_SP_S5", "CD8_Stem_S3", "CD8_TG_S2", "Na_SP_CD8_S33", "Na_TG_CD8_S31","Tet_TG+Br_CD8_S1_L001")


# データを格納するリスト
data_list <- list()

# ファイルを読み込むループ
for (name in namelist) {
  file_path <- paste0(name, "_output.clones_TR", AorB, ".tsv")  # ファイル名を作成
  if (file.exists(file_path)) {  # ファイルが存在するか確認
    temp_data <- read_tsv(file_path, col_types = cols(.default = "c")) %>%
    mutate(Sample = name)  # サンプル名を追加
    data_list[[name]] <- temp_data
  } else {
    warning(paste("ファイルが見つかりません:", file_path))
  }
}

# data_list の各データフレームに対して処理を適用
data_list <- map(data_list, function(df) {
  df %>%
    select(allVHitsWithScore, allJHitsWithScore, aaSeqCDR3, nSeqCDR3,readCount) %>%
    mutate(
      allVHitsWithScore = str_remove(allVHitsWithScore, "\\*.*$"),  # *以降を削除
      allJHitsWithScore = str_remove(allJHitsWithScore, "\\*.*$"),  # *以降を削除
      V_J_CDR3aa = paste(allVHitsWithScore, allJHitsWithScore, aaSeqCDR3, sep = "_"),
      readCount = as.numeric(readCount)
    )%>%
    select(V_J_CDR3aa, allVHitsWithScore, allJHitsWithScore, aaSeqCDR3, nSeqCDR3, readCount)  # V_J_CDR3aa を最前列に移動
})

# すべてのデータを結合
merged_data <- bind_rows(data_list)

# V, J, CDR3の情報を抽出してキーを作成
merged_data <- merged_data %>%
  mutate(
    allVHitsWithScore = str_remove(allVHitsWithScore, "\\*.*$"),  # *以降を削除
    allJHitsWithScore = str_remove(allJHitsWithScore, "\\*.*$"),  # *以降を削除
    V_J_CDR3aa = paste(allVHitsWithScore, allJHitsWithScore, aaSeqCDR3, sep = "_")
  )%>%
　select(V_J_CDR3aa, allVHitsWithScore, allJHitsWithScore, aaSeqCDR3, nSeqCDR3)

# 重複を除去（ユニークなV_J_CDR3リスト）
unique_vj_cdr3 <- distinct(merged_data, V_J_CDR3aa, .keep_all = TRUE)

#重複を排除し、readCount を足し算した値を残す
for (sample in names(data_list)) {
  temp_df <- data_list[[sample]] %>%
    group_by(V_J_CDR3aa) %>%  #`V_J_CDR3aa` ごとにまとめる
    summarise(readCount = sum(readCount, na.rm = TRUE), .groups = "drop") #`readCount` を合計

  # 確認用の出力（必要ならコメントアウト）
  print(head(temp_df))  # 最初の数行を表示
  
  data_list[[sample]] <- temp_df  #集約後のデータを `data_list` に上書き保存
}

# データの確認
head(merged_data)

# 各サンプルの readCount を `unique_vj_cdr3` に追加
for (sample in names(data_list)) {
  temp_df <- data_list[[sample]] %>%
  select(V_J_CDR3aa, readCount) %>%
  rename(!!sample := readCount)  #`readCount` をサンプル名に変更
  
  unique_vj_cdr3 <- unique_vj_cdr3 %>%
  left_join(temp_df, by = "V_J_CDR3aa")
}

# 結果を確認
head(unique_vj_cdr3)

write_csv(unique_vj_cdr3, paste0("output_",CD4orCD8, "_TCR-", AorB, ".csv"))

