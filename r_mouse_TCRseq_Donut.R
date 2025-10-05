# 必要なパッケージをインストール
install.packages("ggplot2")  # グラフ描画用
install.packages("dplyr")    # データ操作用
install.packages("crayon")
install.packages("gridExtra")

# パッケージの読み込み
library(ggplot2)
library(dplyr)
library(crayon)
library(gridExtra)

# ベースディレクトリ（保存したいフォルダのパスを指定）
base_dir <- "/Volumes/AOMINE-E1TB/基盤＿飯島先生/20250520_TCRab-seq_All_Aomine/04R"

#ディレクトリ移動
setwd(base_dir)

#-----変数決定-----
TCR <- "A"
CD <- "CD4"
namelist <- c("CD4_Br_S6", "CD4_SK_S8", "CD4_SP_S1", "CD4_Stem_S7", "CD4_TG_S10", "Na_SP_CD4_S32", "Na_TG_CD4_S30","CD4_Br_TG_S6_L001")
sample3 <- "Na_SP_CD4_S32"
sample4 <- "Na_TG_CD4_S30"
#CD <- "CD8"
#namelist <- c("CD8_Br_S1", "CD8_SK_S4", "CD8_SP_S5", "CD8_Stem_S3", "CD8_TG_S2", "Na_SP_CD8_S33", "Na_TG_CD8_S31","Tet_TGandBr_CD8_S1_L001") #特殊文字は変更、＋→and
#sample3 <- "Na_SP_CD8_S33"
#sample4 <- "Na_TG_CD8_S31"

# データの読み込み
df <- read.csv(paste0("output_",CD,"_TCR-",TCR,".csv"), stringsAsFactors = FALSE)

#------- conc = sample1, sample2 -------
# 保存用フォルダの指定（新規作成するフォルダ名）
output_folder <- file.path(base_dir, paste0(CD, "_TCR-", TCR))

# フォルダが存在しない場合は作成
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

#データを抽出し保存
for (sample1 in namelist) {
  for (sample2 in namelist) {
    
    # sample1 で発現している かつ sample2 で発現している TCRクローン
    df1 <- df %>%
      filter(!is.na(.data[[sample1]]) & .data[[sample1]] > 0, 
             !is.na(.data[[sample2]]) & .data[[sample2]] > 0) %>%
      select(V_J_CDR3aa,all_of(sample2))
    
    # sample1 で発現していない かつ sample2 で発現している TCRクローン
    df2 <- df %>%
      filter((is.na(.data[[sample1]]) | .data[[sample1]] == 0), 
             !is.na(.data[[sample2]]) & .data[[sample2]] > 0) %>%
      select(V_J_CDR3aa,all_of(sample2))
    
    # 結果をCSVに保存
    write.csv(df1, file.path(output_folder, paste0("df1_", sample1, "_exp_--", sample2, "_exp_", TCR, ".csv")), row.names = FALSE)
    write.csv(df2, file.path(output_folder, paste0("df2_", sample1, "_non_--", sample2, "_exp_", TCR, ".csv")), row.names = FALSE)
  }
}

#------- conc = sample1, sample2 -sample3 -------
# 保存用フォルダの指定（新規作成するフォルダ名）
output_folder <- file.path(base_dir, paste0(CD, "_TCR-", TCR, "_-", sample3))

# フォルダが存在しない場合は作成
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

#データを抽出し保存
for (sample1 in namelist) {
  for (sample2 in namelist) {
    
    # sample1 で発現している かつ sample2 で発現している TCRクローン
    df1 <- df %>%
      filter(!is.na(.data[[sample1]]) & .data[[sample1]] > 0, 
             !is.na(.data[[sample2]]) & .data[[sample2]] > 0,
             (is.na(.data[[sample3]]) | .data[[sample3]] == 0)) %>%
      select(V_J_CDR3aa,all_of(sample2))
    
    # sample1 で発現していない かつ sample2 で発現している TCRクローン
    df2 <- df %>%
      filter((is.na(.data[[sample1]]) | .data[[sample1]] == 0), 
             !is.na(.data[[sample2]]) & .data[[sample2]] > 0,
             (is.na(.data[[sample3]]) | .data[[sample3]] == 0)) %>%
      select(V_J_CDR3aa,all_of(sample2))
    
    # 結果をCSVに保存
    write.csv(df1, file.path(output_folder, paste0("df1_", sample1, "_exp_--", sample2, "_exp_", TCR, ".csv")), row.names = FALSE)
    write.csv(df2, file.path(output_folder, paste0("df2_", sample1, "_non_--", sample2, "_exp_", TCR, ".csv")), row.names = FALSE)
  }
}

#------- conc = sample1, sample2 -sample4 -------
# 保存用フォルダの指定（新規作成するフォルダ名）
output_folder <- file.path(base_dir, paste0(CD, "_TCR-", TCR, "_-", sample4))

# フォルダが存在しない場合は作成
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

#データを抽出し保存
for (sample1 in namelist) {
  for (sample2 in namelist) {
    
    # sample1 で発現している かつ sample2 で発現している TCRクローン
    df1 <- df %>%
      filter(!is.na(.data[[sample1]]) & .data[[sample1]] > 0, 
             !is.na(.data[[sample2]]) & .data[[sample2]] > 0,
             (is.na(.data[[sample4]]) | .data[[sample4]] == 0)) %>%
      select(V_J_CDR3aa,all_of(sample2))
    
    # sample1 で発現していない かつ sample2 で発現している TCRクローン
    df2 <- df %>%
      filter((is.na(.data[[sample1]]) | .data[[sample1]] == 0), 
             !is.na(.data[[sample2]]) & .data[[sample2]] > 0,
             (is.na(.data[[sample4]]) | .data[[sample4]] == 0)) %>%
      select(V_J_CDR3aa,all_of(sample2))
    
    # 結果をCSVに保存
    write.csv(df1, file.path(output_folder, paste0("df1_", sample1, "_exp_--", sample2, "_exp_", TCR, ".csv")), row.names = FALSE)
    write.csv(df2, file.path(output_folder, paste0("df2_", sample1, "_non_--", sample2, "_exp_", TCR, ".csv")), row.names = FALSE)
  }
}

#------- conc = sample1, sample2 -sample3 -sample4 -------
# 保存用フォルダの指定（新規作成するフォルダ名）
output_folder <- file.path(base_dir, paste0(CD, "_TCR-", TCR, "_-", sample3, "_-", sample4))

# フォルダが存在しない場合は作成
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

#データを抽出し保存
for (sample1 in namelist) {
  for (sample2 in namelist) {
    
    # sample1 で発現している かつ sample2 で発現している TCRクローン
    df1 <- df %>%
      filter(!is.na(.data[[sample1]]) & .data[[sample1]] > 0, 
             !is.na(.data[[sample2]]) & .data[[sample2]] > 0,
             (is.na(.data[[sample3]]) | .data[[sample3]] == 0),
             (is.na(.data[[sample4]]) | .data[[sample4]] == 0)) %>%
      select(V_J_CDR3aa,all_of(sample2))
    
    # sample1 で発現していない かつ sample2 で発現している TCRクローン
    df2 <- df %>%
      filter((is.na(.data[[sample1]]) | .data[[sample1]] == 0), 
             !is.na(.data[[sample2]]) & .data[[sample2]] > 0,
             (is.na(.data[[sample3]]) | .data[[sample3]] == 0),
             (is.na(.data[[sample4]]) | .data[[sample4]] == 0)) %>%
      select(V_J_CDR3aa,all_of(sample2))
    
    # 結果をCSVに保存
    write.csv(df1, file.path(output_folder, paste0("df1_", sample1, "_exp_--", sample2, "_exp_", TCR, ".csv")), row.names = FALSE)
    write.csv(df2, file.path(output_folder, paste0("df2_", sample1, "_non_--", sample2, "_exp_", TCR, ".csv")), row.names = FALSE)
  }
}




#-----------------------------------
#フォルダとファイルの仕分け

#file.copy(from = "source/file.csv", to = "dest/file.csv", overwrite = TRUE)
#上記のコードなどで移動させても良い。

#-----------------------------------




# グラフ化

# -----変数リストを作成-----
Tissue_list = c("Br", "TG", "Stem", "SK", "SP", "Na_TG", "Na_SP","Tet_TG_Br")
Selection_list = c("Original", "Original-Na_SP", "Original-Na_TG", "Original-Na_SP-Na_TG")
# 変数リストでループ
for (n in Tissue_list) {
  for (m in Selection_list) {
    Tissue <- n           
    Selection <- m　

    # ベースディレクトリ（保存したいフォルダのパスを指定）
    base_dir <- paste0("/Volumes/AOMINE-E1TB/基盤＿飯島先生/20250520_TCRab-seq_All_Aomine/04R/", Selection,"/same_with_",Tissue)
    # ディレクトリ移動
    setwd(base_dir)

    # PDFの出力ファイル名
    output_summarry_pdf <- paste0("same_with_", Tissue, "_",Selection, ".pdf")

    # Tissue に応じた Sample1_list と color_palette の自動選択-----
    Tissue_map <- list(
      "Br" = list(Sample1_list = c("CD4_Br_S6","CD8_Br_S1"), color_palette = "red"),
      "TG" = list(Sample1_list = c("CD4_TG_S10","CD8_TG_S2"), color_palette = "orange"),
      "Stem" = list(Sample1_list = c("CD4_Stem_S7","CD8_Stem_S3"), color_palette = "gold"),
      "SK" = list(Sample1_list = c("CD4_SK_S8","CD8_SK_S4"), color_palette = "green"),
      "SP" = list(Sample1_list = c("CD4_SP_S1","CD8_SP_S5"), color_palette = "blue"),
      "Na_TG" = list(Sample1_list = c("Na_TG_CD4_S30","Na_TG_CD8_S31"), color_palette = "salmon"),
      "Na_SP" = list(Sample1_list = c("Na_SP_CD4_S32","Na_SP_CD8_S33"), color_palette = "slateblue"),
      "Tet_TG_Br" = list(Sample1_list = c("CD4_Br_TG_S6_L001","Tet_TGandBr_CD8_S1_L001"), color_palette = "deepskyblue")
    )
    # Tissue がリストにあるか確認し、値を取得
    if (Tissue %in% names(Tissue_map)) {
      Sample1_list <- Tissue_map[[Tissue]]$Sample1_list
      color_palette <- Tissue_map[[Tissue]]$color_palette
    } else {
      stop("Error: Tissue の値が正しくありません")
    }
    # 確認用の出力
    print(paste("選択されたTissue:", Tissue))
    print(paste("Sample1_list:", paste(Sample1_list, collapse = ", ")))
    print(paste("color_palette:", color_palette))

    # TCRabおよびSample2_listの決定
    TCRab <- c("A","B")
    Sample2_list <- c("CD4_Br_S6", "CD4_Stem_S7", "CD4_SK_S8", "CD4_TG_S10", "CD4_SP_S1", "Na_TG_CD4_S30", 
                      "Na_SP_CD4_S32", "CD4_Br_TG_S6_L001", "CD8_Br_S1", "CD8_Stem_S3", "CD8_SK_S4", "CD8_TG_S2", "CD8_SP_S5", 
                      "Na_TG_CD8_S31", "Na_SP_CD8_S33","Tet_TGandBr_CD8_S1_L001")

    # グラフを保存するリストを作成
    plot_list <- list()

    # グラフを作成するループ
    for (sample1 in Sample1_list) {
      for (TCR in TCRab) {
        for (sample2 in Sample2_list) {
      
          print(paste("Processing: sample1 =", sample1, ", sample2 =", sample2, ", TCR =", TCR))
      
          # 2つのCSVファイルのパスを指定
          file1 <- paste0("df1_", sample1, "_exp_--", sample2, "_exp_", TCR, ".csv")
          file2 <- paste0("df2_", sample1, "_non_--", sample2, "_exp_", TCR, ".csv")

          # **ファイルが両方とも存在するかチェック**
          if (!(file.exists(file1) & file.exists(file2))) {
            print(paste("ファイルが見つからない: ", file1, file2))
            next  # どちらかのファイルが存在しなければスキップ
          }
      
　　    　# 名前を抽出
          sample1_clean <- sub("CD[48]_", "", sample1)
          sample2_clean <- sub("_[^_]*$", "", sample2)

          # ファイル名からタイトルを生成
          file_title <- paste0("c=", sample1_clean, "_", sample2_clean, "_TCR", TCR)
          labs_title <- paste0(sample2_clean, "_TCR", TCR)

          # 2つのCSVファイルを読み込み
          df1 <- read.csv(file1, header = TRUE)
          df2 <- read.csv(file2, header = TRUE)
      
          # カラーマッピングの辞書を作成
          color_map <- list(
            "red" = c("red4", "red3", "red"),
            "orange" = c("orange3", "orange2", "orange"),
            "gold" = c("gold3", "gold2", "gold1"),
            "green" = c("green4", "green3", "green2"),
            "blue" = c("blue4", "blue3", "blue2"),
            "salmon" = c("salmon3", "salmon2", "salmon"),
            "slateblue" = c("slateblue4", "slateblue", "slateblue1"),
            "deepskyblue" = c("deepskyblue", "deepskyblue3", "deepskyblue4")
          )
      
          # どちらかが空の場合、片方の df のみを使用
          print(paste("df1:", nrow(df1), "df2:", nrow(df2)))
          if (nrow(df1) == 0 & nrow(df2) == 0) {
            next  # 両方とも空ならスキップ
          } else if (nrow(df1) == 0) {
            colnames(df2)[2] <- "percentage"
            df2 <- df2 %>%
              mutate(color = rep(c("grey51", "grey77", "grey88"), length.out = n()))
            df_combined <- df2  # df2 のみ使用
          } else if (nrow(df2) == 0) {
            colnames(df1)[2] <- "percentage"
            df1 <- df1 %>%
              mutate(color = rep(color_map[[color_palette]], length.out = n()))
            df_combined <- df1  # df1 のみ使用
          } else {
      
          　# 2列目の名前を 'percentage' に変更
          　colnames(df1)[2] <- "percentage"
          　colnames(df2)[2] <- "percentage"

          　# "percentage" 列を数値の大きい順に並べ替え
          　df1 <- df1 %>% arrange(desc(percentage))
          　df2 <- df2 %>% arrange(desc(percentage))

          　# 'df'列を追加して、どのデータから来ているかを区別
          　df1$df <- "df1"
          　df2$df <- "df2"
      
            # 型を統一（重要）
            df1$V_J_CDR3aa <- as.character(df1$V_J_CDR3aa)
            df2$V_J_CDR3aa <- as.character(df2$V_J_CDR3aa)
      
          　# 'color'列を追加して、色を指定
          　df1 <- df1 %>%
            　mutate(color = rep(color_map[[color_palette]], length.out = n()))
          　df2 <- df2 %>%
            　mutate(color = rep(c("grey51", "grey77", "grey88"), length.out = n()))
      
            # 両方ある場合は結合
            df_combined <- bind_rows(df1, df2)
          }
      
          # "percentage" 列を数値の大きい順に並べ替え
          df_combined <- df_combined %>% arrange(desc(percentage))

          # パーセンテージを計算(割合になっている場合はスキップ)
          df_combined <- df_combined %>%
          　mutate(percentage = percentage / sum(percentage, na.rm = TRUE) * 100)

          # 'color' 列を factor 型もしくはcharacter型に変換
          # もし color 列が "#FF0000FF" や "rgb(1,0,0)" の形式の場合は、factor には変換せずそのまま character 型で使用する
          #df_combined$color <- as.factor(df_combined$color) #to factor
          #df_combined$color <- as.character(df_combined$color) #to character
          #sapply(df_combined, class) #確認

          # カラーマッピングの辞書を作成
          color_map2 <- list(
            "red" = c("red4" = "red4", "red3" = "red3", "red" = "red", "grey51"= "grey51", "grey77" = "grey77", "grey88" = "grey88"),
            "orange" = c("orange3" = "orange3", "orange2" = "orange2", "orange" = "orange", "grey51"= "grey51", "grey77" = "grey77", "grey88" = "grey88"),
            "gold" = c("gold3" = "gold3", "gold2" = "gold2", "gold1" = "gold1", "grey51"= "grey51", "grey77" = "grey77", "grey88" = "grey88"),
            "green" = c("green4" = "green4", "green3" = "green3", "green2" = "green2", "grey51"= "grey51", "grey77" = "grey77", "grey88" = "grey88"),
            "blue" = c("blue4" = "blue4", "blue3" = "blue3", "blue2" = "blue2", "grey51"= "grey51", "grey77" = "grey77", "grey88" = "grey88"),
            "salmon" = c("salmon3" = "salmon3", "salmon2" = "salmon2", "salmon" = "salmon", "grey51"= "grey51", "grey77" = "grey77", "grey88" = "grey88"),
            "slateblue" = c("slateblue4" = "slateblue4", "slateblue" = "slateblue", "slateblue1" = "slateblue1", "grey51"= "grey51", "grey77" = "grey77", "grey88" = "grey88"),
            "deepskyblue" = c("deepskyblue" = "deepskyblue", "deepskyblue3" = "deepskyblue3", "deepskyblue4" = "deepskyblue4", "grey51"= "grey51", "grey77" = "grey77", "grey88" = "grey88")
          )


          # ドーナツチャートの作成
          donut <- ggplot(df_combined, aes(x = 2, y = percentage, fill = color, group = percentage)) +
                     geom_bar(stat = "identity", color = NA, linewidth = 0.1, width = 1) + #color = NAの所でセル間の線の色を指定
                     coord_polar("y", start = 0) +
                     theme_void() +  # 背景や軸を消す
                     xlim(0.5, 2.5) +  # ドーナツにするために中心にスペースを作る
                     theme(legend.position = "none") +  # 凡例を非表示にする
                     #geom_text(aes(label = VJ.CDR3aa), position = position_stack(vjust = 0.5), size = 3) +  # ラベルを追加
                     labs(title = labs_title) + # タイトルをファイル名から自動的に設定
                     scale_fill_manual(values = color_map2[[color_palette]])

          # 保存するPDFファイル名を生成
          output_file <- paste0(file_title, ".pdf")
          # グラフをPDFに保存する
          ggsave(output_file, plot = donut, width = 8, height = 8)
          # 作成したグラフをリストに追加
          plot_list <- append(plot_list, list(donut))
        }
      }
    }

    # PDFに出力
    if (length(plot_list) > 0) {
      pdf(output_summarry_pdf, width = 20, height = 10)  # 幅と高さを調整
      do.call(grid.arrange, c(plot_list, nrow = 4))  # 2行に並べて表示
      dev.off()
      print("PDFファイルが作成されました！")
    } else {
      print("エラー: グラフが作成されていません。PDFは保存されません。")
    }

  }
}




