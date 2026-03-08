# ================================================================
# Modul: Analisis Ekspresi Gen Arabidopsis thaliana
# Dataset: GSE34667
# Platform: Microarray Affymetrix ATH1 Genome Array (GPL198)
# Kondisi: Clean Air vs Ozone Treatment
# Tujuan: Mengidentifikasi Differentially Expressed Genes (DEG)
# ================================================================

# ================================================================
# PART A. PENGANTAR KONSEP
# ================================================================

# Analisis ekspresi gen bertujuan untuk membandingkan tingkat ekspresi gen
# antara dua kondisi biologis (Clean Air vs Ozone treatment)
# Pada modul ini kita menggunakan pendekatan statistik limma (Linear Models
# for Microarray Data), yang merupakan standar emas untuk data microarray.

# ================================================================
# PART B. PERSIAPAN LINGKUNGAN KERJA (INSTALL & LOAD PACKAGE)
# ================================================================

# Apa itu package?
# Package adalah kumpulan fungsi siap pakai di R
# Bioinformatika di R sangat bergantung pada package dari CRAN dan Bioconductor

# 1. Install BiocManager (manajer paket Bioconductor)
# IF adalah struktur logika : "jika kondisi terpenuhi, lakukan aksi"

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 2. Install paket Bioconductor (GEOquery & limma)
# GEOquery: mengambil data dari database GEO
# limma: analisis statistik ekspresi gen

BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update = FALSE)

# ath1121501.db: database anotasi khusus untuk chip Affymetrix ATH1 (Arabidopsis)
BiocManager::install("ath1121501.db", ask = FALSE, update = FALSE)

# 3. Install paket CRAN untuk visualisasi dan manipulasi data
# pheatmap: heatmap ekspresi gen
# ggplot2: grafik (volcano plot)
# dplyr: manipulasi tabel data

install.packages(c("pheatmap", "ggplot2", "dplyr"))

# umap: grafik (plot UMAP)
if (!requireNamespace("umap", quietly = TRUE)) {
  install.packages("umap")
}

# 4. Memanggil library
# library() digunakan agar fungsi di dalam package bisa digunakan

library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(ath1121501.db)
library(AnnotationDbi)
library(umap)

# ================================================================
# PART C. PENGAMBILAN DATA DARI GEO
# ================================================================

# GEO (Gene Expression Omnibus) adalah database publik milik NCBI
# getGEO(): fungsi untuk mengunduh dataset berdasarkan ID GEO
# GSEMatrix = TRUE -> data diambil dalam format ExpressionSet
# AnnotGPL  = TRUE -> anotasi gen (Gene Symbol) ikut diunduh

gset <- getGEO("GSE34667", GSEMatrix = TRUE, AnnotGPL = TRUE)

# Pilih platform yang sesuai (GPL198 = Affymetrix ATH1)
if (length(gset) > 1) {
  idx <- grep("GPL198", attr(gset, "names"))
} else {
  idx <- 1
}
gset <- gset[[idx]]

# ExpressionSet berisi:
# - exprs() : matriks ekspresi gen
# - pData() : metadata sampel
# - fData() : metadata fitur (probe / gen)

# Cek info dataset
cat("================================================\n")
cat("INFO DATASET GSE34667\n")
cat("================================================\n")
cat("Jumlah probe:", nrow(gset), "\n")
cat("Jumlah sampel:", ncol(gset), "\n")
cat("Platform:", annotation(gset), "\n")
cat("================================================\n")

# ================================================================
# PART D. PRE-PROCESSING DATA EKSPRESI
# ================================================================

# exprs(): mengambil matriks ekspresi gen
# Baris  = probe/gen
# Kolom  = sampel
ex <- exprs(gset)

# Mengapa perlu log2 transformasi?
# Data microarray mentah memiliki rentang nilai sangat besar.
# Log2 digunakan untuk:
# 1. Menstabilkan varians
# 2. Mendekati asumsi model linear
# 3. Memudahkan interpretasi log fold change

# quantile(): menghitung nilai kuantil (persentil)
# as.numeric(): mengubah hasil quantile (yang berupa named vector)
# menjadi vektor numerik biasa agar mudah dibandingkan
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))

# LogTransform adalah variabel logika (TRUE / FALSE)
# Operator logika:
# >  : lebih besar dari
# || : OR (atau)
# && : AND (dan)
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

# IF statement:
# Jika LogTransform = TRUE, maka lakukan log2
if (LogTransform) {
  # Nilai <= 0 tidak boleh di-log, maka diubah menjadi NA
  ex[ex <= 0] <- NA
  ex <- log2(ex)
  cat("✓ Data sudah di-log2 transform\n")
} else {
  cat("✓ Data tidak perlu log2 transform\n")
}

# Update expression data di gset
exprs(gset) <- ex

# Hapus baris dengan missing values
gset <- gset[complete.cases(exprs(gset)), ]
ex <- exprs(gset)

cat("✓ Jumlah probe setelah filter:", nrow(gset), "\n")

# ================================================================
# PART E. DEFINISI KELOMPOK SAMPEL
# ================================================================

# Untuk GSE34667, kita menggunakan kode grup dari GEO2R
# "1" = Clean Air (treatment udara bersih)
# "0" = Ozone (treatment ozon)

gsms <- "111000111000111100001110000"
sml <- strsplit(gsms, split = "")[[1]]

# Verifikasi jumlah sampel
cat("\n================================================\n")
cat("DEFINISI KELOMPOK SAMPEL\n")
cat("================================================\n")
cat("Jumlah kode grup:", length(sml), "\n")
cat("Jumlah sampel di dataset:", ncol(gset), "\n")

if (length(sml) != ncol(gset)) {
  stop("ERROR: Jumlah kode grup tidak sama dengan jumlah sampel!")
}

# Buat faktor
gs <- factor(sml)

# Cek urutan level
cat("\nUrutan level faktor:\n")
print(levels(gs))

# Beri nama sesuai urutan level:
# "0" (level pertama) = Ozone
# "1" (level kedua) = Clean_Air
groups <- make.names(c("Ozone", "Clean_Air"))
levels(gs) <- groups

# Tambahkan ke dataset
gset$group <- gs

# Verifikasi hasil
cat("\nNama grup yang terbentuk:\n")
print(levels(gset$group))

cat("\nDistribusi sampel per grup:\n")
print(table(gset$group))

# Lihat detail per sampel
sample_info <- data.frame(
  No = 1:ncol(gset),
  Sample = sampleNames(gset),
  Group = gset$group
)
cat("\nDetail sampel:\n")
print(sample_info)

# Simpan nama grup untuk dipakai nanti
nama_grup <- levels(gset$group)
cat("\nNama grup: ", nama_grup[1], "dan", nama_grup[2], "\n")
cat("================================================\n")

# ================================================================
# PART F. DESIGN MATRIX (KERANGKA STATISTIK)
# ================================================================

# model.matrix():
# Membuat matriks desain untuk model linear
# ~0 berarti TANPA intercept (best practice limma)
design <- model.matrix(~0 + gset$group)

# colnames(): memberi nama kolom agar mudah dibaca
colnames(design) <- levels(gset$group)

# Lihat design matrix
cat("\nDesign matrix:\n")
print(design)

# Menentukan perbandingan biologis
# Bandingkan: Clean_Air vs Ozone
grup_treatment <- "Clean_Air"
grup_control <- "Ozone"

contrast_formula <- paste(grup_treatment, "-", grup_control)
cat("\n================================================\n")
cat("Kontras yang dianalisis:", contrast_formula, "\n")
cat("Interpretasi:\n")
cat("  - logFC > 0 = gen lebih tinggi di Clean Air\n")
cat("  - logFC < 0 = gen lebih tinggi di Ozone\n")
cat("================================================\n")

# ================================================================
# PART G. ANALISIS DIFFERENTIAL EXPRESSION (LIMMA)
# ================================================================

# lmFit():
# Membangun model linear untuk setiap gen
fit <- lmFit(ex, design)

# makeContrasts(): mendefinisikan perbandingan antar grup
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)

# contrasts.fit(): menerapkan kontras ke model
fit2 <- contrasts.fit(fit, contrast_matrix)

# eBayes():
# Empirical Bayes untuk menstabilkan estimasi varians
fit2 <- eBayes(fit2)

# topTable():
# Mengambil hasil akhir DEG
# adjust = "fdr" -> koreksi multiple testing
# p.value = 0.05 -> threshold signifikansi
topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.05
)

cat("\n================================================\n")
cat("HASIL ANALISIS DEG\n")
cat("================================================\n")
cat("Jumlah DEG signifikan (adj.P.Val < 0.05):", nrow(topTableResults), "\n")
cat("\nTop 10 DEG:\n")
print(head(topTableResults, 10))
cat("================================================\n")

# Ambil semua hasil untuk visualisasi
allResults <- topTable(fit2, adjust = "fdr", sort.by = "B", number = Inf)

# ================================================================
# PART H. ANOTASI NAMA GEN
# ================================================================

# Penting:
# Pada data microarray Affymetrix, unit analisis awal adalah PROBE ID,
# bukan gen. Oleh karena itu, anotasi ulang diperlukan menggunakan
# database resmi Bioconductor (ath1121501.db).

cat("\n================================================\n")
cat("ANOTASI GEN\n")
cat("================================================\n")

# Mengambil ID probe dari hasil DEG
probe_ids <- rownames(topTableResults)

# Lihat kolom yang tersedia di database anotasi
cat("Kolom anotasi yang tersedia:\n")
print(columns(ath1121501.db))

# Mapping probe -> gene symbol & gene name
# Gunakan database ath1121501.db untuk Arabidopsis
gene_annotation <- tryCatch({
  AnnotationDbi::select(
    ath1121501.db,
    keys = probe_ids,
    columns = c("SYMBOL", "GENENAME", "TAIR", "ENTREZID"),
    keytype = "PROBEID"
  )
}, error = function(e) {
  # Jika error, coba dengan kolom yang lebih sedikit
  cat("Warning: Menggunakan kolom anotasi terbatas\n")
  AnnotationDbi::select(
    ath1121501.db,
    keys = probe_ids,
    columns = c("SYMBOL", "TAIR"),
    keytype = "PROBEID"
  )
})

# Gabungkan dengan hasil limma
topTableResults$PROBEID <- rownames(topTableResults)

topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)

# Hapus duplikat (jika ada probe yang map ke multiple genes)
topTableResults <- topTableResults[!duplicated(topTableResults$PROBEID), ]

# Cek hasil anotasi
cat("\nContoh hasil anotasi:\n")
print(head(topTableResults[, c("PROBEID", "SYMBOL", "TAIR", "logFC", "adj.P.Val")]))
cat("================================================\n")

# ================================================================
# PART I.1 BOXPLOT DISTRIBUSI NILAI EKSPRESI
# ================================================================

# Boxplot digunakan untuk:
# - Mengecek distribusi nilai ekspresi antar sampel
# - Melihat apakah ada batch effect
# - Mengevaluasi apakah normalisasi/log-transform sudah wajar

cat("\n================================================\n")
cat("MEMBUAT VISUALISASI...\n")
cat("================================================\n")

# Set warna berdasarkan grup
# Hijau untuk Clean Air, Ungu untuk Ozone
group_colors <- ifelse(gset$group == "Clean_Air", "#4DAF4A", "#984EA3")

png("01_boxplot_distribusi.png", width = 1000, height = 600)
par(mar = c(10, 4, 4, 2))
boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Boxplot Distribusi Nilai Ekspresi per Sampel\nGSE34667 - Arabidopsis thaliana",
  ylab = "Expression Value (log2)",
  cex.axis = 0.7
)
legend(
  "topright",
  legend = c("Clean Air", "Ozone"),
  fill = c("#4DAF4A", "#984EA3"),
  cex = 0.8
)
dev.off()
cat("✓ 01_boxplot_distribusi.png tersimpan\n")

# ================================================================
# PART I.2 DISTRIBUSI NILAI EKSPRESI (DENSITY PLOT)
# ================================================================

# Density plot menunjukkan sebaran global nilai ekspresi gen
# Digunakan untuk:
# - Mengecek efek log-transform
# - Membandingkan distribusi antar grup

# Gabungkan ekspresi & grup ke data frame
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)

png("02_density_plot.png", width = 800, height = 600)
print(
  ggplot(expr_long, aes(x = Expression, color = Group)) +
    geom_density(linewidth = 1) +
    scale_color_manual(
      values = c("Clean_Air" = "#4DAF4A", "Ozone" = "#984EA3"),
      labels = c("Clean Air", "Ozone")
    ) +
    theme_minimal() +
    labs(
      title = "Distribusi Nilai Ekspresi Gen",
      subtitle = "GSE34667 - Arabidopsis thaliana (Clean Air vs Ozone)",
      x = "Expression Value (log2)",
      y = "Density"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
)
dev.off()
cat("✓ 02_density_plot.png tersimpan\n")

# ================================================================
# PART I.3 UMAP (VISUALISASI DIMENSI RENDAH)
# ================================================================

# UMAP digunakan untuk:
# - Mereduksi ribuan gen menjadi 2 dimensi
# - Melihat pemisahan sampel secara global

# Hapus NA dari matriks ekspresi
ex_clean <- na.omit(ex)

# Transpose matriks ekspresi:
# UMAP bekerja pada OBSERVATION = sampel
umap_input <- t(ex_clean)

# Jalankan UMAP
set.seed(123) # untuk reproducibility
umap_result <- umap(umap_input, n_neighbors = min(15, ncol(ex) - 1))

# Simpan hasil ke data frame
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group,
  Sample = sampleNames(gset)
)

# Plot UMAP
png("03_UMAP_plot.png", width = 800, height = 600)
print(
  ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
    geom_point(size = 4, alpha = 0.8) +
    scale_color_manual(
      values = c("Clean_Air" = "#4DAF4A", "Ozone" = "#984EA3"),
      labels = c("Clean Air", "Ozone")
    ) +
    theme_minimal() +
    labs(
      title = "UMAP Plot: Clean Air vs Ozone",
      subtitle = "GSE34667 - Arabidopsis thaliana",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
)
dev.off()
cat("✓ 03_UMAP_plot.png tersimpan\n")

# ================================================================
# PART I.4 PCA PLOT
# ================================================================

# PCA untuk melihat pengelompokan sampel
pca <- prcomp(t(ex_clean), scale. = TRUE)
var_explained <- round(100 * summary(pca)$importance[2, 1:2], 1)

pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Group = gset$group,
  Sample = sampleNames(gset)
)

png("04_PCA_plot.png", width = 800, height = 600)
print(
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 4, alpha = 0.8) +
    stat_ellipse(level = 0.95, linetype = "dashed") +
    scale_color_manual(
      values = c("Clean_Air" = "#4DAF4A", "Ozone" = "#984EA3"),
      labels = c("Clean Air", "Ozone")
    ) +
    theme_minimal() +
    labs(
      title = "PCA Plot: Clean Air vs Ozone",
      subtitle = "GSE34667 - Arabidopsis thaliana",
      x = paste0("PC1 (", var_explained[1], "%)"),
      y = paste0("PC2 (", var_explained[2], "%)")
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
)
dev.off()
cat("✓ 04_PCA_plot.png tersimpan\n")

# ================================================================
# PART J.1 VISUALISASI VOLCANO PLOT
# ================================================================

# Volcano plot menggabungkan:
# - Log fold change (efek biologis)
# - Signifikansi statistik

# Anotasi untuk semua gen (untuk volcano plot)
all_probe_ids <- rownames(allResults)

all_gene_annotation <- tryCatch({
  AnnotationDbi::select(
    ath1121501.db,
    keys = all_probe_ids,
    columns = c("SYMBOL", "TAIR"),
    keytype = "PROBEID"
  )
}, error = function(e) {
  data.frame(PROBEID = all_probe_ids, SYMBOL = NA, TAIR = NA)
})

volcano_data <- data.frame(
  PROBEID = rownames(allResults),
  logFC = allResults$logFC,
  adj.P.Val = allResults$adj.P.Val
)

# Tambahkan gene symbol
volcano_data <- merge(
  volcano_data,
  all_gene_annotation[, c("PROBEID", "SYMBOL")],
  by = "PROBEID",
  all.x = TRUE
)
volcano_data <- volcano_data[!duplicated(volcano_data$PROBEID), ]

# Klasifikasi status gen
volcano_data$status <- "Not Significant"
volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val < 0.05] <- "Up in Clean Air"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.05] <- "Up in Ozone"

# Hitung jumlah
cat("\nJumlah gen per kategori:\n")
print(table(volcano_data$status))

# Visualisasi
png("05_volcano_plot.png", width = 900, height = 700)
print(
  ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c(
      "Up in Ozone" = "#984EA3",
      "Not Significant" = "grey70",
      "Up in Clean Air" = "#4DAF4A"
    )) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    theme_minimal() +
    labs(
      title = "Volcano Plot - GSE34667 (Arabidopsis thaliana)",
      subtitle = "Clean Air vs Ozone Treatment",
      x = "Log2 Fold Change (Clean Air / Ozone)",
      y = "-Log10 Adjusted P-value",
      color = "Expression"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    annotate("text", x = 2.5, y = max(-log10(volcano_data$adj.P.Val), na.rm = TRUE) * 0.9,
             label = "↑ Higher in\nClean Air", color = "#4DAF4A", size = 3.5) +
    annotate("text", x = -2.5, y = max(-log10(volcano_data$adj.P.Val), na.rm = TRUE) * 0.9,
             label = "↑ Higher in\nOzone", color = "#984EA3", size = 3.5)
)
dev.off()
cat("✓ 05_volcano_plot.png tersimpan\n")

# ================================================================
# PART J.2 VISUALISASI HEATMAP
# ================================================================

# Heatmap digunakan untuk melihat pola ekspresi gen
# antar sampel berdasarkan gen-gen paling signifikan

# Pilih 50 gen paling signifikan berdasarkan adj.P.Val
topTableResults <- topTableResults[order(topTableResults$adj.P.Val), ]
top50 <- head(topTableResults, 50)

# Ambil matriks ekspresi untuk gen terpilih
mat_heatmap <- ex[top50$PROBEID, ]

# Gunakan Gene Symbol (fallback ke TAIR ID atau Probe ID jika kosong)
gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "",
  ifelse(is.na(top50$TAIR) | top50$TAIR == "",
         top50$PROBEID,
         top50$TAIR),
  top50$SYMBOL
)
rownames(mat_heatmap) <- make.unique(gene_label)

# Pembersihan data (WAJIB agar tidak error hclust)
# Hapus baris dengan NA
mat_heatmap <- mat_heatmap[rowSums(is.na(mat_heatmap)) == 0, ]

# Hapus gen dengan varians nol
gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]

# Anotasi kolom (kelompok sampel)
annotation_col <- data.frame(
  Group = gset$group
)
rownames(annotation_col) <- colnames(mat_heatmap)

# Warna anotasi
annotation_colors <- list(
  Group = c("Clean_Air" = "#4DAF4A", "Ozone" = "#984EA3")
)

# Visualisasi heatmap
png("06_heatmap_top50.png", width = 900, height = 1200)
pheatmap(
  mat_heatmap,
  scale = "row",
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 8,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  main = "Top 50 DEGs - GSE34667\nClean Air vs Ozone (Arabidopsis thaliana)"
)
dev.off()
cat("✓ 06_heatmap_top50.png tersimpan\n")

# ================================================================
# PART J.3 BARPLOT JUMLAH DEG
# ================================================================

# Hitung Up dan Down regulated
up_clean_air <- sum(volcano_data$status == "Up in Clean Air", na.rm = TRUE)
up_ozone <- sum(volcano_data$status == "Up in Ozone", na.rm = TRUE)

deg_summary <- data.frame(
  Direction = c("Higher in Clean Air", "Higher in Ozone"),
  Count = c(up_clean_air, up_ozone)
)

png("07_barplot_DEG.png", width = 700, height = 500)
print(
  ggplot(deg_summary, aes(x = Direction, y = Count, fill = Direction)) +
    geom_bar(stat = "identity", width = 0.5) +
    geom_text(aes(label = Count), vjust = -0.5, size = 6) +
    scale_fill_manual(values = c(
      "Higher in Clean Air" = "#4DAF4A",
      "Higher in Ozone" = "#984EA3"
    )) +
    theme_minimal() +
    labs(
      title = "Jumlah Differentially Expressed Genes",
      subtitle = "GSE34667 - Clean Air vs Ozone (|logFC| > 1, adj.P.Val < 0.05)",
      x = "",
      y = "Jumlah Gen"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none"
    ) +
    ylim(0, max(deg_summary$Count) * 1.2)
)
dev.off()
cat("✓ 07_barplot_DEG.png tersimpan\n")

# ================================================================
# PART J.4 HISTOGRAM P-VALUE
# ================================================================

png("08_histogram_pvalue.png", width = 800, height = 600)
hist(
  allResults$adj.P.Val,
  col = "grey70",
  border = "white",
  breaks = 50,
  xlab = "Adjusted P-value",
  ylab = "Number of Genes",
  main = "Distribution of Adjusted P-values\nGSE34667 - Clean Air vs Ozone"
)
abline(v = 0.05, col = "red", lty = 2, lwd = 2)
legend("topright", legend = "p = 0.05", col = "red", lty = 2, lwd = 2)
dev.off()
cat("✓ 08_histogram_pvalue.png tersimpan\n")

# ================================================================
# PART J.5 MA PLOT
# ================================================================

png("09_MA_plot.png", width = 800, height = 600)
allResults$Significance <- ifelse(allResults$adj.P.Val < 0.05, "Significant", "Not Significant")

print(
  ggplot(allResults, aes(x = AveExpr, y = logFC, color = Significance)) +
    geom_point(size = 1, alpha = 0.5) +
    scale_color_manual(values = c("Significant" = "#E41A1C", "Not Significant" = "grey70")) +
    geom_hline(yintercept = c(-1, 0, 1),
               linetype = c("dashed", "solid", "dashed"),
               color = c("#984EA3", "black", "#4DAF4A")) +
    theme_minimal() +
    labs(
      title = "MA Plot - GSE34667",
      subtitle = "Clean Air vs Ozone (Arabidopsis thaliana)",
      x = "Average Expression (log2)",
      y = "Log2 Fold Change"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
)
dev.off()
cat("✓ 09_MA_plot.png tersimpan\n")

# ================================================================
# PART K. MENYIMPAN HASIL
# ================================================================

# write.csv(): menyimpan hasil analisis ke file CSV
write.csv(topTableResults, "Hasil_GSE34667_DEG_signifikan.csv", row.names = FALSE)
cat("✓ Hasil_GSE34667_DEG_signifikan.csv tersimpan\n")

# Simpan semua hasil (termasuk non-signifikan)
allResults$PROBEID <- rownames(allResults)
allResults_annotated <- merge(
  allResults,
  all_gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)
allResults_annotated <- allResults_annotated[!duplicated(allResults_annotated$PROBEID), ]
write.csv(allResults_annotated, "Hasil_GSE34667_ALL.csv", row.names = FALSE)
cat("✓ Hasil_GSE34667_ALL.csv tersimpan\n")

# ================================================================
# RINGKASAN HASIL
# ================================================================

cat("\n")
cat("════════════════════════════════════════════════════════\n")
cat("           RINGKASAN ANALISIS GSE34667                  \n")
cat("════════════════════════════════════════════════════════\n")
cat("Dataset        : GSE34667\n")
cat("Organisme      : Arabidopsis thaliana\n")
cat("Platform       : Affymetrix ATH1 Genome Array (GPL198)\n")
cat("Kondisi        : Clean Air vs Ozone Treatment\n")
cat("────────────────────────────────────────────────────────\n")
cat("Jumlah sampel  : ", ncol(gset), "\n")
cat("  - Clean Air  : ", sum(gset$group == "Clean_Air"), " sampel\n")
cat("  - Ozone      : ", sum(gset$group == "Ozone"), " sampel\n")
cat("────────────────────────────────────────────────────────\n")
cat("Jumlah probe   : ", nrow(gset), "\n")
cat("DEG signifikan : ", nrow(topTableResults), " (adj.P.Val < 0.05)\n")
cat("  - Higher in Clean Air : ", up_clean_air, " gen (logFC > 1)\n")
cat("  - Higher in Ozone     : ", up_ozone, " gen (logFC < -1)\n")
cat("════════════════════════════════════════════════════════\n")
cat("\nFile output tersimpan di:\n")
cat(getwd(), "\n")
cat("\nDaftar file yang dihasilkan:\n")
cat("  📊 01_boxplot_distribusi.png\n")
cat("  📊 02_density_plot.png\n")
cat("  📊 03_UMAP_plot.png\n")
cat("  📊 04_PCA_plot.png\n")
cat("  📊 05_volcano_plot.png\n")
cat("  📊 06_heatmap_top50.png\n")
cat("  📊 07_barplot_DEG.png\n")
cat("  📊 08_histogram_pvalue.png\n")
cat("  📊 09_MA_plot.png\n")
cat("  📄 Hasil_GSE34667_DEG_signifikan.csv\n")
cat("  📄 Hasil_GSE34667_ALL.csv\n")
cat("════════════════════════════════════════════════════════\n")

getwd()

message("\n✓ Analisis selesai! Semua file telah disimpan.")