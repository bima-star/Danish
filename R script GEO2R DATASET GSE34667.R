# ================================================================
# GEO2R Analysis - GSE34667 
# ================================================================

library(GEOquery)
library(limma)
library(umap)

# Load data
gset <- getGEO("GSE34667", GSEMatrix = TRUE, AnnotGPL = TRUE)
if (length(gset) > 1) idx <- grep("GPL198", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

fvarLabels(gset) <- make.names(fvarLabels(gset))

# Group membership
gsms <- "111000111000111100001110000"
sml <- strsplit(gsms, split="")[[1]]

# Log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { 
  ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) 
}

# Design matrix
gs <- factor(sml)
groups <- make.names(c("group 1","group 2"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ]

# Fit model
fit <- lmFit(gset, design)
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

# Get results
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

# =============================================
# FIX: Pilih kolom yang ADA saja
# =============================================
cat("Kolom yang tersedia:\n")
print(colnames(tT))

wanted_cols <- c("ID", "adj.P.Val", "P.Value", "t", "B", "logFC", 
                 "ORF", "SPOT_ID", "Gene.Symbol", "Gene.symbol", "Gene.title")
available_cols <- wanted_cols[wanted_cols %in% colnames(tT)]

cat("\nKolom yang diambil:\n")
print(available_cols)

tT <- tT[, available_cols]

write.table(tT, file="DEG_results.tsv", row.names=FALSE, sep="\t")
cat("\nHasil tersimpan di: DEG_results.tsv\n")

# =============================================
# VISUALISASI
# =============================================

# Full results for visualization
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

# 1. Histogram P-value
png("01_histogram_pvalue.png", width=800, height=600)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")
dev.off()
cat("✓ 01_histogram saved\n")

# 2. Decide tests
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=0)

# 3. Venn diagram
png("02_venn_diagram.png", width=700, height=600)
vennDiagram(dT, circle.col=palette())
dev.off()
cat("✓ 02_venn_diagram saved\n")

# 4. QQ plot
png("03_QQ_plot.png", width=700, height=600)
t.good <- which(!is.na(fit2$F))
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")
dev.off()
cat("✓ 03_QQ_plot saved\n")

# 5. Volcano plot
ct <- 1
png("04_volcano_plot.png", width=800, height=600)
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))
dev.off()
cat("✓ 04_volcano_plot saved\n")

# 6. MD plot
png("05_MD_plot.png", width=800, height=600)
plotMD(fit2, column=ct, status=dT[,ct], legend=FALSE, pch=20, cex=1)
abline(h=0)
dev.off()
cat("✓ 05_MD_plot saved\n")

# 7. Box plot
ex <- exprs(gset)
png("06_boxplot.png", width=900, height=600)
ord <- order(gs)
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste("GSE34667", "/", annotation(gset), sep="")
boxplot(ex[,ord], boxwex=0.6, notch=TRUE, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")
dev.off()
cat("✓ 06_boxplot saved\n")

# 8. Density plot
png("07_density_plot.png", width=800, height=600)
par(mar=c(4,4,2,1))
title <- paste("GSE34667", "/", annotation(gset), " value distribution", sep="")
plotDensities(ex, group=gs, main=title, legend="topright")
dev.off()
cat("✓ 07_density_plot saved\n")

# 9. UMAP plot (FIXED - tanpa maptools)
png("08_UMAP_plot.png", width=800, height=600)
ex_clean <- na.omit(ex)
ex_clean <- ex_clean[!duplicated(ex_clean), ]
ump <- umap(t(ex_clean), n_neighbors = 11, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=11", xlab="", ylab="", 
     col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
text(ump$layout, labels = rownames(ump$layout), cex=0.5, pos=3)
dev.off()
cat("✓ 08_UMAP_plot saved\n")

# 10. Mean-variance trend (FIXED)
png("09_mean_variance.png", width=800, height=600)
tryCatch({
  limma::plotSA(fit2, main="Mean variance trend, GSE34667")
}, error = function(e) {
  plot(fit2$Amean, sqrt(fit2$sigma),
       pch = 16, cex = 0.3, col = "darkgrey",
       main = "Mean-Variance Trend, GSE34667",
       xlab = "Average Log Expression",
       ylab = "sqrt(sigma)")
  lines(lowess(fit2$Amean, sqrt(fit2$sigma)), col = "red", lwd = 2)
})
dev.off()
cat("✓ 09_mean_variance saved\n")

# =============================================
# SELESAI!
# =============================================
cat("\n====================================\n")
cat("Semua visualisasi berhasil disimpan!\n")
cat("====================================\n")
cat("Lokasi file:", getwd(), "\n")
