# ------------------------------
# 1. Data Collection:
# ------------------------------
library(TCGAbiolinks)

query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")

GDCdownload(query)

# ------------------------------
# 2. Data Preparation:
# ------------------------------
prep <- GDCprepare(query)

# ------------------------------
# 3. Preprocessing:
# ------------------------------
library(TCGAbiolinks)

data_prepro <- TCGAanalyze_Preprocessing(prep, cor.cut = 0.6)

# ------------------------------
# 4. Normalization:
# ------------------------------
library(TCGAbiolinks)

dataNorm <- TCGAanalyze_Normalization(data_prepro, geneInfo = geneInfoHT, method = "geneLength")

# ------------------------------
# 5. Filtering:
# ------------------------------
library(TCGAbiolinks)

dataFilt <- TCGAanalyze_Filtering(dataNorm, method = "quantile", qnt.cut = 0.25)

# ------------------------------
# 6. Voom Transformation:
# ------------------------------
library(limma)

v <- voom(dataFilt, plot = TRUE)
E <- v$E
ex <- 2^E
mrnaEx <- log2(ex + 1)

# ------------------------------
# 7. Sample Annotation:
# ------------------------------
samples <- data.frame(colnames(mrnaEx))
colnames(samples) <- "barcode"
samples$substr <- substr(samples$barcode, 14, 15)
samples$Type[samples$substr > 10] <- "Normal"
samples$Type[samples$substr < 9] <- "Tumor"

write.table(samples, "D:/my freelancing/Tahmine/tcga/samples.txt", quote = FALSE, sep = "\t")

# ------------------------------
# 8. Group Selection:
# ------------------------------
library(TCGAbiolinks)

bar <- samples  
n <- bar[bar$substr == 11, ]
t <- bar[bar$substr == "01", ]

Tsample <- TCGAquery_SampleTypes(barcode = t$barcode, typesample = "TP")
Nsample <- TCGAquery_SampleTypes(barcode = n$barcode, typesample = "NT")

# ------------------------------
# 9. Differential Expression Analysis (DEA):
# ------------------------------
library(TCGAbiolinks)

dataDEGs <- TCGAanalyze_DEA(
  mat1 = mrnaEx[, Tsample],
  mat2 = mrnaEx[, Nsample],
  Cond1type = "Tumor",
  Cond2type = "Normal", method = "glmLRT")

write.table(dataDEGs, "D:/my freelancing/Tahmine/tcga/deg.txt", quote = FALSE, sep = "\t")

# ------------------------------
# 10. lncRNA Extraction:
# ------------------------------
lnc <- dataDEGs[dataDEGs$gene_type == "lncRNA", ]
lncEx <- mrnaEx[rownames(lnc), ]
rownames(lncEx) <- lnc$gene_name
siglnc <- lnc[lnc$FDR < 0.01, ]

write.table(lnc, "D:/my freelancing/Tahmine/tcga/lncdeg.txt", quote = FALSE, sep = "\t")
write.table(siglnc, "D:/my freelancing/Tahmine/tcga/siglncdeg.txt", quote = FALSE, sep = "\t")
write.table(lncEx, "D:/my freelancing/Tahmine/tcga/lncEx.txt", quote = FALSE, sep = "\t")

# ------------------------------
# 11. Volcano Plot:
# ------------------------------
library(EnhancedVolcano)
library(ggplot2)

keyvals <- ifelse(
  lnc$logFC < 0 & lnc$FDR < 0.01, 'royalblue4',
  ifelse(lnc$logFC > 0 & lnc$FDR < 0.01, 'green4', no = "cyan3"))

names(keyvals)[keyvals == 'royalblue4'] <- 'low'
names(keyvals)[keyvals == 'green4'] <- 'high'
names(keyvals)[keyvals == 'cyan3'] <- 'Not Signiificance'

EnhancedVolcano(lnc, x = "logFC", y = "FDR", lab = lnc$gene_name,
                pCutoff = 0.01, FCcutoff = 0, pointSize = 3,
                colCustom = keyvals,
                selectLab = c("CDKN2B-AS1", "LINC02418", "HAGLR", "FEZF1-AS1"),
                labSize = 3,
                colAlpha = 0.8, caption = "", border = "full",
                title = "Volcano Plot", subtitle = "", drawConnectors = TRUE,
                widthConnectors = 0.5, endsConnectors = "both",
                arrowheads = FALSE,
                titleLabSize = 15, borderWidth = 0.7, legendPosition = "top") +
  theme_light() + theme(legend.title = element_blank())

# ------------------------------
# 12. Expression Labeling:
# ------------------------------
lnc$ExpressionType[lnc$logFC > 0] <- "up_regulated"
lnc$ExpressionType[lnc$logFC < 0] <- "down_regulated"
siglnc <- lnc[lnc$FDR < 0.01, ]

write.table(lnc, "D:/my freelancing/Tahmine/tcga/lncdeg.txt", quote = FALSE, sep = "\t")
write.table(siglnc, "D:/my freelancing/Tahmine/tcga/siglncdeg.txt", quote = FALSE, sep = "\t")

# ------------------------------
# 13. ID Conversion:
# ------------------------------
library(org.Hs.eg.db)
library(AnnotationDbi)

siglnc$entrezid <- mapIds(org.Hs.eg.db, keys = rownames(siglnc),
                          column = "ENTREZID", keytype = "ENSEMBL")

# ------------------------------
# 14. KEGG Enrichment:
# ------------------------------
library(clusterProfiler)
library(DOSE)
library(enrichplot)

KEGG <- enrichKEGG(gene = siglnc$entrezid, organism = "hsa", keyType = "kegg",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                   pAdjustMethod = "BH", universe = NULL)

dotplot(KEGG, color = "p.adjust")

# ------------------------------
# 15. GO Enrichment:
# ------------------------------
go <- enrichGO(siglnc$entrezid, OrgDb = org.Hs.eg.db, ont = "all",
               pvalueCutoff = 0.05, pAdjustMethod = "BH", 
               qvalueCutoff = 0.05)

dotplot(go, split = "ONTOLOGY", color = "p.adjust") + 
  facet_grid(ONTOLOGY ~ ., scales = "free") +
  theme(axis.text.y = element_text(size = 9, colour = "black"),
        text = element_text(size = 8))
