# 1. Data Collection:

# Query the TCGA dataset for Colon Adenocarcinoma (COAD) gene expression data using the STAR counts method
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")

# Download the queried data
GDCdownload(query)

# 2. Data Preparation:

# Prepare the data for further analysis
prep <- GDCprepare(query)

# 3. Preprocessing:

# Preprocess the data by removing genes with low correlation across samples (correlation threshold of 0.6)
data_prepro <- TCGAanalyze_Preprocessing(prep, cor.cut = 0.6)

# 4. Normalization:

# Normalize the gene expression data using gene length normalization
dataNorm <- TCGAanalyze_Normalization(data_prepro, geneInfo = geneInfoHT, method = "geneLength")


# 5. Filtering:

# Apply quantile filtering to remove genes with low expression across samples (cutoff at 25th percentile)
dataFilt <- TCGAanalyze_Filtering(dataFilt, method = "quantile", qnt.cut = 0.25)


# 6. Voom Transformation:

# Apply the voom transformation to the filtered data (variance modeling) and plot the results
v <- voom(dataFilt, plot = T)

# Extract the transformed expression values from the voom object
E <- v$E

# Exponentiate the transformed values and log-transform to prepare the data for analysis
ex <- 2^E
mrnaEx <- log2(ex + 1)

# 7. Sample Annotation:

# Create a data frame for sample barcodes and assign sample type (Normal or Tumor) based on the barcode
samples <- data.frame(colnames(mrnaEx))
colnames(samples) <- "barcode"
samples$substr <- substr(samples$barcode, 14, 15)
samples$Type[samples$substr > 10] <- "Normal"
samples$Type[samples$substr < 9] <- "Tumor"

# Write the sample information to a file
write.table(samples, "D:/my freelancing/Tahmine/tcga/samples.txt", quote = F, sep = "\t")


# 8. Group Selection:

# Extract tumor and normal sample subsets
n <- bar[bar$substr == 11, ]
t <- bar[bar$substr == "01", ]

# Retrieve sample types for tumor (TP) and normal (NT) samples
Tsample <- TCGAquery_SampleTypes(barcode = t$barcode, typesample = "TP")
Nsample <- TCGAquery_SampleTypes(barcode = n$barcode, typesample = "NT")

# 9. Differential Expression Analysis (DEA):

#Perform differential expression analysis (DEA) between tumor and normal samples using GLM-based testing
dataDEGs <- TCGAanalyze_DEA(
  mat1 = mrnaEx[, Tsample],
  mat2 = mrnaEx[, Nsample],
  Cond1type = "Tumor",
  Cond2type = "Normal", method = "glmLRT")

# Write the results of differential expression analysis to a file
write.table(dataDEGs, "D:/my freelancing/Tahmine/tcga/deg.txt", quote = F, sep = "\t")


# 10. lncRNA Extraction:

# Filter for lncRNA genes and extract their expression data
lnc <- dataDEGs[dataDEGs$gene_type == "lncRNA", ]
lncEx <- mrnaEx[rownames(lnc), ]

# Rename rows in the expression matrix with gene names
rownames(lncEx) <- lnc$gene_name

# Filter significant lncRNAs based on FDR < 0.01
siglnc <- lnc[lnc$FDR < 0.01, ]

# Write filtered lncRNAs and significant lncRNAs to files
write.table(lnc, "D:/my freelancing/Tahmine/tcga/lncdeg.txt", quote = F, sep = "\t")
write.table(siglnc, "D:/my freelancing/Tahmine/tcga/siglncdeg.txt", quote = F, sep = "\t")
write.table(lncEx, "D:/my freelancing/Tahmine/tcga/lncEx.txt", quote = F, sep = "\t")


# 11. Volcano Plot:

# Generate a volcano plot to visualize differentially expressed lncRNAs with color coding based on logFC and FDR
keyvals <- ifelse(
  lnc$logFC < 0 & lnc$FDR < 0.01, 'royalblue4',
  ifelse(lnc$logFC > 0 & lnc$FDR < 0.01, 'green4', no = "cyan3"))

names(keyvals)[keyvals == 'royalblue4'] <- 'low'
names(keyvals)[keyvals == 'green4'] <- 'high'
names(keyvals)[keyvals == 'cyan3'] <- 'Not Signiificance'

# Create the volcano plot using the EnhancedVolcano function
EnhancedVolcano(lnc, x = "logFC", y = "FDR", lab = lnc$gene_name,
                pCutoff = 0.01, FCcutoff = 0, pointSize = 3,
                colCustom = keyvals,
                selectLab = c("CDKN2B-AS1", "LINC02418", "HAGLR", "FEZF1-AS1"),
                labSize = 3,
                colAlpha = 0.8, caption = "", border = "full",
                title = "Volcano Plot", subtitle = "", drawConnectors = TRUE,
                widthConnectors = 0.5, endsConnectors = "both",
                arrowheads = F,
                titleLabSize = 15, borderWidth = 0.7, legendPosition = "top") +
  theme_light() + theme(legend.title = element_blank())

# 12. Expression Labeling:

# Classify lncRNAs as up-regulated or down-regulated based on logFC
lnc$ExpressionType[lnc$logFC > 0] <- "up_regulated"
lnc$ExpressionType[lnc$logFC < 0] <- "down_regulated"

# Filter significant lncRNAs and write results to files
siglnc <- lnc[lnc$FDR < 0.01, ]
write.table(lnc, "D:/my freelancing/Tahmine/tcga/lncdeg.txt", quote = F, sep = "\t")
write.table(siglnc, "D:/my freelancing/Tahmine/tcga/siglncdeg.txt", quote = F, sep = "\t")

##############
# 13. ID Conversion:

# Enrichment analysis for KEGG pathways and GO terms
library(DOSE)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# Map lncRNA genes to their corresponding Entrez IDs
siglnc$entrezid <- mapIds(org.Hs.eg.db, keys = rownames(siglnc), column = "ENTREZID", 
                          keytype = "ENSEMBL")


# 14. KEGG Enrichment:

# Perform KEGG pathway enrichment analysis for significant lncRNAs
KEGG <- enrichKEGG(gene = siglnc$entrezid, organism = "hsa", keyType = "kegg",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                   pAdjustMethod = "BH", universe = NULL)

# Plot the results of the KEGG enrichment analysis
dotplot(KEGG, color = "p.adjust")

# 15. GO Enrichment:

# Perform Gene Ontology (GO) enrichment analysis for significant lncRNAs
go <- enrichGO(siglnc$entrezid, OrgDb = org.Hs.eg.db, ont = "all",
               pvalueCutoff = 0.05, pAdjustMethod = "BH", 
               qvalueCutoff = 0.05)

# Plot the results of the GO enrichment analysis
dotplot(go, split = "ONTOLOGY", color = "p.adjust") + 
  facet_grid(ONTOLOGY ~ ., scales = "free") +
  theme(axis.text.y = element_text(size = 9, colour = "black"),
        text = element_text(size = 8))
