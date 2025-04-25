# Diagnostic-lncRNA-Biomarkers-in-Colorectal-Cancer-A-TCGA-RNA-Seq-Analysis
This pipeline focuses on identifying differentially expressed lncRNAs between tumor and normal samples from the TCGA-COAD cohort. 

1. Data Collection
RNA-Seq data for TCGA-COAD is retrieved from the GDC portal in raw count format (STAR aligner output).

2. Data Preparation
Data is processed into an R-compatible format using GDCprepare for further analysis.

3. Preprocessing
Samples with low correlation (cor < 0.6) are filtered out to reduce noise and improve data quality.

4. Normalization
Gene expression is normalized based on gene length to allow accurate comparison across genes.

5. Filtering
Genes with expression levels below the first quartile (25%) are removed to exclude low-expressed, non-informative genes.

6. Voom Transformation
Voom transformation is applied to stabilize variance and prepare data for linear modeling.

7. Sample Annotation
Sample barcodes are parsed to classify samples as tumor or normal based on barcode structure.

8. Group Selection
Tumor and normal samples are selected for differential expression analysis.

9. Differential Expression Analysis (DEA)
A generalized linear model (glmLRT) is applied to identify differentially expressed genes between tumor and normal samples.

10. lncRNA Extraction
Significant lncRNAs (FDR < 0.01) are extracted based on gene annotations from the list of DEGs.

11. Volcano Plot
A volcano plot is generated using EnhancedVolcano to visualize upregulated and downregulated lncRNAs.

12. Expression Labeling
lncRNAs are labeled as "upregulated" or "downregulated" based on log fold change (logFC) values.

13. ID Conversion
Ensembl IDs of significant lncRNAs are converted to Entrez IDs using the org.Hs.eg.db database.

14. KEGG Enrichment
KEGG pathway enrichment analysis is performed to identify relevant signaling and metabolic pathways.

15. GO Enrichment
Gene Ontology enrichment is conducted across Biological Process (BP), Cellular Component (CC), and Molecular Function (MF) categories.
