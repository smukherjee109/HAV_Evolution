# ####################################################################
# DESeq2 Pipeline: GSE234478 (HAV-Infected vs Mock Huh7.5 Cells)
# ####################################################################

# --- 1. Set Working Directory ---
# Update this to where your featureCounts output is saved
setwd("/mnt/c/College/Kusuma School of Biological Sciences (PhD)/Project")

# --- 2. Install and Load Packages ---
packages <- c("DESeq2", "tidyverse", "pheatmap", "ggrepel", 
              "RColorBrewer", "AnnotationDbi", "org.Hs.eg.db")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    message(paste("Installing", pkg, "..."))
    tryCatch({
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }, error = function(e) {
      install.packages(pkg)
    })
    library(pkg, character.only = TRUE)
  }
}

# --- 3. Load Raw Count Data ---
# Replace with the name of your featureCounts output file
counts_file <- "GSE234478_featureCounts.txt" 
countData <- read.delim(counts_file, header = TRUE, row.names = 1, sep = "\t")

# Ensure count data are integers
countData <- round(countData)

# --- 4. Prepare Metadata (colData) ---
# IMPORTANT: Replace these with the EXACT column names from your countData file
sample_names <- c(
  "Mock_Rep1", "Mock_Rep2", "Mock_Rep3", 
  "HAV_Rep1",  "HAV_Rep2",  "HAV_Rep3"
)

conditions <- c(
  "Mock", "Mock", "Mock",
  "HAV",  "HAV",  "HAV"
)

colData <- data.frame(
  row.names = sample_names,
  condition = factor(conditions)
)

# Set "Mock" as the reference baseline. 
# Positive fold changes will mean UPREGULATED in HAV infection.
colData$condition <- relevel(colData$condition, ref = "Mock")

# Ensure count data columns match metadata rows perfectly
countData <- countData[, rownames(colData)]

# --- 5. Create DESeqDataSet Object & Filter ---
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)

# Filter low counts (genes must have at least 10 reads total)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# --- 6. Run DESeq2 ---
dds <- DESeq(dds)

# Create output directory
output_dir <- "DESeq2_Results_GSE234478"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# --- 7. Extract Results, Annotate, and Save ---
res_hav_vs_mock <- results(dds, contrast = c("condition", "HAV", "Mock"))
res_ordered <- res_hav_vs_mock[order(res_hav_vs_mock$padj), ]

res_df <- as.data.frame(res_ordered) %>%
  rownames_to_column(var="entrez") %>% # Assuming featureCounts used Entrez IDs
  mutate(significance = ifelse(padj < 0.05, "Significant", "Not Significant"))

# Annotate Entrez IDs to Symbols and Descriptions
res_df$symbol <- mapIds(org.Hs.eg.db, keys = res_df$entrez, 
                        column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
res_df$description <- mapIds(org.Hs.eg.db, keys = res_df$entrez, 
                             column = "GENENAME", keytype = "ENTREZID", multiVals = "first")

write.csv(res_df, file = file.path(output_dir, "HAV_vs_Mock_annotated_results.csv"), row.names = FALSE)

# --- 8. Generate Volcano Plot ---
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.5) +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw() + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "Volcano Plot: HAV Infection vs Mock (GSE234478)",
       x = expression(log[2]("Fold Change")),
       y = expression(-log[10]("Adjusted p-value")))

ggsave(file.path(output_dir, "Volcano_HAV_vs_Mock.png"), plot = volcano_plot, width = 8, height = 6)


# --- 9. Export Data for GSEA Desktop Software ---
gsea_output_dir <- file.path(output_dir, "GSEA_Files")
if (!dir.exists(gsea_output_dir)) {
  dir.create(gsea_output_dir)
}

# 9a. Create Expression Dataset File (.txt)
vst_counts <- vst(dds, blind = FALSE)
vst_matrix <- assay(vst_counts)

gsea_expression_df <- as.data.frame(vst_matrix) %>%
  rownames_to_column(var = "original_id")

gsea_expression_df$row_var <- rowVars(vst_matrix)

# Map Entrez to Symbol for GSEA
gsea_expression_df$NAME <- mapIds(org.Hs.eg.db, keys = gsea_expression_df$original_id, 
                                  column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")

gsea_expression_df <- gsea_expression_df %>%
  filter(!is.na(NAME)) %>%          # Remove genes that didn't map to a symbol
  mutate(DESCRIPTION = NAME) %>%
  dplyr::select(NAME, DESCRIPTION, everything(), -row_var, -original_id)

write.table(gsea_expression_df,
            file = file.path(gsea_output_dir, "GSE234478_expression_data.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 9b. Create Phenotype Label File (.cls)
phenotypes <- dds$condition
class_levels <- levels(phenotypes)
num_samples <- length(phenotypes)
num_classes <- length(class_levels)

cls_file_path <- file.path(gsea_output_dir, "GSE234478_phenotypes.cls")
file_conn <- file(cls_file_path, "w")

writeLines(paste(num_samples, num_classes, 1, sep = " "), file_conn)
writeLines(paste("#", paste(class_levels, collapse = " ")), file_conn)
writeLines(paste(as.character(phenotypes), collapse = " "), file_conn)

close(file_conn)

message("GSE234478 Pipeline Complete. Results, plots, and GSEA files are saved in the output directory.")