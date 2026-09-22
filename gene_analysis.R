## gene_analysis.R
## Simple gene expression comparison: control vs treated
## Uses only base R: named vectors, matrix, vectorized arithmetic,
## conditional subsetting, ifelse(), and factors. No loops needed.

# ---- 1. Raw data as named vectors ---------------------------------------
gene_names <- c("BRCA1", "TP53", "EGFR", "MYC", "PTEN", "KRAS")

control <- c(5.2, 7.8, 3.1, 9.4, 6.0, 4.5)
treated <- c(8.9, 7.6, 6.7, 12.1, 2.3, 9.8)

names(control) <- gene_names
names(treated) <- gene_names

# ---- 2. Fold change (vector arithmetic, element-wise) --------------------
fold_change <- treated / control
names(fold_change) <- gene_names   # keep names after arithmetic

cat("=== Fold change per gene ===\n")
print(round(fold_change, 2))

# ---- 3. Conditional subsetting: genes with fold_change > 1.5 -------------
big_movers <- fold_change[fold_change > 1.5]

cat("\n=== Genes with fold change > 1.5 ===\n")
print(round(big_movers, 2))

# ---- 4. Combine into one matrix (genes x conditions) ----------------------
expr_matrix <- matrix(
  c(control, treated),
  nrow = length(gene_names),
  ncol = 2,
  dimnames = list(gene_names, c("control", "treated"))
)

cat("\n=== Expression matrix (genes x conditions) ===\n")
print(expr_matrix)

# ---- 5. Row-wise mean expression, no loop ---------------------------------
gene_avg <- rowMeans(expr_matrix)

cat("\n=== Average expression per gene (rowMeans) ===\n")
print(round(gene_avg, 2))

# ---- 6. Classify genes with ifelse(), store as a factor -------------------
# vectorized, nested ifelse -- no loop, no per-gene if/else
call_label <- ifelse(
  fold_change > 1.2, "upregulated",
  ifelse(fold_change < 0.8, "downregulated", "stable")
)

gene_call <- factor(
  call_label,
  levels = c("upregulated", "stable", "downregulated")
)
names(gene_call) <- gene_names

cat("\n=== Classification (factor) ===\n")
print(gene_call)

cat("\n=== How many genes per class ===\n")
print(table(gene_call))

# ---- 7. Assemble a tidy summary table for output --------------------------
summary_df <- data.frame(
  gene          = gene_names,
  control       = control,
  treated       = treated,
  fold_change   = round(fold_change, 3),
  classification = gene_call,
  row.names     = NULL
)

cat("\n=== Summary table ===\n")
print(summary_df)

# ---- 8. Save everything to results.txt -------------------------------------
sink("results.txt")

cat("Gene Expression Analysis: control vs treated\n")
cat("=============================================\n\n")

cat("Fold change per gene:\n")
print(round(fold_change, 2))

cat("\nGenes with fold change > 1.5:\n")
print(round(big_movers, 2))

cat("\nExpression matrix (genes x conditions):\n")
print(expr_matrix)

cat("\nAverage expression per gene (rowMeans):\n")
print(round(gene_avg, 2))

cat("\nClassification counts:\n")
print(table(gene_call))

cat("\nFull summary table:\n")
print(summary_df, row.names = FALSE)

sink()

cat("\nSaved results to results.txt\n")
