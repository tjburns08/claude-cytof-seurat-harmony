##############################################################################
# CyTOF Harmony Batch Correction Analysis within Seurat
##############################################################################

library(flowCore)
library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

set.seed(42)

# Directories
data_dir <- "data"
output_dir <- "output"

##############################################################################
# Helper: Shannon entropy (normalized to [0,1] given n_categories)
##############################################################################
shannon_entropy <- function(x) {
  p <- x / sum(x)
  p <- p[p > 0]
  -sum(p * log2(p))
}

##############################################################################
# Helper: Per-cell KNN entropy
# For each cell, look at its k-nearest neighbors, tally up the center labels,
# compute Shannon entropy. Returns a vector of per-cell entropy values.
##############################################################################
compute_knn_entropy <- function(seurat_obj, nn_name = "CyTOF_nn", meta_col = "center") {
  nn_idx <- seurat_obj@graphs[[nn_name]]
  centers <- seurat_obj@meta.data[[meta_col]]
  n_cells <- length(centers)
  cell_entropy <- numeric(n_cells)

  for (i in seq_len(n_cells)) {
    # Get neighbor indices from the NN graph (nonzero entries in row i)
    neighbors_i <- nn_idx[i, ]
    neighbor_idx <- which(neighbors_i > 0)
    neighbor_centers <- centers[neighbor_idx]
    center_counts <- table(factor(neighbor_centers, levels = unique(centers)))
    cell_entropy[i] <- shannon_entropy(center_counts)
  }
  return(cell_entropy)
}

##############################################################################
# 1. Read FCS files, asinh transform, add metadata, subsample
##############################################################################

fcs_files <- list.files(data_dir, pattern = "\\.fcs$", full.names = TRUE)
cat("Found", length(fcs_files), "FCS files\n")

# Read one file to identify markers
ff_example <- read.FCS(fcs_files[1], truncate_max_range = FALSE)
params <- pData(parameters(ff_example))
cat("All channels:\n")
print(params[, c("name", "desc")])

# Define type (surface) markers
exclude_desc <- c("beads1", "beads2", "beads3", "beads4", "beads5",
                   "cells1", "cells2", "livedead", "Event #", "Event_length", "Time")
type_marker_idx <- which(!params$desc %in% exclude_desc)
type_marker_names <- params$desc[type_marker_idx]
cat("\nType markers (", length(type_marker_names), "):\n")
print(type_marker_names)

# Read all files
all_data <- list()
all_meta <- list()

for (f in fcs_files) {
  ff <- read.FCS(f, truncate_max_range = FALSE)
  mat <- exprs(ff)

  fname <- basename(f)
  parts <- strsplit(fname, "_")[[1]]
  center <- paste(parts[1], parts[2], sep = "_")
  sample_id <- paste(parts[1], parts[2], parts[3], parts[4], sep = "_")

  # asinh transform (cofactor 5 for CyTOF)
  mat_transformed <- asinh(mat / 5)

  # Subsample to 2000 cells
  n_cells <- nrow(mat_transformed)
  if (n_cells > 2000) {
    idx <- sample(n_cells, 2000)
    mat_transformed <- mat_transformed[idx, , drop = FALSE]
  }

  colnames(mat_transformed) <- params$desc
  mat_type <- mat_transformed[, type_marker_names, drop = FALSE]

  meta <- data.frame(
    center = rep(center, nrow(mat_type)),
    sample_id = rep(sample_id, nrow(mat_type)),
    stringsAsFactors = FALSE
  )

  all_data[[f]] <- mat_type
  all_meta[[f]] <- meta
  cat("Read", fname, ":", n_cells, "cells -> subsampled to", nrow(mat_type), "\n")
}

# Combine
combined_data <- do.call(rbind, all_data)
combined_meta <- do.call(rbind, all_meta)
rownames(combined_data) <- paste0("cell_", seq_len(nrow(combined_data)))
rownames(combined_meta) <- rownames(combined_data)

cat("\nCombined data:", nrow(combined_data), "cells x", ncol(combined_data), "markers\n")
cat("Centers:", unique(combined_meta$center), "\n")
cat("Samples:", length(unique(combined_meta$sample_id)), "\n")

##############################################################################
# 2. Create Seurat object (no normalization, no scaling, no PCA)
##############################################################################

seurat_obj <- CreateSeuratObject(
  counts = t(combined_data),
  meta.data = combined_meta,
  assay = "CyTOF"
)

LayerData(seurat_obj, assay = "CyTOF", layer = "data") <- LayerData(seurat_obj, assay = "CyTOF", layer = "counts")

# Store marker expression as a dimensional reduction
marker_matrix <- t(as.matrix(GetAssayData(seurat_obj, layer = "data")))
colnames(marker_matrix) <- paste0("marker_", seq_len(ncol(marker_matrix)))
seurat_obj[["markers"]] <- CreateDimReducObject(
  embeddings = marker_matrix, key = "marker_", assay = "CyTOF"
)

cat("\nSeurat object created:", ncol(seurat_obj), "cells,", nrow(seurat_obj), "features\n")

##############################################################################
# 3. Pre-Harmony: UMAP and clustering on type markers directly
##############################################################################

cat("\n=== PRE-HARMONY RUN ===\n")

# Subsample 10000 cells for UMAP
n_total <- ncol(seurat_obj)
if (n_total > 10000) {
  umap_cells <- sample(colnames(seurat_obj), 10000)
} else {
  umap_cells <- colnames(seurat_obj)
}

# Find neighbors and cluster using marker space
seurat_obj <- FindNeighbors(seurat_obj, reduction = "markers",
                            dims = 1:ncol(marker_matrix))
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)

# Store pre-harmony cluster labels before they get overwritten
seurat_obj$clusters_pre <- Idents(seurat_obj)

# Run UMAP on 10,000-cell subset
seurat_obj_sub <- subset(seurat_obj, cells = umap_cells)
seurat_obj_sub <- RunUMAP(seurat_obj_sub, reduction = "markers",
                          dims = 1:ncol(marker_matrix),
                          reduction.name = "umap_pre",
                          reduction.key = "UMAPpre_")

cat("Pre-Harmony UMAP computed on", length(umap_cells), "cells\n")
cat("Pre-Harmony clusters found:", length(levels(Idents(seurat_obj))), "\n")

##############################################################################
# 3b. Pre-Harmony: KNN-based per-cell entropy
##############################################################################

cat("\n--- Batch Effect Assessment (Pre-Harmony) ---\n")

max_entropy <- log2(length(unique(seurat_obj$center)))

# Per-cell KNN entropy
cat("Computing per-cell KNN entropy (pre-Harmony)...\n")
knn_entropy_pre <- compute_knn_entropy(seurat_obj, nn_name = "CyTOF_nn", meta_col = "center")
seurat_obj$knn_entropy_pre <- knn_entropy_pre

cat("Per-cell KNN entropy (pre-Harmony):\n")
cat("  Mean:", round(mean(knn_entropy_pre), 4), "\n")
cat("  Median:", round(median(knn_entropy_pre), 4), "\n")
cat("  SD:", round(sd(knn_entropy_pre), 4), "\n")
cat("  Max possible:", round(max_entropy, 4), "\n")

# Per-cluster entropy (from cluster x center contingency table)
cluster_center_table_pre <- table(seurat_obj$clusters_pre, seurat_obj$center)
cluster_entropies_pre <- apply(cluster_center_table_pre, 1, shannon_entropy)
cat("\nPer-cluster center entropy (pre-Harmony):\n")
print(round(cluster_entropies_pre, 3))
cat("Mean cluster entropy:", round(mean(cluster_entropies_pre), 3), "/", round(max_entropy, 3), "\n")

# Transfer KNN entropy to the UMAP subset for plotting
seurat_obj_sub$knn_entropy_pre <- seurat_obj$knn_entropy_pre[umap_cells]

##############################################################################
# 3c. Pre-Harmony plots
##############################################################################

# Seurat DimPlots
p_pre_center <- DimPlot(seurat_obj_sub, reduction = "umap_pre", group.by = "center") +
  ggtitle("Pre-Harmony: Center")
p_pre_sample <- DimPlot(seurat_obj_sub, reduction = "umap_pre", group.by = "sample_id") +
  ggtitle("Pre-Harmony: Sample ID") +
  theme(legend.text = element_text(size = 6))
p_pre_cluster <- DimPlot(seurat_obj_sub, reduction = "umap_pre", group.by = "seurat_clusters") +
  ggtitle("Pre-Harmony: Clusters")

# UMAP colored by KNN entropy
umap_pre_coords <- Embeddings(seurat_obj_sub, "umap_pre")
entropy_df_pre <- data.frame(
  UMAP_1 = umap_pre_coords[, 1],
  UMAP_2 = umap_pre_coords[, 2],
  knn_entropy = seurat_obj_sub$knn_entropy_pre
)
p_pre_entropy_umap <- ggplot(entropy_df_pre, aes(x = UMAP_1, y = UMAP_2, color = knn_entropy)) +
  geom_point(size = 0.3) +
  scale_color_viridis_c(limits = c(0, max_entropy), option = "inferno") +
  theme_minimal() +
  ggtitle("Pre-Harmony: KNN Neighborhood Entropy") +
  labs(color = "Entropy")

# Save pre-harmony plots
png(file.path(output_dir, "pre_harmony_center.png"), width = 800, height = 600)
print(p_pre_center)
dev.off()

png(file.path(output_dir, "pre_harmony_sample.png"), width = 1000, height = 600)
print(p_pre_sample)
dev.off()

png(file.path(output_dir, "pre_harmony_clusters.png"), width = 800, height = 600)
print(p_pre_cluster)
dev.off()

png(file.path(output_dir, "pre_harmony_knn_entropy_umap.png"), width = 800, height = 600)
print(p_pre_entropy_umap)
dev.off()

##############################################################################
# 4. Run Harmony
##############################################################################

cat("\n=== RUNNING HARMONY ===\n")

seurat_obj <- RunHarmony(seurat_obj,
                         group.by.vars = "center",
                         reduction = "markers",
                         reduction.save = "harmony",
                         assay.use = "CyTOF",
                         project.dim = FALSE)

cat("Harmony completed\n")

##############################################################################
# 5. Post-Harmony: Re-run clustering and UMAP
##############################################################################

cat("\n=== POST-HARMONY RE-RUN ===\n")

n_harmony_dims <- ncol(Embeddings(seurat_obj, "harmony"))
seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony",
                            dims = 1:n_harmony_dims)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)

# Store post-harmony cluster labels
seurat_obj$clusters_post <- Idents(seurat_obj)

# Run UMAP on same subset
seurat_obj_sub_post <- subset(seurat_obj, cells = umap_cells)
seurat_obj_sub_post <- RunUMAP(seurat_obj_sub_post, reduction = "harmony",
                               dims = 1:n_harmony_dims,
                               reduction.name = "umap_post",
                               reduction.key = "UMAPpost_")

cat("Post-Harmony UMAP computed on", length(umap_cells), "cells\n")
cat("Post-Harmony clusters found:", length(levels(Idents(seurat_obj))), "\n")

##############################################################################
# 5b. Post-Harmony: KNN-based per-cell entropy
##############################################################################

cat("\n--- Batch Effect Assessment (Post-Harmony) ---\n")

cat("Computing per-cell KNN entropy (post-Harmony)...\n")
knn_entropy_post <- compute_knn_entropy(seurat_obj, nn_name = "CyTOF_nn", meta_col = "center")
seurat_obj$knn_entropy_post <- knn_entropy_post

cat("Per-cell KNN entropy (post-Harmony):\n")
cat("  Mean:", round(mean(knn_entropy_post), 4), "\n")
cat("  Median:", round(median(knn_entropy_post), 4), "\n")
cat("  SD:", round(sd(knn_entropy_post), 4), "\n")
cat("  Max possible:", round(max_entropy, 4), "\n")

# Per-cluster entropy
cluster_center_table_post <- table(seurat_obj$clusters_post, seurat_obj$center)
cluster_entropies_post <- apply(cluster_center_table_post, 1, shannon_entropy)
cat("\nPer-cluster center entropy (post-Harmony):\n")
print(round(cluster_entropies_post, 3))
cat("Mean cluster entropy:", round(mean(cluster_entropies_post), 3), "/", round(max_entropy, 3), "\n")

# Transfer to UMAP subset
seurat_obj_sub_post$knn_entropy_post <- seurat_obj$knn_entropy_post[umap_cells]

##############################################################################
# 5c. Post-Harmony plots
##############################################################################

p_post_center <- DimPlot(seurat_obj_sub_post, reduction = "umap_post", group.by = "center") +
  ggtitle("Post-Harmony: Center")
p_post_sample <- DimPlot(seurat_obj_sub_post, reduction = "umap_post", group.by = "sample_id") +
  ggtitle("Post-Harmony: Sample ID") +
  theme(legend.text = element_text(size = 6))
p_post_cluster <- DimPlot(seurat_obj_sub_post, reduction = "umap_post", group.by = "seurat_clusters") +
  ggtitle("Post-Harmony: Clusters")

# UMAP colored by KNN entropy
umap_post_coords <- Embeddings(seurat_obj_sub_post, "umap_post")
entropy_df_post <- data.frame(
  UMAP_1 = umap_post_coords[, 1],
  UMAP_2 = umap_post_coords[, 2],
  knn_entropy = seurat_obj_sub_post$knn_entropy_post
)
p_post_entropy_umap <- ggplot(entropy_df_post, aes(x = UMAP_1, y = UMAP_2, color = knn_entropy)) +
  geom_point(size = 0.3) +
  scale_color_viridis_c(limits = c(0, max_entropy), option = "inferno") +
  theme_minimal() +
  ggtitle("Post-Harmony: KNN Neighborhood Entropy") +
  labs(color = "Entropy")

# Save post-harmony plots
png(file.path(output_dir, "post_harmony_center.png"), width = 800, height = 600)
print(p_post_center)
dev.off()

png(file.path(output_dir, "post_harmony_sample.png"), width = 1000, height = 600)
print(p_post_sample)
dev.off()

png(file.path(output_dir, "post_harmony_clusters.png"), width = 800, height = 600)
print(p_post_cluster)
dev.off()

png(file.path(output_dir, "post_harmony_knn_entropy_umap.png"), width = 800, height = 600)
print(p_post_entropy_umap)
dev.off()

##############################################################################
# 6. Comparison plots
##############################################################################

# --- KNN Entropy UMAP side by side ---
png(file.path(output_dir, "knn_entropy_umap_comparison.png"), width = 1400, height = 600, res = 120)
print(p_pre_entropy_umap + ggtitle("Pre-Harmony") | p_post_entropy_umap + ggtitle("Post-Harmony"))
dev.off()

# --- KNN Entropy histogram (overlaid) ---
entropy_hist_df <- data.frame(
  entropy = c(knn_entropy_pre, knn_entropy_post),
  stage = rep(c("Pre-Harmony", "Post-Harmony"), each = length(knn_entropy_pre))
)
entropy_hist_df$stage <- factor(entropy_hist_df$stage, levels = c("Pre-Harmony", "Post-Harmony"))

p_entropy_hist <- ggplot(entropy_hist_df, aes(x = entropy, fill = stage)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 50) +
  scale_fill_manual(values = c("Pre-Harmony" = "#E74C3C", "Post-Harmony" = "#2ECC71")) +
  geom_vline(xintercept = max_entropy, linetype = "dashed", color = "black") +
  annotate("text", x = max_entropy - 0.15, y = Inf, label = "Max entropy", vjust = 2, size = 3) +
  theme_minimal() +
  ggtitle("Per-Cell KNN Neighborhood Entropy Distribution") +
  xlab("Shannon Entropy") + ylab("Count") +
  labs(fill = "Stage")

png(file.path(output_dir, "knn_entropy_histogram.png"), width = 900, height = 500, res = 120)
print(p_entropy_hist)
dev.off()

# --- KNN Entropy bar plot (mean +/- SD) ---
entropy_summary <- data.frame(
  stage = factor(c("Pre-Harmony", "Post-Harmony"), levels = c("Pre-Harmony", "Post-Harmony")),
  mean = c(mean(knn_entropy_pre), mean(knn_entropy_post)),
  sd = c(sd(knn_entropy_pre), sd(knn_entropy_post))
)

p_entropy_bar <- ggplot(entropy_summary, aes(x = stage, y = mean, fill = stage)) +
  geom_col(width = 0.6, alpha = 0.8) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  geom_hline(yintercept = max_entropy, linetype = "dashed", color = "black") +
  annotate("text", x = 2.4, y = max_entropy + 0.03, label = paste0("Max = ", round(max_entropy, 2)), size = 3) +
  scale_fill_manual(values = c("Pre-Harmony" = "#E74C3C", "Post-Harmony" = "#2ECC71")) +
  theme_minimal() +
  ggtitle("Mean Per-Cell KNN Entropy") +
  xlab("") + ylab("Shannon Entropy (mean +/- SD)") +
  theme(legend.position = "none") +
  ylim(0, max_entropy + 0.3)

png(file.path(output_dir, "knn_entropy_barplot.png"), width = 500, height = 500, res = 120)
print(p_entropy_bar)
dev.off()

# --- Per-cluster entropy bar plot (pre vs post side by side) ---
# Align cluster labels: use cluster number as factor
cluster_ent_df <- data.frame(
  cluster = c(names(cluster_entropies_pre), names(cluster_entropies_post)),
  entropy = c(cluster_entropies_pre, cluster_entropies_post),
  stage = c(rep("Pre-Harmony", length(cluster_entropies_pre)),
            rep("Post-Harmony", length(cluster_entropies_post)))
)
cluster_ent_df$stage <- factor(cluster_ent_df$stage, levels = c("Pre-Harmony", "Post-Harmony"))
# Order clusters numerically
cluster_ent_df$cluster <- factor(cluster_ent_df$cluster,
                                  levels = sort(unique(as.numeric(as.character(cluster_ent_df$cluster)))))

p_cluster_entropy <- ggplot(cluster_ent_df, aes(x = cluster, y = entropy, fill = stage)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.8) +
  geom_hline(yintercept = max_entropy, linetype = "dashed", color = "black") +
  annotate("text", x = Inf, y = max_entropy + 0.05, label = paste0("Max = ", round(max_entropy, 2)),
           hjust = 1.1, size = 3) +
  scale_fill_manual(values = c("Pre-Harmony" = "#E74C3C", "Post-Harmony" = "#2ECC71")) +
  theme_minimal() +
  ggtitle("Per-Cluster Center Entropy: Pre vs Post Harmony") +
  xlab("Cluster") + ylab("Shannon Entropy") +
  labs(fill = "Stage") +
  theme(axis.text.x = element_text(size = 7))

png(file.path(output_dir, "cluster_entropy_barplot.png"), width = 1200, height = 500, res = 120)
print(p_cluster_entropy)
dev.off()

# --- Big combined comparison ---
pdf(file.path(output_dir, "harmony_comparison.pdf"), width = 16, height = 20)
print(
  (p_pre_center + ggtitle("PRE: Center") | p_post_center + ggtitle("POST: Center")) /
  (p_pre_entropy_umap + ggtitle("PRE: KNN Entropy") | p_post_entropy_umap + ggtitle("POST: KNN Entropy")) /
  (p_pre_cluster + ggtitle("PRE: Cluster") | p_post_cluster + ggtitle("POST: Cluster")) /
  (p_entropy_bar | p_entropy_hist)
)
dev.off()

png(file.path(output_dir, "harmony_comparison.png"), width = 1600, height = 1800, res = 120)
print(
  (p_pre_center + ggtitle("PRE: Center") | p_post_center + ggtitle("POST: Center")) /
  (p_pre_entropy_umap + ggtitle("PRE: KNN Entropy") | p_post_entropy_umap + ggtitle("POST: KNN Entropy")) /
  (p_pre_cluster + ggtitle("PRE: Cluster") | p_post_cluster + ggtitle("POST: Cluster")) /
  (p_entropy_bar | p_entropy_hist)
)
dev.off()

##############################################################################
# 7. Summary
##############################################################################

cat("\n\n========================================\n")
cat("SUMMARY: BATCH EFFECT COMPARISON\n")
cat("========================================\n")
cat("Max possible entropy (6 centers):", round(max_entropy, 4), "\n\n")
cat("--- Per-Cell KNN Entropy ---\n")
cat("Pre-Harmony  mean:", round(mean(knn_entropy_pre), 4), " SD:", round(sd(knn_entropy_pre), 4), "\n")
cat("Post-Harmony mean:", round(mean(knn_entropy_post), 4), " SD:", round(sd(knn_entropy_post), 4), "\n\n")
cat("--- Per-Cluster Entropy ---\n")
cat("Pre-Harmony  mean:", round(mean(cluster_entropies_pre), 4), "\n")
cat("Post-Harmony mean:", round(mean(cluster_entropies_post), 4), "\n")

# Save the Seurat object
saveRDS(seurat_obj, file.path(output_dir, "seurat_harmony.rds"))
cat("\nSeurat object saved to output/seurat_harmony.rds\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
