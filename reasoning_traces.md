# Reasoning Traces: CyTOF Harmony Batch Correction in Seurat

## Decision 1: Identifying Type Markers

**Context**: CyTOF panels distinguish between "type" markers (surface markers used for cell type identification) and "state" markers (intracellular markers reflecting cell function/activation). The prompt specifies using only "type" markers for clustering and UMAP.

**Reasoning**: Examined the FCS file parameter descriptions. The panel contains 30 channels. I excluded:
- Bead channels: Ce140Di (beads1), Ce142Di (beads2), Eu153Di (beads3), Lu175Di (beads4), Lu176Di (beads5)
- Cell identification channels: Ir191Di (cells1), Ir193Di (cells2)
- Viability: Pt195Di (livedead)
- Housekeeping: Event #, Event_length, Time

This left 19 surface/type markers: CD27, CD3, CD123, CD33, CD14, CD61, CD19, CD4, CD8A, CD16, CD235, CD20, CD66, CD45, CD11C, CD45RA, CD38, HLADR, CD56. All canonical surface markers for immunophenotyping. The panel does not contain intracellular/state markers (no phospho-proteins, no transcription factors), so all non-housekeeping channels are type markers.

## Decision 2: Asinh Transformation with Cofactor 5

**Context**: CyTOF data requires transformation before analysis. The standard approach is the inverse hyperbolic sine (asinh) transformation.

**Reasoning**: Used cofactor of 5, the standard for CyTOF data (as opposed to cofactor 150 for fluorescence flow cytometry). This compresses the dynamic range while preserving the ability to distinguish populations, and is the de facto standard in the field (Bendall et al., 2011).

## Decision 3: Subsampling Strategy

**Context**: 36 FCS files with 50,000-133,000 cells each. Need to subsample to 2000 per file.

**Reasoning**: Random subsampling of 2000 cells per file gives 72,000 total cells (36 x 2000). This ensures equal representation across samples. For UMAP, an additional subsample of 10,000 cells from the 72,000 was used, as specified in the prompt.

## Decision 4: Seurat Object Construction Without Standard scRNA-seq Pipeline

**Context**: CyTOF data is not gene expression data. The standard Seurat workflow (LogNormalize -> ScaleData -> PCA -> FindNeighbors -> UMAP) is inappropriate because:
1. Data is already transformed (asinh)
2. There is no gene expression to normalize
3. PCA is unnecessary with only 19 markers (vs. thousands of genes)

**Reasoning**:
- Created a Seurat object with the asinh-transformed data as both the `counts` and `data` layers (skipping normalization)
- Stored the 19 marker expression values directly as a dimensional reduction ("markers")
- Ran FindNeighbors and UMAP directly on this marker space (dims 1:19)
- This approach treats each marker as an independent dimension, which is appropriate for the low-dimensional CyTOF setting

## Decision 5: Batch Effect Assessment -- KNN Neighborhood Entropy

**Context**: Need a clear, quantitative way to detect and visualize batch effects.

**Reasoning**: Chose KNN neighborhood entropy as the primary metric. For each cell:
1. Look at its k-nearest neighbors (from Seurat's FindNeighbors graph)
2. Tally the center labels among those neighbors
3. Compute Shannon entropy of the center distribution

**Why this approach**:
- It's **per-cell**, giving spatial resolution -- you can see exactly where on the UMAP batch effects concentrate
- It produces a **distribution** that can be compared as overlaid histograms
- It collapses to a clean **summary statistic** (mean +/- SD)
- Maximum entropy = log2(6) = 2.585 for 6 centers (perfectly mixed neighborhoods)
- It leverages the same KNN graph that Seurat already computes for clustering

**Pre-Harmony results**:
- Mean per-cell KNN entropy: 1.51 (SD 0.73) out of max 2.585
- This means the average cell's neighborhood is only 58.4% of maximally mixed
- Median was 1.72, indicating a right-skewed distribution with a long low-entropy tail (many cells in batch-dominated neighborhoods)

**Additionally**: Computed per-cluster center entropy from the cluster x center contingency table. Pre-Harmony, clusters 13 (entropy 0.53), 19 (0.83), and 24 (0.07) were almost entirely from Center 4.

## Decision 6: Harmony Integration

**Context**: Need to correct batch effects using Harmony within Seurat.

**Reasoning**:
- Ran Harmony on the marker embedding space (the 19-dimensional type marker space), correcting for "center"
- Used `project.dim = FALSE` because we are not using PCA (Harmony's default projection step requires scale.data which we intentionally skip)
- Harmony converged in 3 iterations

## Decision 7: Post-Harmony Re-evaluation

**Context**: After Harmony, need to re-run everything on the corrected embeddings and re-assess.

**Reasoning**: Re-ran FindNeighbors, FindClusters, and UMAP on the Harmony-corrected 19-dimensional embedding. Then recomputed KNN neighborhood entropy using the new neighbor graph.

**Post-Harmony results**:
- Mean per-cell KNN entropy: 1.83 (SD 0.49) out of max 2.585
- That's 70.9% of maximally mixed (up from 58.4%)
- The SD narrowed from 0.73 to 0.49, meaning fewer outlier cells with very low entropy
- Per-cluster entropy: mean increased from 2.19 to 2.51. Only 1 cluster (cluster 25, 37 cells) has entropy below 1.5

## Decision 8: Visualization Strategy

**Context**: Need to make the batch effect and its correction clearly visible.

**Reasoning**: Produced multiple complementary visualizations:
1. **UMAP colored by center**: Classic view, shows spatial segregation
2. **UMAP colored by KNN entropy**: Heat map showing where batch effects concentrate spatially
3. **Entropy histogram**: Overlaid pre/post distributions make the shift obvious
4. **Mean entropy bar plot with error bars**: Clean summary comparison
5. **Per-cluster entropy bar plot**: Side-by-side pre/post for every cluster, making it clear which clusters improved
6. **Combined panel**: All key views in one figure for presentation

The inferno colormap (dark = low entropy, bright = high entropy) was chosen for the UMAP entropy plots because it makes low-entropy (batch-affected) regions visually prominent.
