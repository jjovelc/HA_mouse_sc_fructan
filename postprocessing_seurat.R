# Load required libraries
# Load required libraries
library(Seurat)    # Main library for scRNAseq analysis
library(tximport)  # Used to import text files produced by salmon alevin
library(tidyverse) # Includes several R packages needed
library(patchwork)

# Set the working directory
setwd('/Users/juanjovel/OneDrive/jj/UofC/data_analysis/heatherArmstrong/pediatric_MS/sc_mice/matrices')
dir <- getwd()

metadata <- read.table('metadata.txt', sep = '\t', header = T, row.names = 1)

# Capture all .h5 files in the "cellRangerResults" directory
matrices_dirs <- basename(list.dirs(dir, recursive = FALSE))

# Initialize lists to store Seurat gut_mergedects for spine and gut samples
seurat_list_spine <- list()
seurat_list_gut   <- list()

# Iterate over each matrix dir and import each dataset
for (i in seq_along(matrices_dirs)) {
  # Extract the file name without the path for easier checking
  dir_name <- matrices_dirs[i]
  sample <- gsub("_.*$", "", dir_name)
  
  # Read the .h5 file into R
  data <- Read10X(data.dir = dir_name)
  
  # Check the prefix and assign to the appropriate list
  if (startsWith(dir_name, "B")) {
    # Spine sample (starting with "B")
    if (metadata[sample, "Diet"] == 'B-fructan'){
      obj_name <- paste0(gsub("_.*$", "", dir_name), "_fructan")
      seurat_object <- CreateSeuratObject(counts = data, project = obj_name)
      seurat_list_spine[[obj_name]] <- seurat_object
    } else if (metadata[sample, "Diet"] == 'Control'){
      obj_name <- paste0(gsub("_.*$", "", dir_name), "_control")
      seurat_object <- CreateSeuratObject(counts = data, project = obj_name)
      seurat_list_spine[[obj_name]] <- seurat_object
    }
  } else if (startsWith(dir_name, "C")) {
    # Gut sample (starting with "C")
    if (metadata[sample, "Diet"] == 'B-fructan'){
      obj_name <- paste0(gsub("_.*$", "", dir_name), "_fructan")
      seurat_object <- CreateSeuratObject(counts = data, project = obj_name)
      seurat_list_gut[[obj_name]] <- seurat_object
    } else if (metadata[sample, "Diet"] == 'Control'){
      obj_name <- paste0(gsub("_.*$", "", dir_name), "_control")
      seurat_object <- CreateSeuratObject(counts = data, project = obj_name)
      seurat_list_gut[[obj_name]] <- seurat_object
    }
  }
}

# Function to generate QC violin plots
generateQCPlot <- function(obj, prefix, suffix = "") {
  filename <- paste0(prefix, '_seurat_QCplot', suffix, '.png')
  png(filename, width = 12, height = 9, units = 'in', pointsize = 24, res = 300)
  p <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, cols = 'dodgerblue')
  print(p)
  dev.off()
}

# Function to generate feature relationship scatter plots
generateFeatureRelationshipPlot <- function(obj, prefix, suffix = "") {
  filename <- paste0(prefix, '_seurat_feature-relationship', suffix, '.png')
  png(filename, width = 12, height = 9, units = 'in', pointsize = 24, res = 300)
  plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(plot1 + plot2)
  dev.off()
}

# Preprocess individual libraries
preprocess_object <- function(seurat_list){
  for (i in seq_along(seurat_list)) {
    # Access the Seurat gut_mergedect using the index
    obj <- seurat_list[[i]]
    prefix <- gsub("_.*$", "", names(seurat_list)[i])
    
    ### QC ###
    # Ensure "nFeature_RNA" and "nCount_RNA" are calculated from the counts layer explicitly
    if (!"nFeature_RNA" %in% colnames(obj@meta.data)) {
      obj <- AddMetaData(obj, metadata = Matrix::colSums(obj@assays$RNA$counts > 0), col.name = "nFeature_RNA")
    }
    if (!"nCount_RNA" %in% colnames(obj@meta.data)) {
      obj <- AddMetaData(obj, metadata = Matrix::colSums(obj@assays$RNA$counts), col.name = "nCount_RNA")
    }
    
    # Calculate the percentage of mitochondrial reads (explicitly specifying the assay)
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-", assay = "RNA")
    # Generate initial QC plot and feature relationship plot
    
    generateQCPlot(obj, prefix)
    generateFeatureRelationshipPlot(obj, prefix)
    
    # Filter cells based on QC criteria, ensuring no missing values
    if (startsWith(prefix, "B")) {
    
      obj <- subset(obj, subset = !is.na(nFeature_RNA) & !is.na(nCount_RNA) &
                    nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5 & 
                    nCount_RNA > 500 & nCount_RNA < 25000)
    } else if (startsWith(prefix, "C")) {
    
      obj <- subset(obj, subset = !is.na(nFeature_RNA) & !is.na(nCount_RNA) &
                      nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5 & 
                      nCount_RNA > 500 & nCount_RNA < 12000)
    }
      
    
    # Generate post-filtering QC plot and feature relationship plot
    generateQCPlot(obj, prefix, suffix = "_AF")
    generateFeatureRelationshipPlot(obj, prefix, suffix = "_AF")
    
    # Normalize counts in each of the datasets
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj)
    obj <- ScaleData(obj)
    obj <- RunPCA(obj, npcs=30)
    obj <- FindNeighbors(obj, dims = 1:30)
    obj <- FindClusters(obj, resolution = 0.5) 
    obj <- RunUMAP(obj, reduction="pca", dims = 1:30)
    p <- DimPlot(obj, reduction="umap", group.by = "seurat_clusters")
    file_name <- paste0(prefix, "_UMAPplot_seurat_clusters.png")
    png(file_name , width = 10, height = 10, units = "in", res = 300)
    print(p)
    dev.off()
    # Save changes back to the list
    seurat_list[[i]] <- obj
  }
}

# Apply pre-proprocessing to each seurat gut_mergedect
preprocess_object(seurat_list_spine)
preprocess_object(seurat_list_gut)

# Merge spine controls and spine fructan
spine_control <- merge(seurat_list_spine[["B10_control"]], 
                       y = c(seurat_list_spine[["B11_control"]], seurat_list_spine[["B12_control"]]), 
                       add.cell.ids = c("B10", "B11", "B12"), project = "Spine_Control")

spine_fructan <- merge(seurat_list_spine[["B7_fructan"]], 
                       y = c(seurat_list_spine[["B8_fructan"]], seurat_list_spine[["B9_fructan"]]), 
                       add.cell.ids = c("B7", "B8", "B9"), project = "Spine_Fructan")

# Merge gut controls and gut fructan
gut_control <- merge(seurat_list_gut[["C4_control"]], 
                     y = c(seurat_list_gut[["C5_control"]], seurat_list_gut[["C6_control"]]), 
                     add.cell.ids = c("C4", "C5", "C6"), project = "Gut_Control")

gut_fructan <- merge(seurat_list_gut[["C1_fructan"]], 
                     y = c(seurat_list_gut[["C2_fructan"]], seurat_list_gut[["C3_fructan"]]), 
                     add.cell.ids = c("C1", "C2", "C3"), project = "Gut_Fructan")

# Assign identities to the merged gut_mergedects
Idents(spine_control) <- "Control"
Idents(spine_fructan) <- "B-fructan"
Idents(gut_control) <- "Control"
Idents(gut_fructan) <- "B-fructan"

# Merge spine gut_mergedects and gut gut_mergedects, adding identities
spine_merged <- merge(spine_control, y = spine_fructan, add.cell.ids = c("Control", "B-fructan"), project = "Spine")
gut_merged   <- merge(gut_control, y = gut_fructan, add.cell.ids = c("Control", "B-fructan"), project = "Gut")

# Add metadata to merged gut_mergedects
spine_merged$Diet <- Idents(spine_merged)
gut_merged$Diet <- Idents(gut_merged)

# Process spine merged
spine_merged <- NormalizeData(spine_merged)
spine_merged <- FindVariableFeatures(spine_merged)
spine_merged <- ScaleData(spine_merged)
spine_merged <- RunPCA(spine_merged, npcs = 30)
spine_merged <- FindNeighbors(spine_merged, dims = 1:30)
spine_merged <- FindClusters(spine_merged, resolution = 0.5) 
spine_merged <- RunUMAP(spine_merged, reduction = "pca", dims = 1:30)

# Process gut_merged
gut_merged <- NormalizeData(gut_merged)
gut_merged <- FindVariableFeatures(gut_merged)
gut_merged <- ScaleData(gut_merged)
gut_merged <- RunPCA(gut_merged, npcs = 30)
gut_merged <- FindNeighbors(gut_merged, dims = 1:30)
gut_merged <- FindClusters(gut_merged, resolution = 0.5) 
gut_merged <- RunUMAP(gut_merged, reduction = "pca", dims = 1:30)


plot_UMAP <- function(gut_merged, file_name, coloring){
  p <- DimPlot(gut_merged, reduction = "umap", group.by = coloring)
  png(file_name , width = 10, height = 10, units = "in", res = 300)
  print(p)
  dev.off()
}

plot_UMAP(spine_merged, "merged_datasets_spine_UMAP_by_orig_ident.png", "orig.ident")
plot_UMAP(spine_merged, "merged_datasets_spine_UMAP_by_seurat_clusters.png", "seurat_clusters")
plot_UMAP(spine_merged, "merged_datasets_spine_UMAP_by_Diet.png", "Diet")
plot_UMAP(gut_merged, "merged_datasets_gut_UMAP_by_orig_ident.png", "orig.ident")
plot_UMAP(gut_merged, "merged_datasets_gut_UMAP_by_Diet.png", "Diet")
plot_UMAP(gut_merged, "merged_datasets_gut_UMAP_by_seurat_clusters.png", "seurat_clusters")



# Exploration of marker genes
plot_markers <- function(obj, prefix){
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25,
                            logfc.threshold = 0.25)
  
  # Extract top 10 markers per cluster
  top_markers <- markers %>% group_by("seurat_clusters") %>% top_n(n = 10, wt = avg_log2FC)
  clusters    <- obj@meta.data$seurat_clusters 
  
  # Visualize top markers using DoHeatmap
  file_name <- paste0(prefix, "_markers_heatmap.png")
  png(file_name, width = 15, height = 10, units = "in", res = 300)
  p <- DoHeatmap(obj, features = top_markers$gene)
  print(p)
  dev.off()
  
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  top3_markers <- markers %>%
    group_by("seurat_clusters") %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 3) %>%
    ungroup()
  marker_genes_top3 <- top3_markers$gene
  
  # Create a unique list of marker genes
  unique_marker_genes <- unique(marker_genes_top3)
  
  # Create the DotPlot
  dotplot <- DotPlot(obj,
                     features = unique_marker_genes,
                     group.by = "seurat_clusters",
                     cols = c("lightgrey", "blue")) +  # Customize colors as needed
    RotatedAxis() +                       # Rotate x-axis labels for better readability
    ggtitle("DotPlot of Top 3 Marker Genes per Cluster") +
    theme(plot.title = element_text(hjust = 0.5))  # Center the title
  
  # Display the plot
  file_name <- paste0(prefix, "_markers_dotplot.png")
  print(file_name)
  png(file_name, width = 15, height = 10, units = 'in', res = 300)
  print(dotplot)
  dev.off()
  
}

spine_merged <- JoinLayers(spine_merged)
gut_merged <- JoinLayers(gut_merged)

plot_markers(spine_merged, "spine_merged")
plot_markers(gut_merged, "gut_merged")

# Save merged gut_mergedects for future analysis
saveRDS(spine_merged, file = "spine_merged_seurat.rds")
saveRDS(gut_merged, file = "gut_merged_seurat.rds")

