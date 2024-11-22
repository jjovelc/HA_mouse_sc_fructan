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

# Initialize lists to store Seurat objects for spine and gut samples
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
    # Access the Seurat object using the index
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
    obj <- subset(obj, subset = !is.na(nFeature_RNA) & !is.na(nCount_RNA) &
                    nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5 & 
                    nCount_RNA > 500 & nCount_RNA < 10000)
    
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

# Apply pre-proprocessing to each seurat object
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

# Assign identities to the merged objects
Idents(spine_control) <- "Control"
Idents(spine_fructan) <- "B-fructan"
Idents(gut_control) <- "Control"
Idents(gut_fructan) <- "B-fructan"

# Merge spine objects and gut objects, adding identities
spine_merged <- merge(spine_control, y = spine_fructan, add.cell.ids = c("Control", "B-fructan"), project = "Spine")
gut_merged <- merge(gut_control, y = gut_fructan, add.cell.ids = c("Control", "B-fructan"), project = "Gut")

# Add metadata to merged objects
spine_merged$Diet <- Idents(spine_merged)
gut_merged$Diet <- Idents(gut_merged)


process_and_plot <- function(obj, file_name, coloring){
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, npcs = 30)
  obj <- FindNeighbors(obj, dims = 1:30)
  obj <- FindClusters(obj, resolution = 0.5) 
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)
  p <- DimPlot(obj, reduction = "umap", group.by = coloring)
  png(file_name , width = 10, height = 10, units = "in", res = 300)
  print(p)
  dev.off()
}

process_and_plot(spine_merged, "merged_datasets_spine_UMAP_by_ident.png", "ident")
process_and_plot(spine_merged, "merged_datasets_spine_UMAP_by_seurat_clusters.png", "seurat_clusters")
process_and_plot(gut_merged, "merged_datasets_gut_UMAP_by_ident.png", "ident")
process_and_plot(gut_merged, "merged_datasets_gut_UMAP_by_seurat_clusters.png", "seurat_clusters")

# Save merged objects for future analysis
saveRDS(spine_merged, file = "spine_merged_seurat.rds")
saveRDS(gut_merged, file = "gut_merged_seurat.rds")

