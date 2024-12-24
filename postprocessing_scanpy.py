import os
import scanpy as sc
import anndata as ad

# Path to the parent directory containing the sample directories
parent_dir = "/Users/juanjovel/OneDrive/jj/UofC/data_analysis/heatherArmstrong/pediatric_MS/sc_mice/matrices"

# List all subdirectories matching the naming pattern
sample_dirs = [os.path.join(parent_dir, d) for d in os.listdir(parent_dir) if d.endswith("_filtered_feature_bc_matrix")]

# Dictionary to store the AnnData objects with sample names as keys
adata_dict = {}

# Load each sample and make cell names unique by adding sample prefix
for sample_dir in sample_dirs:
    sample_name = os.path.basename(sample_dir).replace("_filtered_feature_bc_matrix", "")
    print(f"Loading data for {sample_name}...")
    adata = sc.read_10x_mtx(sample_dir, var_names="gene_symbols", cache=True)
    
    # Add sample prefix to cell names to ensure uniqueness across samples
    adata.obs_names = [f"{sample_name}_{cell}" for cell in adata.obs_names]
    adata.obs['sample'] = sample_name  # Add sample information
    adata_dict[sample_name] = adata

print("Now all samples will be merged")

# Merge all samples into a single AnnData object
adata_combined = sc.concat(adata_dict, axis=0, label="sample")

# Verify uniqueness in combined dataset
n_duplicates = sum(adata_combined.obs_names.duplicated())
if n_duplicates > 0:
    print(f"Warning: {n_duplicates} duplicate observations found in combined dataset")
    adata_combined.obs_names_make_unique()

# Determine some QC parameters
# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("mt-")

# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))

# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^Hb[^(P)]")

sc.pp.calculate_qc_metrics(
adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)

# Generate and save violin plot
output_file = "merged_violin_plot.png"
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
    save=output_file
)

output_file = "merged_scatter_plot.png"
sc.pl.scatter(
    adata, 
    "total_counts", 
    "n_genes_by_counts", 
    color="pct_counts_mt",
    save=output_file)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Run the doublets detection algorithm
sc.pp.scrublet(adata, batch_key="sample")

# Saving count data
adata.layers["counts"] = adata.X.copy()

# Normalizing to median total counts
sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)

# Feature selection
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")

sc.pl.highly_variable_genes(adata)

# Run PCA
sc.tl.pca(adata)

# Elbow plot
output_file = "merged_elbow_plot.png"
sc.pl.pca_variance_ratio(
    adata, 
    n_pcs=50, 
    log=True,
    save=output_file)

# Plot PCAs
output_file = "merged_PCAs_plot.png"
sc.pl.pca(
    adata,
    color=["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
    dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
    ncols=2,
    size=2,
    save=output_file
)

# Neirest neighbor graph construction and visualization
# Construct neighbors graph using PCA representation of data
sc.pp.neighbors(adata)

# Embed the graph into UMAP
sc.tl.umap(adata)

# Visualize UMAP according to sample
output_file = "merged_UMAP_plot.png"
sc.pl.umap(
    adata,
    color="sample",
    # Setting a smaller point size to get prevent overlap
    size=2,
    save=output_file
)

# Clustering
# As with Seurat and many other frameworks, we recommend the Leiden graph-clustering 
# method (community detection based on optimizing modularity) [Traag et al., 2019]. 
# Note that Leiden clustering directly clusters the neighborhood graph of cells, which 
# we already computed in the previous section.

# Using the igraph implementation and a fixed number of iterations can be significantly faster, especially for larger datasets
sc.tl.leiden(adata, flavor="igraph", n_iterations=2)

output_file = "merged_UMAP_afterClustering_plot.png"
sc.pl.umap(
    adata, 
    color=["leiden"],
    save=output_file)

# Re-assess quality control and cell filtering
# As indicated before, we will now re-assess our filtering strategy 
# by visualizing different QC metrics using UMAP.
output_file = "merged_UMAP_doublets_plot.png"
sc.pl.umap(
    adata,
    color=["leiden", "predicted_doublet", "doublet_score"],
    # increase horizontal space between panels
    wspace=0.5,
    size=3,
    save=output_file
)

# QC plots
output_file = "merged_UMAP_QC_plots.png"
sc.pl.umap(
    adata,
    color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
    wspace=0.5,
    ncols=2,
    save=output_file
)






# Save the combined AnnData object
adata_combined.write_h5ad("combined_samples.h5ad")
