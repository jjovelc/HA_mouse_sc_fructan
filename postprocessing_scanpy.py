import os
import scanpy as sc
import anndata as ad

# Path to the parent directory containing the sample directories
parent_dir = "/home/dayhoff/sdd/juan/projects/heatherArmstrong/SC_study/cellRangerResults/matrices"

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

    # Determine some QC parameters
    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith("mt-")

    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )

    # Generate and save violin plot
    output_file = f"{sample_name}_violin_plot.png"
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        save=output_file
    )


print("Now all samples will be merged")

# Merge all samples into a single AnnData object
adata_combined = sc.concat(adata_dict, axis=0, label="sample")

# Verify uniqueness in combined dataset
n_duplicates = sum(adata_combined.obs_names.duplicated())
if n_duplicates > 0:
    print(f"Warning: {n_duplicates} duplicate observations found in combined dataset")
    adata_combined.obs_names_make_unique()

# Save the combined AnnData object
adata_combined.write_h5ad("combined_samples.h5ad")
