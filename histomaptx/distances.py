import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Point
from shapely.ops import nearest_points
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.interpolate import interp1d
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from statsmodels.nonparametric.smoothers_lowess import lowess
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.ndimage.filters import gaussian_filter

def plot_cluster_expression_distance_heatmap(histomap, cluster_df, cluster_id, distance_col, bins=50, interpolate=False, max_distance=None, smooth=0):
    """
    Plots a heatmap of gene expression for a specific cluster over distance bins, with an optional max_distance cutoff.
    
    Args:
        histomap: histomap containing AnnData object containing gene expression data.
        cluster_df: DataFrame with 'Gene' and 'Cluster' columns.
        cluster_id: The cluster number to visualize.
        distance_col: Column name in adata.obs containing distance values.
        bins: Number of bins for the distance axis.
        interpolate: Whether to interpolate missing values (NaN) (default False).
        max_distance: Optional maximum distance threshold to include in the analysis.
    """
    # Get genes belonging to the selected cluster
    cluster_genes = cluster_df[cluster_df["Cluster"] == cluster_id]["Gene"].tolist()
    if not cluster_genes:
        raise ValueError(f"No genes found for cluster {cluster_id}.")
    adata = histomap.to_anndata()
    # Extract distance values
    distances = adata.obs[distance_col]
    
    # Apply max_distance filter if provided
    if max_distance is not None:
        mask = distances <= max_distance
        adata = adata[mask]  # Subset AnnData according to max_distance
        distances = distances[mask]  # Filter distances
    
    # Define bins and bin centers
    bin_edges = np.linspace(distances.min(), distances.max(), bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Compute average gene expression per bin for each gene

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    heatmap_data = []
    for gene in cluster_genes:
        if gene not in adata.var_names:
            continue  # Skip missing genes

        gene_expression = adata[:, gene].X.toarray().flatten() if hasattr(adata[:, gene].X, "toarray") else np.array(adata[:, gene].X).flatten()
        
        # Remove NaN values
        mask = ~np.isnan(gene_expression) & ~np.isnan(distances)
        gene_expression, distances_filtered = gene_expression[mask], distances[mask]

        # Compute binned averages
        bin_indices = np.digitize(distances_filtered, bin_edges) - 1
        avg_expression = [
            gene_expression[bin_indices == i].mean() if np.any(bin_indices == i) else np.nan
            for i in range(bins)
        ]
        
        if interpolate:
            # Interpolate to fill missing values (NaN) with linear interpolation
            avg_expression = np.array(avg_expression)
            nan_mask = np.isnan(avg_expression)
            
            # Perform linear interpolation to fill NaNs
            interp = interp1d(bin_centers[~nan_mask], avg_expression[~nan_mask], kind='linear', fill_value="extrapolate")
            avg_expression[nan_mask] = interp(bin_centers[nan_mask])
        
        heatmap_data.append(avg_expression)
    
    # Convert to DataFrame for seaborn heatmap
    heatmap_df = pd.DataFrame(heatmap_data, index=cluster_genes, columns=bin_centers)
    smoothed_values = gaussian_filter(heatmap_df.values, sigma=smooth)

    # Create a new DataFrame with the same index and columns
    heatmap_df_smoothed = pd.DataFrame(smoothed_values, index=heatmap_df.index, columns=heatmap_df.columns)

    # Plot heatmap
    plt.figure(figsize=(10, len(cluster_genes) * 0.5 + 3))
    sns.clustermap(heatmap_df_smoothed, cmap="viridis", xticklabels=10, yticklabels=True, cbar=True, col_cluster=False)
    plt.xlabel("Distance to Annotation Spot")
    plt.ylabel("Gene")
    plt.title(f"Expression Heatmap - Cluster {cluster_id}")
    plt.show()

def compute_loess_trends(adata, genes, distance_col, frac=0.3, grid_size=100, max_distance=None):
    """
    Computes LOESS-smoothed gene expression trends over distance, with an optional max distance cutoff.
    The expression values are standardized before LOESS smoothing, unless the standard deviation is zero.

    Args:
        adata: AnnData object containing gene expression data.
        genes: List of genes to analyze.
        distance_col: Column name in adata.obs containing distance values.
        frac: Smoothing parameter for LOESS (0 to 1, higher means more smoothing).
        grid_size: Number of points to interpolate for alignment.
        max_distance: Optional maximum distance threshold to include in the analysis.

    Returns:
        DataFrame with genes as rows and interpolated LOESS values as columns.
    """
    distances = adata.obs[distance_col].values

    # Apply distance filter if max_distance is provided
    if max_distance is not None:
        valid_idx = distances <= max_distance
        distances = distances[valid_idx]
        adata = adata[valid_idx]  # Subset AnnData accordingly

    common_grid = np.linspace(distances.min(), distances.max(), grid_size)
    trends = {}

    for gene in genes:
        if gene not in adata.var_names:
            continue  # Skip missing genes
        
        # Get expression values for the gene
        expr_values = adata[:, gene].X.toarray().flatten() if hasattr(adata[:, gene].X, "toarray") else np.array(adata[:, gene].X).flatten()
        distances_filtered = adata.obs[distance_col].values  # Use filtered distances

        # Remove NaN values
        mask = ~np.isnan(expr_values) & ~np.isnan(distances_filtered)
        expr_values, distances_filtered = expr_values[mask], distances_filtered[mask]

        # Check if there is any data left after filtering
        if len(expr_values) == 0 or len(distances_filtered) == 0:
            print(f"Warning: No data left for gene {gene} after filtering.")
            continue  # Skip this gene if no data left

        # Standardize the gene expression (subtract mean, divide by standard deviation)
        if np.std(expr_values) > 0:
            expr_values_standardized = (expr_values - np.mean(expr_values)) / np.std(expr_values)
        else:
            expr_values_standardized = expr_values  # Skip standardization if standard deviation is zero

        # Compute LOESS smoothing on standardized data
        smoothed = lowess(expr_values_standardized, distances_filtered, frac=frac, return_sorted=True)
        
        # Debugging: Check if smoothed values are empty
        if smoothed.shape[0] == 0:
            print(f"Warning: LOESS smoothing failed for gene {gene}.")
            continue

        # Interpolate LOESS-smoothed values on the common grid
        try:
            interp_values = np.interp(common_grid, smoothed[:, 0], smoothed[:, 1])  # Interpolate on grid
            trends[gene] = interp_values
        except ValueError as e:
            print(f"Error: Interpolation failed for gene {gene}. Details: {e}")
            continue  # Skip this gene if interpolation fails

    return pd.DataFrame(trends, index=common_grid)


def cluster_genes(trend_df, method="ward", num_clusters=3, min_obs=2):
    """
    Clusters genes based on LOESS-smoothed trends, independent of expression levels.
    Ensures that no cluster contains fewer than 'min_obs' genes.
    
    Args:
        trend_df: DataFrame with genes as columns and LOESS values as rows.
        method: Clustering method (default: 'ward').
        num_clusters: Number of clusters.
        min_obs: Minimum number of genes required in each cluster.
    
    Returns:
        DataFrame with gene cluster assignments, ensuring minimum observations per cluster.
    """
    # Perform hierarchical clustering
    linkage_matrix = linkage(trend_df.T, method=method)
    
    # Assign clusters
    clusters = fcluster(linkage_matrix, num_clusters, criterion='maxclust')
    
    # Create cluster assignments DataFrame
    cluster_df = pd.DataFrame({"Gene": trend_df.columns, "Cluster": clusters})
    
    # Ensure minimum observations per cluster
    cluster_sizes = cluster_df.groupby("Cluster").size()
    small_clusters = cluster_sizes[cluster_sizes < min_obs].index.tolist()
    
    # Handle small clusters (e.g., merge them with larger clusters)
    if small_clusters:
        for small_cluster in small_clusters:
            # Find the largest cluster to merge with
            larger_cluster = cluster_sizes[cluster_sizes > min_obs].idxmax()
            
            # Reassign genes in the small cluster to the larger cluster
            cluster_df.loc[cluster_df["Cluster"] == small_cluster, "Cluster"] = larger_cluster
    
    return cluster_df, linkage_matrix


def plot_cluster_avg_profiles(trend_df, cluster_df):
    """
    Plots the average LOESS profile for each gene cluster.

    Args:
        trend_df: DataFrame where rows are distances and columns are genes.
        cluster_df: DataFrame with 'Gene' and 'Cluster' columns.
    """
    trend_df_T = trend_df.T  # Transpose to align with cluster_df
    trend_df_T["Cluster"] = cluster_df.set_index("Gene")["Cluster"]

    # Compute mean LOESS profile for each cluster
    cluster_means = trend_df_T.groupby("Cluster").mean().T

    # Plot
    plt.figure(figsize=(8, 6))
    for cluster in cluster_means.columns:
        plt.plot(trend_df.index, cluster_means[cluster], label=f"Cluster {cluster}")

    plt.xlabel("Distance to Annotation Spot")
    plt.ylabel("Average Smoothed Gene Expression")
    plt.title("Average LOESS Profiles by Cluster")
    plt.legend()
    plt.show()

def identify_gene_distance_trend(histomap, annotation, ncluster=3, grid_size=100, max_distance=False, n_top_genes=1000, frac=0.3, plot_trend_per_cluster=False):
    """
    Identify variable genes based on annotation, compute LOESS trends for genes according to distance,
    cluster the genes based on smoothed expression profiles, and return a clustered average expression profile plot.
    
    Args:
        histomap: HistoMap object containing gene expression data.
        annotation: The annotation (e.g., "Tumor") to calculate distance from and to subset.
        ncluster: Number of clusters to form based on LOESS smoothed expression profiles.
        max_distance: Optional maximum distance threshold to include in the analysis.
        frac: The fraction for LOESS smoothing.
        plot_trend_per_cluster: Boolean flag to plot average trends per cluster.

    Returns:
        cluster_df: DataFrame with gene cluster assignments based on trend.
    """
    # Convert HistoMap object to AnnData
    adata = histomap.to_anndata()
    # Construct the positive and negative column names based on annotation
    positive_col = annotation + "_positive"
    
    # Validate that the annotation column exists
    if positive_col not in adata.obs.columns:
        raise ValueError(f"Column '{positive_col}' not found in adata.obs.")
    
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata_subset = adata[~adata.obs["distance_to_" + annotation].isna(), :]

    print('Get variable genes.')
    sc.pp.highly_variable_genes(adata_subset, n_top_genes=n_top_genes)
    highly_variable_genes = adata_subset.var[adata_subset.var['highly_variable']].index.tolist()
    print('Compute trends according to the distance for the variable genes')
    
    trend_df = compute_loess_trends(adata_subset, highly_variable_genes, "distance_to_" + annotation, max_distance=max_distance)

    cluster_df, linkage_matrix = cluster_genes(trend_df, num_clusters=ncluster)
    cluster_df['Max_distance'] = max_distance
    if plot_trend_per_cluster:
        plot_cluster_avg_profiles(trend_df, cluster_df)
    
    return cluster_df


def compute_nearest_annotation_distance(self, annotation):
    """
    For each spot that does not have the given annotation, compute the distance to the closest spot with the annotation.
    
    Args:
        self: HistoMap object containing spot_geodata.
        annotation: Annotation type to use as reference (e.g., "Tumor").
    """
    # Filter spots with the given annotation
    annotated_spots = self.spot_geodata[self.spot_geodata['Annotation_map'] == annotation]
    other_spots = self.spot_geodata[self.spot_geodata['Annotation_map'] != annotation]
    
    # Ensure there are annotated spots to compute distance from
    if annotated_spots.empty:
        raise ValueError(f"No spots found with annotation '{annotation}'.")
    
    # Extract geometries
    annotated_points = annotated_spots.geometry.centroid
    other_points = other_spots.geometry.centroid
    
    # Compute nearest distances
    distances = []
    for point in other_points:
        nearest_point = min(annotated_points, key=lambda p: point.distance(p))
        distances.append(point.distance(nearest_point))
    
    # Assign distances to the non-annotated spots
    self.spot_geodata.loc[other_spots.index, f'distance_to_{annotation}'] = distances

def plot_distance_overlay(histomap, annotation, display_image=True, max_cutoff=None):
    """
    Plots the distances to the closest spot with the given annotation.
    
    Args:
        histomap: HistoMap object containing spot_geodata.
        annotation: Annotation type to compute distances to.
        display_image: Whether to display the underlying histology image.
        max_cutoff: Optional max cutoff for the color scale.
    """
    gdf = histomap.spot_geodata.copy()
    distance_col = f'distance_to_{annotation}'
    
    if distance_col not in gdf.columns:
        raise ValueError(f"Column '{distance_col}' not found. Run compute_nearest_annotation_distance() first.")
    
    # Convert distance column to float if necessary
    gdf[distance_col] = gdf[distance_col].astype(float)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Display the image beneath the annotations
    if display_image and hasattr(histomap, 'plotting_image') and histomap.plotting_image is not None:
        extent = [0, histomap.full_res_width, histomap.full_res_height, 0]
        ax.imshow(histomap.plotting_image.values.transpose(1, 2, 0), extent=extent, origin='upper', cmap='gray')
    
    # Define color scale limits
    vmin = gdf[distance_col].min()
    vmax = max_cutoff if max_cutoff is not None else gdf[distance_col].max()
    
    # Plot the polygons with a colormap based on distance
    gdf.plot(ax=ax, column=distance_col, cmap='plasma', legend=True,
             legend_kwds={'label': f"Distance to {annotation} Spot",
                          'orientation': "vertical"},
             vmin=vmin, vmax=vmax)
    
    # Set title and axis labels
    ax.set_title(f'Distance to {annotation}', fontsize=15)
    ax.set_xlabel('X Coordinate', fontsize=12)
    ax.set_ylabel('Y Coordinate', fontsize=12)
    
    plt.tight_layout()
    plt.show()

def calculate_optimal_bins(distances, rule='sturges'):
    """
    Calculates the optimal number of bins for a histogram using Sturges' Rule or Freedman-Diaconis Rule.
    
    Args:
        distances: The array of distances to calculate bins for.
        rule: The rule to use ('sturges' or 'freedman-diaconis').
    
    Returns:
        The optimal number of bins.
    """
    n = len(distances)
    
    if rule == 'sturges':
        # Sturges' Rule
        return int(np.ceil(np.log2(n) + 1))
    
    elif rule == 'freedman-diaconis':
        # Freedman-Diaconis Rule
        iqr = np.percentile(distances, 75) - np.percentile(distances, 25)
        bin_width = 2 * iqr * n ** (-1 / 3)
        return int(np.ceil((distances.max() - distances.min()) / bin_width))
    
    else:
        raise ValueError("Unknown rule. Choose either 'sturges' or 'freedman-diaconis'.")

def plot_gene_expression_distance_heatmap(histomap, genes, annotation, bins=None, interpolate=True, normalize=True, log1p=True, bin_rule='freedman-diaconis', cmap="rocket", cluster_gene=True):
    """
    Plots gene expression as a heatmap based on distance.
    Each row represents a gene, with average expression across distance bins.
    
    Args:
        histomap: HistoMap object containing gene expression data.
        genes: List of gene names to plot.
        annotation: The annotation (e.g., "Tumor") to calculate distance from.
        bins: Number of bins for the distance axis. If None, automatically determined.
        interpolate: Whether to interpolate missing values (default True).
        normalize: Whether to normalize the expression data (default True).
        log1p: Whether to apply log1p transformation to expression data (default True).
        bin_rule: Rule for determining bin size ('sturges' or 'freedman-diaconis').
    """
    # Convert HistoMap object to AnnData
    adata = histomap.to_anndata()

    # Construct the distance column name based on annotation
    distance_col = "distance_to_" + annotation

    # Validate that the distance column exists in adata
    if distance_col not in adata.obs.columns:
        raise ValueError(f"Column '{distance_col}' not found in adata.obs.")

    # Apply normalization if requested
    if normalize:
        sc.pp.normalize_total(adata)

    # Apply log1p transformation if requested
    if log1p:
        sc.pp.log1p(adata)

    # Check if all genes are present in adata
    for gene in genes:
        if gene not in adata.var_names:
            raise ValueError(f"Gene '{gene}' not found in AnnData object.")

    # Extract distance values
    distances = adata.obs[distance_col]

    # Determine bin size if not provided
    if bins is None:
        bins = calculate_optimal_bins(distances, rule=bin_rule)

    # Define bins and bin centers
    bin_edges = np.linspace(distances.min(), distances.max(), bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Compute average gene expression per bin for each gene
    heatmap_data = []
    for gene in genes:
        gene_expression = adata[:, gene].X.toarray().flatten() if hasattr(adata[:, gene].X, "toarray") else np.array(adata[:, gene].X).flatten()

        # Remove NaN values
        mask = ~np.isnan(gene_expression) & ~np.isnan(distances)
        gene_expression = gene_expression[mask]
        distances_filtered = distances[mask]

        # Compute binned averages
        bin_indices = np.digitize(distances_filtered, bin_edges) - 1
        avg_expression = [
            gene_expression[bin_indices == i].mean() if np.any(bin_indices == i) else np.nan
            for i in range(bins)
        ]

        if interpolate:
            # Interpolate to fill missing values (NaN)
            avg_expression = np.array(avg_expression)
            nan_mask = np.isnan(avg_expression)

            # Perform linear interpolation to fill NaNs
            interp = interp1d(bin_centers[~nan_mask], avg_expression[~nan_mask], kind='linear', fill_value="extrapolate")
            avg_expression[nan_mask] = interp(bin_centers[nan_mask])

        heatmap_data.append(avg_expression)

    # Convert to DataFrame for seaborn heatmap
    heatmap_df = pd.DataFrame(heatmap_data, index=genes, columns=bin_centers)

    # Plot heatmap
    plt.figure(figsize=(10, len(genes) * 0.5 + 3))
    sns.clustermap(heatmap_df, cmap=cmap,  yticklabels=True, cbar=True, row_cluster=cluster_gene, col_cluster=False, dendrogram_ratio=0.1)
    plt.xlabel("Distance to Annotation Spot")
    plt.ylabel("Gene")
    plt.title(f"Gene Expression Heatmap vs. Distance - Annotation: {annotation}")
    plt.show()
