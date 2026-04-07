import matplotlib.pyplot as plt
import geopandas as gpd

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
