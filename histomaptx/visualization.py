import geopandas as gpd
import pandas as pd
import ast
import gzip
import io
import zipfile
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import warnings
import spatialdata
import spatialdata_io
import geopandas as gpd
from shapely.geometry import Polygon, Point
from rtree import index
import numpy as np
from . import histomap_utils
import seaborn as sns 

from spatialdata import bounding_box_query
import spatialdata_plot
from matplotlib.patches import PathPatch
from matplotlib.path import Path
import numpy as np

def polygon_to_path(polygon):
    """Convert a Shapely Polygon (with holes) to a Matplotlib Path."""
    vertices = []
    codes = []

    # Exterior
    exterior_coords = np.asarray(polygon.exterior.coords)
    vertices.extend(exterior_coords)
    codes.extend([Path.MOVETO] + [Path.LINETO] * (len(exterior_coords) - 2) + [Path.CLOSEPOLY])

    # Interiors (holes)
    for interior in polygon.interiors:
        interior_coords = np.asarray(interior.coords)
        vertices.extend(interior_coords)
        codes.extend([Path.MOVETO] + [Path.LINETO] * (len(interior_coords) - 2) + [Path.CLOSEPOLY])

    return Path(vertices, codes)



def plot_cells(histomap, fill=True, contour='k', display_image=True):
    """
    Plots segmented cells based on a GeoDataFrame and optionally displays an image beneath the cells.

    Parameters:
    - histomap: HistoMap object containing the segmentation dataframe and spot geodata.
    - fill: False, True, or a list of colors.
        - False: No fill (only contours).
        - True: Uses a default colormap ('tab20') for fill.
        - List of colors: Specifies the fill color for each cell.
    - contour: Color or list of colors for the contours. If None, the default colormap is used.
    - display_image: Boolean, whether to display the image beneath the cells (True or False).
    """
    # Check if 'n_cell' exists in spot_geodata
    if 'n_cell' not in histomap.spot_geodata.columns:
        raise ValueError("'n_cell' column does not exist in spot_geodata. Please ensure segmentation has been added correctly.")
    
    fig, ax = plt.subplots(figsize=(10, 10))
    gdf = histomap.segmentation_dataframe
    
    # Default colormap
    cmap = plt.colormaps.get_cmap('tab20')
    unique_cells = len(gdf)

    # Handle fill colors
    if fill is True:
        fill = [cmap(idx % 20) for idx in range(unique_cells)]
    elif isinstance(fill, list):
        if len(fill) != unique_cells:
            raise ValueError("The length of the 'fill' list must match the number of cells.")
    elif fill is False:
        fill = [None] * unique_cells  # Ensure correct indexing

    # Handle contour colors
    if contour is None:
        contour = [cmap(idx % 20) for idx in range(unique_cells)]
    elif isinstance(contour, str):
        contour = [contour] * unique_cells  # Convert single color to list
    elif isinstance(contour, list) and len(contour) != unique_cells:
        raise ValueError("The length of the 'contour' list must match the number of cells.")

    # Display the image beneath the cells (if requested and available)
    if display_image and hasattr(histomap, 'plotting_image') and histomap.plotting_image is not None:
        extent = [0, histomap.full_res_width, histomap.full_res_height, 0]
        ax.imshow(histomap.plotting_image.values.transpose(1, 2, 0), extent=extent, origin='upper', cmap='gray')

    # Plot cells
    for idx, (geom, fill_color, contour_color) in enumerate(zip(gdf.geometry, fill, contour)):
        if geom.is_empty or not geom.is_valid:
            continue

        if geom.geom_type == 'Polygon':
            x, y = geom.exterior.xy
            if fill_color:
                ax.fill(x, y, color=fill_color, edgecolor=contour_color)
            else:
                ax.plot(x, y, color=contour_color)
        elif geom.geom_type == 'MultiPolygon':
            for polygon in geom.geoms:
                x, y = polygon.exterior.xy
                if fill_color:
                    ax.fill(x, y, color=fill_color, edgecolor=contour_color)
                else:
                    ax.plot(x, y, color=contour_color)

    ax.set_title('Segmented Cells')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_aspect('equal', adjustable='box')

    # Adjust limits
    bounds = gdf.total_bounds
    ax.set_xlim(bounds[0], bounds[2])
    ax.set_ylim(bounds[1], bounds[3])

    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.show()


def plot_combined_annotation_overlap(histomap, annotation1, annotation2):
    """
    Identifies spots that overlap with both specified annotations and plots them.
    Spots overlapping both annotations are colored red, all other spots are colored blue.
    
    Parameters:
    - annotation1: First annotation to check for overlap
    - annotation2: Second annotation to check for overlap
    
    Returns:
    - A matplotlib figure showing the spots with appropriate coloring
    """

    # Check if the annotations exist in the data
    available_annotations = histomap.data_exploded['Annotation'].unique()
    if annotation1 not in available_annotations or annotation2 not in available_annotations:
        raise ValueError(f"One or both annotations not found. Available annotations: {', '.join(available_annotations)}")
    
    # Create a copy of the spot data
    gdf = histomap.spot_geodata.copy()
    
    # Create a new column indicating whether a spot overlaps with both annotations
    overlap_col1 = str(annotation1) + '_overlap'
    overlap_col2 = str(annotation2) + '_overlap'
    
    # Check if the overlap columns exist
    if overlap_col1 not in gdf.columns or overlap_col2 not in gdf.columns:
        raise ValueError(f"Overlap data for one or both annotations not found. Please ensure compute_overlap_annotation() was called.")
    
    # Create a boolean column indicating spots that overlap with both annotations
    gdf['dual_overlap'] = (gdf[overlap_col1] > 0) & (gdf[overlap_col1] < 100) & (gdf[overlap_col2] > 0)& (gdf[overlap_col2] < 100)
    
    # Plot the results
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Plot spots that don't overlap with both annotations in blue
    non_overlapping = gdf[~gdf['dual_overlap']]
    if not non_overlapping.empty:
        non_overlapping.plot(ax=ax, color='blue', alpha=0.5, label='Non-overlapping spots')
    
    # Plot spots that overlap with both annotations in red
    overlapping = gdf[gdf['dual_overlap']]
    if not overlapping.empty:
        overlapping.plot(ax=ax, color='red', alpha=0.7, label=f'Spots overlapping {annotation1} & {annotation2}')
    
    # Add legend, title and labels
    ax.legend()
    ax.set_title(f'Spots Overlapping Both {annotation1} and {annotation2}', fontsize=15)
    ax.set_xlabel('X Coordinate', fontsize=12)
    ax.set_ylabel('Y Coordinate', fontsize=12)
    
    # Add annotation statistics
    spot_count = len(gdf)
    overlap_count = len(overlapping)
    overlap_percentage = (overlap_count / spot_count) * 100 if spot_count > 0 else 0
    
    stats_text = (
        f"Total spots: {spot_count}\n"
        f"Spots overlapping both annotations: {overlap_count} ({overlap_percentage:.2f}%)"
    )
    
    # Place the statistics text box in the upper right corner
    ax.text(0.98, 0.98, stats_text, 
            transform=ax.transAxes, 
            horizontalalignment='right',
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.show()
    
    return fig, ax, overlapping

def plot_annotation_overlay(
    histomap,
    annotation,
    resolution="lowres",
    disabled_annotations=None,
    xcoords=None,
    ycoords=None,
    figsize=(10, 10),
    display_image=True,
    save=None,
    cmap=None,
):

    if isinstance(annotation, str):
        annotation = [annotation]

    if disabled_annotations is None:
        disabled_annotations = []

    annotation = [ann for ann in annotation if ann not in disabled_annotations]

    if len(annotation) == 0:
        raise ValueError("No annotations to plot after removing disabled annotations.")

    # ---- check overlap columns ----
    missing = [
        ann + "_overlap"
        for ann in annotation
        if ann + "_overlap" not in histomap.spot_geodata.columns
    ]

    if missing:
        raise ValueError(
            f"Missing columns {missing}. Run compute_overlap_annotation() first."
        )

    for ann in annotation:

        overlap_col = ann + "_overlap"

        if histomap.visium_type == "visium_hd":

            sdata = histomap.visium_spatialdata
            bin_str = f"{int(histomap.bin_size):03d}um"
            table_key = next(k for k in sdata.tables.keys() if bin_str in k)

            if overlap_col not in sdata.tables[table_key].obs:
                sdata.tables[table_key].obs[overlap_col] = (
                    histomap.spot_geodata[overlap_col].values
                )

            plot_visium_hd(
                histomap,
                resolution=resolution,
                xcoords=xcoords,
                ycoords=ycoords,
                annotation_key=overlap_col,
                figsize=figsize,
                display_image=display_image,
                save=save,
                cmap=cmap,
            )

        elif histomap.visium_type == "visium":

            sdata = histomap.visium_spatialdata
            table_key = list(sdata.tables.keys())[0]

            if overlap_col not in sdata.tables[table_key].obs:
                sdata.tables[table_key].obs[overlap_col] = (
                    histomap.spot_geodata[overlap_col].values
                )

            plot_visium(
                sdata,
                resolution=resolution,
                xcoords=xcoords,
                ycoords=ycoords,
                annotation_key=overlap_col,
                figsize=figsize,
                display_image=display_image,
                save=save,
                cmap=cmap,
            )

        else:
            raise ValueError(
                f"Unsupported visium_type: {histomap.visium_type}"
            )

def plot_annotation_overlay_old(
    histomap,
    annotation,
    disabled_annotations=None,
    display_image=True,
    max_cutoff=None,
    save=None,
    raster_threshold=5000,  
    dpi=300,
    cmap="viridis",
    **kwargs
):
    """
    Optimized plotting for large Visium / Visium HD datasets.
    Automatically switches to raster mode when number of spots
    exceeds `raster_threshold`.
    """

    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from rasterio.features import rasterize

    if isinstance(annotation, str):
        annotation = [annotation]

    if disabled_annotations is None:
        disabled_annotations = []

    annotation = [ann for ann in annotation if ann not in disabled_annotations]

    if len(annotation) == 0:
        raise ValueError("No annotations to plot after removing disabled annotations.")

    missing_columns = [
        ann + "_overlap"
        for ann in annotation
        if ann + "_overlap" not in histomap.spot_geodata.columns
    ]

    if missing_columns:
        raise ValueError(
            f"Missing columns {missing_columns}. "
            "Run compute_overlap_annotation() first."
        )

    num_annotations = len(annotation)
    fig, axes = plt.subplots(1, num_annotations,
                             figsize=(5 * num_annotations, 8))

    if num_annotations == 1:
        axes = [axes]

    for ax, ann in zip(axes, annotation):

        gdf = histomap.spot_geodata.copy()
        overlap_col = ann + "_overlap"
        gdf[overlap_col] = gdf[overlap_col].astype(float)

        # ---- BACKGROUND IMAGE ----
        if (
            display_image
            and hasattr(histomap, "plotting_image")
            and histomap.plotting_image is not None
        ):
            extent = [0, histomap.full_res_width,
                      histomap.full_res_height, 0]

            ax.imshow(
                histomap.plotting_image.values.transpose(1, 2, 0),
                extent=extent,
                origin="upper",
                cmap="gray"
            )

        vmin = gdf[overlap_col].min()
        vmax = max_cutoff if max_cutoff is not None else gdf[overlap_col].max()

        # ================================
        # AUTO-SWITCH: VECTOR vs RASTER
        # ================================
        if len(gdf) > raster_threshold:

            print(f"Large dataset detected ({len(gdf)} spots). Using raster mode.")

            shapes = (
                (geom, value)
                for geom, value in zip(gdf.geometry, gdf[overlap_col])
            )

            raster = rasterize(
                shapes=shapes,
                out_shape=(
                    int(histomap.full_res_height),
                    int(histomap.full_res_width)
                ),
                fill=np.nan,
                dtype="float32"
            )

            im = ax.imshow(
                raster,
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                origin="upper"
            )

            cbar = plt.colorbar(
                im,
                ax=ax,
                orientation="horizontal",
                pad=0.05
            )
            cbar.set_label(f"{ann} Positive Overlap (%)")

        else:
            # Small dataset → normal GeoPandas plot
            gdf.plot(
                ax=ax,
                column=overlap_col,
                legend=True,
                legend_kwds={
                    'label': f"{ann} Positive Overlap (%)",
                    'orientation': "horizontal"
                },
                vmin=vmin,
                vmax=vmax,
                rasterized=True,  # helps for PDFs
                **kwargs
            )

        ax.set_title(f"{ann} Positive Overlap", fontsize=14)
        ax.set_xlabel("X Coordinate")
        ax.set_ylabel("Y Coordinate")

    plt.tight_layout()

    # ---- SAVE ----
    if save is not None:
        _, ext = os.path.splitext(save)
        if ext.lower() not in [".png", ".pdf", ".svg", ".jpg", ".jpeg", ".tiff"]:
            raise ValueError("Unsupported file extension.")
        fig.savefig(save, bbox_inches="tight", dpi=dpi)
        print(f"Figure saved to: {save}")

    plt.show()




def plot_annotations(histomap, fill=False, contour=None, annotation=None,
                     display_image=False, alpha=1, save=None):
    """Plots the annotations based on the DataFrame, respecting the plot order.
    
    Parameters:
    - histomap: HistoMap object containing annotation data.
    - fill: False, True, or a list of colors. 
            False means no fill (only contours), 
            True uses colors from histomap.annotation_colors, 
            a list of colors specifies the fill color for each annotation.
    - contour: a color or list of colors for the contours. 
               If None, uses colors from histomap.annotation_colors or default colormap.
    - annotation: a specific annotation (string) or a list of annotations to plot. If None, all annotations are plotted.
    - display_image: bool, whether to display `histomap.plotting_image` beneath the annotations.
    - alpha: float between 0 and 1, transparency of the fill color (0 is completely transparent, 1 is opaque).
             Only applies when fill is not False.
    - save: str or None. File path to save the figure. Format inferred from extension
            (e.g., .png, .pdf). If None, the figure is not saved.
    """
    import os

    # Validate alpha parameter
    if not (0 <= alpha <= 1):
        raise ValueError("Alpha must be between 0 and 1")
    
    fig, ax = plt.subplots(figsize=(10, 10))

    # Display the image beneath the annotations
    if display_image and hasattr(histomap, 'plotting_image') and histomap.plotting_image is not None:
        extent = [0, histomap.full_res_width, histomap.full_res_height, 0]
        ax.imshow(histomap.plotting_image.values.transpose(1, 2, 0),
                  extent=extent, origin='upper', cmap='gray')

    cmap = plt.cm.get_cmap('tab20')
    unique_annotations = list(histomap.data_exploded['Annotation'].unique())

    # Validate annotation input
    if annotation is not None:
        if isinstance(annotation, str):
            annotation = [annotation]
        elif not isinstance(annotation, list) or not all(isinstance(a, str) for a in annotation):
            raise TypeError("Annotation must be a string, a list of strings, or None.")

        missing_annotations = [ann for ann in annotation if ann not in unique_annotations]
        if missing_annotations:
            raise ValueError(
                f"Annotations {missing_annotations} not found. "
                f"Available annotations are: {sorted(unique_annotations)}"
            )
        annotations_to_plot = annotation
    else:
        annotations_to_plot = unique_annotations

    annotations_to_plot = [
        ann for ann in annotations_to_plot
        if ann not in histomap.disabled_annotations
    ]

    if histomap.disabled_annotations:
        print(f"Skipping disabled annotations: {', '.join(histomap.disabled_annotations)}")

    if len(histomap.activated_annotations) == 0:
        raise ValueError("No activated annotations to compute overlap.")

    # Process fill colors
    if fill is True and hasattr(histomap, 'annotation_colors'):
        fill_dict = dict(zip(histomap.annotation_colors['annotation'],
                             histomap.annotation_colors['color']))
        fill = [fill_dict.get(ann, cmap(i % 20))
                for i, ann in enumerate(annotations_to_plot)]
    elif fill is True:
        fill = [cmap(idx % 20) for idx in range(len(annotations_to_plot))]
    elif isinstance(fill, list):
        if len(fill) != len(annotations_to_plot):
            raise ValueError(
                "The length of the 'fill' list must match the number of unique annotations."
            )
    elif fill is False:
        fill = [None] * len(annotations_to_plot)

    # Process contour colors
    if contour is None and hasattr(histomap, 'annotation_colors'):
        contour_dict = dict(zip(histomap.annotation_colors['annotation'],
                                histomap.annotation_colors['color']))
        contour = [contour_dict.get(ann, cmap(i % 20))
                   for i, ann in enumerate(annotations_to_plot)]
    elif contour is None:
        contour = [cmap(idx % 20) for idx in range(len(annotations_to_plot))]
    elif isinstance(contour, list):
        if len(contour) != len(annotations_to_plot):
            raise ValueError(
                "The length of the 'contour' list must match the number of unique annotations."
            )
    elif isinstance(contour, str):
        contour = [contour] * len(annotations_to_plot)

    # Sort annotations based on plot_order (0 on top)
    data_sorted = histomap.data_exploded.sort_values(by="plot_order", ascending=False)

    # Plot annotations
    for _, row in data_sorted.iterrows():
        geom = row['geometry']
        ann = row['Annotation']

        if ann not in annotations_to_plot or geom.is_empty or not geom.is_valid:
            continue

        ann_idx = annotations_to_plot.index(ann)
        fill_color = fill[ann_idx] if isinstance(fill, list) else None
        contour_color = contour[ann_idx]

        def add_polygon(ax, polygon):
            path = polygon_to_path(polygon)

            patch = PathPatch(
                path,
                facecolor=fill_color if fill_color else 'none',
                edgecolor=contour_color,
                alpha=alpha if fill_color else 1,
                linewidth=1
            )
            ax.add_patch(patch)

        if geom.geom_type == 'Polygon':
            add_polygon(ax, geom)

        elif geom.geom_type == 'MultiPolygon':
            for polygon in geom:
                add_polygon(ax, polygon)

    ax.set_title('Annotations', fontsize=14)
    ax.set_ylim(ax.get_ylim()[::-1])
    ax.set_aspect('equal', adjustable='box')

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1, 1))

    bounds = histomap.data_exploded.total_bounds
    ax.set_xlim(bounds[0], bounds[2])
    ax.set_ylim(bounds[1], bounds[3])
    plt.gca().invert_yaxis()
    plt.tight_layout()

    # ---- SAVE FIGURE ----
    if save is not None:
        if not isinstance(save, str):
            raise TypeError("save must be a string filepath or None.")
        _, ext = os.path.splitext(save)
        if ext.lower() not in [".png", ".pdf", ".svg", ".jpg", ".jpeg", ".tiff"]:
            raise ValueError(
                "Unsupported file extension. Use .png, .pdf, .svg, .jpg, .jpeg, or .tiff"
            )
        fig.savefig(save, bbox_inches='tight', dpi=300)
        print(f"Figure saved to: {save}")

    plt.show()




def plot_tissue_overlap(histomap, cmap="coolwarm", figsize=(10, 10),
                        display_image=True, save=None):
    """
    Plot the tissue detection overlap using spot_geodata geometries.

    Parameters:
    - histomap (HistoMap): HistoMap object with the geospatial data containing the 'tissue_detection' column.
    - cmap (str): Colormap for visualization.
    - figsize (tuple): Figure size.
    - display_image (bool): Whether to display the underlying histology image.
    - save (str or None): File path to save the figure. Format inferred from extension
                          (e.g., .png, .pdf). If None, the figure is not saved.

    Returns:
    - Matplotlib plot showing the tissue detection values.
    """

    import os

    if "tissue_detection" not in histomap.spot_geodata.columns:
        raise ValueError(
            "Column 'tissue_detection' not found in spot_geodata. "
            "Run compute_tissue_overlap first."
        )

    # Plot
    fig, ax = plt.subplots(figsize=figsize)

    # Display the image beneath the annotations
    if display_image and hasattr(histomap, 'plotting_image') and histomap.plotting_image is not None:
        extent = [0, histomap.full_res_width, histomap.full_res_height, 0]
        ax.imshow(
            histomap.plotting_image.values.transpose(1, 2, 0),
            extent=extent,
            origin='upper',
            cmap='gray'
        )

    histomap.spot_geodata.plot(
        column="tissue_detection",
        cmap=cmap,
        linewidth=0.1,
        edgecolor="black",
        legend=True,
        ax=ax
    )

    # Formatting
    ax.set_title("Tissue Detection Overlap", fontsize=14)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)
    plt.tight_layout()

    # ---- SAVE FIGURE ----
    if save is not None:
        if not isinstance(save, str):
            raise TypeError("save must be a string filepath or None.")

        _, ext = os.path.splitext(save)
        if ext.lower() not in [".png", ".pdf", ".svg", ".jpg", ".jpeg", ".tiff"]:
            raise ValueError(
                "Unsupported file extension. Use .png, .pdf, .svg, .jpg, .jpeg, or .tiff"
            )

        fig.savefig(save, bbox_inches="tight", dpi=300)
        print(f"Figure saved to: {save}")

    plt.show()

def violin_tissue_overlap(histomap, figsize=(8, 6)):
    """
    Plot a violin + boxplot of the tissue detection values, with individual spot points.

    Parameters:
    - histomap (HistoMap): histomap with the geospatial data containing the 'tissue_detection' column.
    - figsize (tuple): Figure size.

    Returns:
    - Violin + boxplot + strip plot visualization of tissue detection values.
    """

    if "tissue_detection" not in histomap.spot_geodata.columns:
        raise ValueError("Column 'tissue_detection' not found in spot_geodata. Run compute_tissue_overlap first.")

    # Set seaborn style
    sns.set_style("whitegrid")

    # Create figure and axis
    fig, ax = plt.subplots(figsize=figsize)

    # Violin plot with refined colors and transparency
    sns.violinplot(
        y=histomap.spot_geodata["tissue_detection"], 
        inner=None, 
        color="lightblue", 
        linewidth=0.8, 
        alpha=0.7, 
        ax=ax
    )

    # Boxplot overlay with stronger visibility
    sns.boxplot(
        y=histomap.spot_geodata["tissue_detection"], 
        width=0.15, 
        boxprops={"facecolor": "white", "edgecolor": "black", "linewidth": 1.2}, 
        medianprops={"color": "black", "linewidth": 1.5},
        whiskerprops={"color": "black", "linewidth": 1.2},
        capprops={"color": "black", "linewidth": 1.2},
        flierprops={"marker": "o", "markerfacecolor": "red", "markeredgecolor": "black", "markersize": 4},
        ax=ax
    )

    # Strip plot to add individual data points
    sns.stripplot(
        y=histomap.spot_geodata["tissue_detection"], 
        color="black", 
        size=3, 
        alpha=0.5, 
        jitter=True, 
        ax=ax
    )

    # Formatting
    ax.set_title("Distribution of Tissue Detection Overlap", fontsize=14, fontweight="bold")
    ax.set_ylabel("Tissue Detection (%)", fontsize=12)
    ax.set_xlabel("")  # No x-label needed for a single variable
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.show()



def plot_positive_spots(histomap, annotation, display_image=True, save=None):
    """Plot spots colored by whether they are positive for the given annotation(s).
    
    Parameters:
    - histomap: HistoMap object containing the spot geodata with annotation_positive columns.
    - annotation: A string or list of strings representing the annotation(s) to plot.
    - display_image: Boolean, whether to display the image beneath the cells (True or False).
    - save: str or None. File path to save the figure. Format inferred from extension
            (e.g., .png, .pdf). If None, the figure is not saved.
    """

    import os

    # Ensure annotation is always a list
    if isinstance(annotation, str):
        annotation = [annotation]

    # Ensure the required columns exist
    missing_cols = []
    for ann in annotation:
        positive_col = ann + "_positive"
        if positive_col not in histomap.spot_geodata.columns:
            missing_cols.append(positive_col)
    
    if missing_cols:
        raise ValueError(
            f"Error: Missing columns {', '.join(missing_cols)}. "
            "Compute positive spots first using set_positive()."
        )

    # Create figure and axis
    if len(annotation) > 1:
        fig, axes = plt.subplots(1, len(annotation), figsize=(5 * len(annotation), 10))
    else:
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        axes = [ax]

    # Get the custom color mapping if available
    color_dict = {}
    if hasattr(histomap, 'annotation_colors'):
        color_dict = dict(zip(histomap.annotation_colors['annotation'],
                              histomap.annotation_colors['color']))

    for ax, ann in zip(axes, annotation):
        gdf = histomap.spot_geodata.copy()
        positive_col = ann + "_positive"

        # Display the image beneath the spots
        if display_image and hasattr(histomap, 'plotting_image') and histomap.plotting_image is not None:
            extent = [0, histomap.full_res_width, histomap.full_res_height, 0]
            ax.imshow(
                histomap.plotting_image.values.transpose(1, 2, 0),
                extent=extent,
                origin='upper',
                cmap='gray'
            )

        negative_color = 'lightgrey'
        positive_color = color_dict.get(ann, 'red')

        positive_spots = gdf[gdf[positive_col] == True].copy()
        negative_spots = gdf[gdf[positive_col] == False].copy()

        negative_spots.plot(ax=ax, color=negative_color,
                            label='Negative', alpha=0.6)
        positive_spots.plot(ax=ax, color=positive_color,
                            label=f'{ann} Positive', alpha=0.8)

        ax.legend(title=f'{ann}')
        ax.set_title(f'Spot Plot for {ann} (Positive/Negative)', fontsize=15)
        ax.set_xlabel('X Coordinate', fontsize=12)
        ax.set_ylabel('Y Coordinate', fontsize=12)

    plt.tight_layout()

    # ---- SAVE FIGURE ----
    if save is not None:
        if not isinstance(save, str):
            raise TypeError("save must be a string filepath or None.")

        _, ext = os.path.splitext(save)
        if ext.lower() not in [".png", ".pdf", ".svg", ".jpg", ".jpeg", ".tiff"]:
            raise ValueError(
                "Unsupported file extension. Use .png, .pdf, .svg, .jpg, .jpeg, or .tiff"
            )

        fig.savefig(save, bbox_inches="tight", dpi=300)
        print(f"Figure saved to: {save}")

    plt.show()


def plot_annotation_order(histomap, fill=False, contour=None, annotation=None, display_image=False, figsize=(8, 8)):
    """Plots the annotations in 3D based on the DataFrame, respecting the plot order.
    
    Parameters:
    - histomap: HistoMap object containing annotation data.
    - fill: False, True, or a list of colors. 
            False means no fill (only contours), 
            True uses the default colormap for fill, 
            a list of colors specifies the fill color for each annotation.
            If None, uses colors from histomap.annotation_colors.
    - contour: a color or list of colors for the contours. 
               If None, uses colors from histomap.annotation_colors or default colormap.
    - annotation: a specific annotation (string) or a list of annotations to plot. If None, all annotations are plotted.
    - display_image: bool, whether to display `histomap.plotting_image` beneath the annotations.
    """
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')

    # Display the image beneath the annotations
    if display_image and hasattr(histomap, 'plotting_image') and histomap.plotting_image is not None:
        extent = [0, histomap.full_res_width, histomap.full_res_height, 0]
        ax.imshow(histomap.plotting_image.values.transpose(1, 2, 0), extent=extent, origin='upper', cmap='gray')

    cmap = plt.cm.get_cmap('tab20')
    unique_annotations = list(histomap.data_exploded['Annotation'].unique())

    # Validate annotation input
    if annotation is not None:
        if isinstance(annotation, str):
            annotation = [annotation]
        elif not isinstance(annotation, list) or not all(isinstance(a, str) for a in annotation):
            raise TypeError("Annotation must be a string, a list of strings, or None.")

        # Check for missing annotations
        missing_annotations = [ann for ann in annotation if ann not in unique_annotations]
        if missing_annotations:
            raise ValueError(
                f"Annotations {missing_annotations} not found. Available annotations are: {sorted(unique_annotations)}"
            )
        annotations_to_plot = annotation
    else:
        annotations_to_plot = unique_annotations

    annotations_to_plot = [ann for ann in annotations_to_plot if ann not in histomap.disabled_annotations]

    if histomap.disabled_annotations:
        print(f"Skipping disabled annotations: {', '.join(histomap.disabled_annotations)}")

    if len(histomap.activated_annotations) == 0:
        raise ValueError("No activated annotations to compute overlap.")

    ####
    xmax, ymax = histomap.full_res_width, histomap.full_res_height
    # Add a rectable to the size of the image on each annotation 
    rectangle = Polygon([(0, 0), (xmax, 0), (xmax, ymax), (0, ymax)])
    tmp_data_exploded = histomap.data_exploded
    
    # subset only activated annotations
    tmp_data_exploded = tmp_data_exploded[tmp_data_exploded['Annotation'].isin(histomap.activated_annotations)]
    
    # Get unique annotations
    unique_annotations = tmp_data_exploded['Annotation'].unique()
    
    # Create a dictionary mapping annotations to their plot_order
    annotation_plot_order = tmp_data_exploded.groupby("Annotation")["plot_order"].first().to_dict()
    
    # Create a new GeoDataFrame with the rectangle for each unique annotation
    rectangle_gdf = gpd.GeoDataFrame({
        'geometry': [rectangle] * len(unique_annotations),  # List of rectangles
        'Annotation': unique_annotations,  # Corresponding annotations
        'plot_order': [annotation_plot_order[ann] for ann in unique_annotations]  # Matching plot_order
    })
    
    # Concatenate the new GeoDataFrame (rectangle_gdf) with the existing one (tmp_data_exploded)
    gdf = pd.concat([tmp_data_exploded, rectangle_gdf], ignore_index=True)
    
    # Process fill colors
    if fill is None and hasattr(histomap, 'annotation_colors'):
        # Use the annotation_colors dataframe
        fill_dict = dict(zip(histomap.annotation_colors['annotation'], 
                            histomap.annotation_colors['color']))
        fill = [fill_dict.get(ann, cmap(i % 20)) for i, ann in enumerate(annotations_to_plot)]
    elif fill is True:
        fill = [cmap(idx % 20) for idx in range(len(annotations_to_plot))]
    elif isinstance(fill, list):
        if len(fill) != len(annotations_to_plot):
            raise ValueError("The length of the 'fill' list must match the number of unique annotations.")
    elif fill is False:
        fill = [None] * len(annotations_to_plot)
    
    # Process contour colors
    if contour is None and hasattr(histomap, 'annotation_colors'):
        # Use the annotation_colors dataframe
        contour_dict = dict(zip(histomap.annotation_colors['annotation'], 
                               histomap.annotation_colors['color']))
        contour = [contour_dict.get(ann, cmap(i % 20)) for i, ann in enumerate(annotations_to_plot)]
    elif contour is None:
        contour = [cmap(idx % 20) for idx in range(len(annotations_to_plot))]
    elif isinstance(contour, list):
        if len(contour) != len(annotations_to_plot):
            raise ValueError("The length of the 'contour' list must match the number of unique annotations.")
    elif isinstance(contour, str):
        contour = [contour] * len(annotations_to_plot)

    # Plot annotations in 3D
    for idx, row in gdf.iterrows():
        geom = row['geometry']
        ann = row['Annotation']
        plot_order = row['plot_order']

        if ann not in annotations_to_plot or geom.is_empty or not geom.is_valid:
            continue
        
        ann_idx = annotations_to_plot.index(ann)
        fill_color = fill[ann_idx] if isinstance(fill, list) else None
        contour_color = contour[ann_idx]

        if geom.geom_type == 'Polygon':
            x, y = geom.exterior.xy
            z = np.full_like(x, plot_order)  # Create a constant z based on plot_order
            # Swap x and z
            if fill_color:
                ax.fill(z, x, y, color=fill_color, edgecolor=contour_color, alpha=0.3, label=ann)
            else:
                ax.plot(z, x, y, color=contour_color, label=ann)
        elif geom.geom_type == 'MultiPolygon':
            for polygon in geom:
                x, y = polygon.exterior.xy
                z = np.full_like(x, plot_order)  # Create a constant z based on plot_order
                # Swap x and z
                if fill_color:
                    ax.fill(z, x, y, color=fill_color, edgecolor=contour_color, alpha=0.3, label=ann)
                else:
                    ax.plot(z, x, y, color=contour_color, label=ann)

    ax.set_xlabel('Plot Order')
    ax.set_ylabel('')
    ax.set_zlabel('')

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1, 1))

    # Hide grid lines
    ax.grid(False)
    ax.set_zticks([])
    ax.set_yticks([])
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # Set the ticks on the z-axis to only show the actual plot_order values
    ax.set_xticks(sorted(set(gdf['plot_order'])))
    plt.tight_layout()
    plt.gca().invert_xaxis()
    plt.gca().invert_zaxis()
    plt.show()



def plot_annotation_map(histomap, display_image=True, save=None, resolution="lowres", xcoords=None,ycoords=None,figsize=(10, 10) ):
    """Plot spots colored by their assigned annotations in the Annotation_map.
    Uses colors from histomap.annotation_colors if available.
    
    Parameters:
    - histomap: HistoMap object containing the spot geodata with 'Annotation_map' column.
    - display_image: Boolean, whether to display the image beneath the cells (True or False).
    - save: str or None. File path to save the figure (format inferred from extension).
    """
    
    import os

    # Ensure the required column exists
    if 'Annotation_map' not in histomap.spot_geodata.columns:
        raise ValueError(
            "Error: 'Annotation_map' column not found. "
            "Please generate the annotation map first."
        )
    if histomap.visium_type == 'visium_hd':    
        plot_visium_hd(histomap, resolution=resolution, xcoords=xcoords,ycoords=ycoords,figsize=figsize, display_image=display_image, annotation_key="Annotation_map")

    elif histomap.visium_type == 'visium':

        sdata = histomap.visium_spatialdata
        sdata = histomap.visium_spatialdata
        table_key = list(sdata.tables.keys())[0]

        if "Annotation_map" not in sdata.tables[table_key].obs:
            sdata.tables[table_key].obs["Annotation_map"] = (
                histomap.spot_geodata["Annotation_map"].values
            )
        groups, palette = get_palette_and_groups(
            histomap,         
            sdata,
            table_key,
            "Annotation_map"
        )
        plot_visium(
            sdata,
            resolution=resolution,
            xcoords=xcoords,
            ycoords=ycoords,
            figsize=figsize,
            display_image=display_image,
            annotation_key="Annotation_map",
            palette=palette,
            groups=groups) 
        
    else:
        # Create figure and axis
        fig, ax = plt.subplots(figsize=(10, 10))
        
        # Display the image beneath the spots
        if display_image and hasattr(histomap, 'plotting_image') and histomap.plotting_image is not None:
            extent = [0, histomap.full_res_width, histomap.full_res_height, 0]
            ax.imshow(
                histomap.plotting_image.values.transpose(1, 2, 0),
                extent=extent,
                origin='upper',
                cmap='gray'
            )

        gdf = histomap.spot_geodata.copy()
        
        unique_annotations = gdf['Annotation_map'].dropna().unique()
        if len(unique_annotations) == 0:
            raise ValueError("No annotations found in the 'Annotation_map' column.")
        
        # Get colors from annotation_colors if available
        color_dict = {}
        if hasattr(histomap, 'annotation_colors'):
            color_dict = dict(zip(histomap.annotation_colors['annotation'],
                                histomap.annotation_colors['color']))
            
        # Default fallback colormap
        cmap = plt.cm.get_cmap('tab20', len(unique_annotations))
        
        for i, annotation in enumerate(unique_annotations):
            annotation_spots = gdf[gdf['Annotation_map'] == annotation]
            
            if annotation_spots.empty:
                continue

            # Handle geometry type
            if annotation_spots.geometry.geom_type.iloc[0] == 'Point':
                x_coords = annotation_spots.geometry.x
                y_coords = annotation_spots.geometry.y
            else:
                x_coords = annotation_spots.geometry.centroid.x
                y_coords = annotation_spots.geometry.centroid.y

            spot_color = color_dict.get(annotation, cmap(i))
            
            ax.scatter(
                x_coords,
                y_coords,
                label=annotation,
                color=spot_color,
                s=20,
                edgecolors='k',
                alpha=0.6
            )
        
        ax.set_title('Spots Colored by Annotations in Annotation_map', fontsize=15)
        ax.set_xlabel('X Coordinate', fontsize=12)
        ax.set_ylabel('Y Coordinate', fontsize=12)
        
        ax.legend(title='Annotations')
        plt.tight_layout()

        # ---- SAVE FIGURE ----
        if save is not None:
            if not isinstance(save, str):
                raise TypeError("save must be a string filepath or None.")

            _, ext = os.path.splitext(save)
            if ext.lower() not in [".png", ".pdf", ".svg", ".jpg", ".jpeg", ".tiff"]:
                raise ValueError(
                    "Unsupported file extension. Use .png, .pdf, .svg, .jpg, .jpeg, or .tiff"
                )

            fig.savefig(save, bbox_inches="tight", dpi=300)
            print(f"Figure saved to: {save}")

        plt.show()


def plot_annotation_map_proportions(histomap, save=None):
    """
    Plots the proportion of each annotation present in
    histomap.spot_geodata['Annotation_map'] as a bar plot.
    Uses colors from histomap.annotation_colors if available.

    Parameters:
    - histomap: An instance of the histomap object.
    - save: str or None. File path to save the figure
            (format inferred from extension).
    """

    import os

    # Ensure required column exists
    if 'Annotation_map' not in histomap.spot_geodata.columns:
        raise ValueError(
            "The histomap object does not contain 'Annotation_map'. "
            "Run generate_annotation_map after computing overlap."
        )

    annotations = histomap.spot_geodata['Annotation_map']

    annotation_counts = annotations.value_counts()
    annotation_proportions = annotation_counts / annotation_counts.sum()

    # Default colors
    colors = plt.cm.Paired.colors

    # Use custom annotation colors if available
    if hasattr(histomap, 'annotation_colors'):
        color_dict = dict(zip(histomap.annotation_colors['annotation'],
                              histomap.annotation_colors['color']))

        custom_colors = []
        for i, annotation in enumerate(annotation_proportions.index):
            custom_colors.append(
                color_dict.get(annotation, plt.cm.Paired(i % 10))
            )

        if custom_colors:
            colors = custom_colors

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))
    annotation_proportions.plot(kind='bar', color=colors, ax=ax)

    ax.set_xlabel('Annotations')
    ax.set_ylabel('Proportion')
    ax.set_title('Proportions of Annotations in Annotation_map')
    ax.set_xticklabels(annotation_proportions.index, rotation=45)

    # Add percentage labels
    for i, v in enumerate(annotation_proportions):
        ax.text(i, v + 0.01, f"{v:.1%}", ha='center', fontsize=9)

    plt.tight_layout()

    # ---- SAVE FIGURE ----
    if save is not None:
        if not isinstance(save, str):
            raise TypeError("save must be a string filepath or None.")

        _, ext = os.path.splitext(save)
        if ext.lower() not in [".png", ".pdf", ".svg", ".jpg", ".jpeg", ".tiff"]:
            raise ValueError(
                "Unsupported file extension. Use .png, .pdf, .svg, .jpg, .jpeg, or .tiff"
            )

        fig.savefig(save, bbox_inches="tight", dpi=300)
        print(f"Figure saved to: {save}")

    plt.show()

def plot_cell_density(histomap, display_image=True, max_cutoff=None,
                      save=None, **plot_kwargs):
    """
    Plots the cell density (n_cell) for each spot in the histomap.

    Parameters:
    - histomap: HistoMap object.
    - display_image: bool, whether to display background image.
    - max_cutoff: float or None, maximum value for colormap scaling.
    - save: str or None. File path to save the figure (format inferred from extension).
    - **plot_kwargs: additional arguments forwarded to GeoDataFrame.plot().
    """

    import os

    # Check if 'n_cell' exists
    if 'n_cell' not in histomap.spot_geodata.columns:
        raise ValueError("'n_cell' column does not exist in spot_geodata.")

    fig, ax = plt.subplots(figsize=(10, 10))
    gdf = histomap.spot_geodata.copy()
    gdf['n_cell'] = gdf['n_cell'].astype(float)

    # Background image
    if display_image and hasattr(histomap, 'plotting_image') and histomap.plotting_image is not None:
        extent = [0, histomap.full_res_width, histomap.full_res_height, 0]
        ax.imshow(
            histomap.plotting_image.values.transpose(1, 2, 0),
            extent=extent,
            origin='upper',
            cmap='gray'
        )

    # Color scale
    vmin = gdf['n_cell'].min()
    vmax = max_cutoff if max_cutoff is not None else gdf['n_cell'].max()

    # Default plot arguments
    default_plot_kwargs = dict(
        column='n_cell',
        cmap='viridis',
        legend=True,
        legend_kwds={'label': 'Number of Cells', 'orientation': "horizontal"},
        vmin=vmin,
        vmax=vmax
    )
    default_plot_kwargs.update(plot_kwargs)  # allow user overrides

    # Plot
    gdf.plot(ax=ax, **default_plot_kwargs)

    ax.set_title('Spot Plot Colored by Cell Density (n_cell)', fontsize=15)
    ax.set_xlabel('X Coordinate', fontsize=12)
    ax.set_ylabel('Y Coordinate', fontsize=12)

    plt.tight_layout()

    # ---- SAVE FIGURE ----
    if save is not None:
        if not isinstance(save, str):
            raise TypeError("save must be a string filepath or None.")

        _, ext = os.path.splitext(save)
        if ext.lower() not in [".png", ".pdf", ".svg", ".jpg", ".jpeg", ".tiff"]:
            raise ValueError(
                "Unsupported file extension. Use .png, .pdf, .svg, .jpg, .jpeg, or .tiff"
            )

        fig.savefig(save, bbox_inches="tight", dpi=300)
        print(f"Figure saved to: {save}")

    plt.show()



def plot_visium_hd(
    histomap,
    resolution="lowres",          # "lowres" | "hires" | "global"
    xcoords=None,
    ycoords=None,
    annotation_key="Annotation_map",
    figsize=(10, 10),
    display_image=True,           # True → only shapes, no background image
    save=None,                    # None or path to save figure
    cmap=None
):
    sdata = histomap.visium_spatialdata
    bin_size = histomap.bin_size
    bin_str = f"{int(bin_size):03d}um"

    # ---- 1️⃣ Map resolution ----
    coord_map = {
        "lowres": "downscaled_lowres",
        "hires": "downscaled_hires",
        "global": "global",
    }

    if resolution not in coord_map:
        raise ValueError("resolution must be 'lowres', 'hires', or 'global'")

    coord_system = coord_map[resolution]

    # ---- 2️⃣ Select image by name pattern ----
    image_key = next(
        (k for k in sdata.images.keys() if resolution in k.lower()),
        None,
    )

    # ---- 3️⃣ Select shape via bin size ----
    shape_key = next(
        (k for k in sdata.shapes.keys() if bin_str in k),
        None,
    )
    if shape_key is None:
        raise ValueError(f"No shape with bin size {bin_str} found.")

    # ---- 4️⃣ Match table ----
    table_key = next(k for k in sdata.tables.keys() if bin_str in k)

    # ---- 5️⃣ Attach annotation if missing ----
    if annotation_key not in sdata.tables[table_key].obs:
        sdata.tables[table_key].obs[annotation_key] = (
            histomap.spot_geodata[annotation_key].values
        )

    # ---- 6️⃣ Optional crop ----
    if xcoords is not None and ycoords is not None:
        sdata = bounding_box_query(
            sdata,
            min_coordinate=[xcoords[0], ycoords[0]],
            max_coordinate=[xcoords[1], ycoords[1]],
            axes=("x", "y"),
            target_coordinate_system=coord_system,
        )

    # ---- 7️⃣ Prepare rendering ----
    if display_image:
        p = sdata.pl.render_images(image_key).pl.render_shapes(
            shape_key,
            color=annotation_key,
            datashader_reduction="mean",
            cmap=cmap,
        )
    else:
        p = sdata.pl.render_shapes(
            shape_key,
            color=annotation_key,
            datashader_reduction="mean",
            cmap=cmap,
        )

    # ---- 8️⃣ Show ----
    ax = p.pl.show(
        coordinate_systems=coord_system,
        title=f"bin_size={bin_str}",
        figsize=figsize,
        return_ax=True,
    )

    # ---- 9️⃣ Save ----
    if save is not None:
        ax.figure.savefig(save, bbox_inches="tight")
        plt.close(ax.figure)

    return p



def plot_visium(
    sdata,
    resolution="lowres",          # "lowres" | "hires" | "global"
    xcoords=None,
    ycoords=None,
    annotation_key="Annotation_map",
    figsize=(10, 10),
    display_image=True,
    save=None,
    cmap=None, 
    palette=None, groups=None
    ):

    # ---- 1️⃣ Map resolution ----
    coord_map = {
        "lowres": "downscaled_lowres",
        "hires": "downscaled_hires",
        "global": "global",
    }

    if resolution not in coord_map:
        raise ValueError("resolution must be 'lowres', 'hires', or 'global'")

    coord_system = coord_map[resolution]

    # ---- 2️⃣ Select image ----
    image_key = next(
        (k for k in sdata.images.keys() if resolution in k.lower()),
        None,
    )

    # ---- 3️⃣ Select shapes (spots) ----
    shape_key = list(sdata.shapes.keys())[0]

    # ---- 4️⃣ Select table ----
    table_key = list(sdata.tables.keys())[0]

    # ---- 5️⃣ Optional crop ----
    if xcoords is not None and ycoords is not None:
        sdata = bounding_box_query(
            sdata,
            min_coordinate=[xcoords[0], ycoords[0]],
            max_coordinate=[xcoords[1], ycoords[1]],
            axes=("x", "y"),
            target_coordinate_system=coord_system,
        )

    # ---- 6️⃣ Prepare figure ----
    fig, ax = plt.subplots(figsize=figsize)

    if display_image:
        p = sdata.pl.render_images(image_key).pl.render_shapes(
            shape_key,
            color=annotation_key,
            cmap=cmap,
            palette=palette,groups=groups,
        )
    else:
        p = sdata.pl.render_shapes(
            shape_key,
            color=annotation_key,
            cmap=cmap,
            palette=palette,groups=groups,
        )

    # ---- 7️⃣ Save / show ----
    if save is not None:
        p.pl.show(
            coordinate_systems=coord_system,
            
            ax=ax,
            return_ax=False,
        )
        fig.savefig(save, bbox_inches="tight")
        plt.close(fig)

    p.pl.show(
        coordinate_systems=coord_system,
        
        ax=ax,
    )

    return p


def get_palette_and_groups(histomap, sdata, table_key, annotation_key, default_color="#cccccc"):
    """
    Robust version:
    - Uses histomap colors when available
    - Falls back to default_color if missing
    """

    # Full mapping from histomap
    color_dict = dict(zip(
        histomap.annotation_colors['annotation'],
        histomap.annotation_colors['color']
    ))

    obs = sdata.tables[table_key].obs[annotation_key].astype("category")
    groups = list(obs.cat.categories)

    # robust mapping
    missing = [g for g in groups if g not in color_dict]
    if missing:
        print(f" Missing colors for: {missing} → using fallback '{default_color}'")

    palette = [color_dict.get(g, default_color) for g in groups]

    return groups, palette