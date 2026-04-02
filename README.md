# HistoMapTx
<img src="figures/logo.png" align="right" width="150px" />
HistoMap is a Python library for analyzing and visualizing histological annotations alongside spatially resolved transcriptomics data (Visium). It provides tools for processing, analyzing, and visualizing GeoJSON-based tissue annotations with spatial transcriptomics spot data.  
It is integrated with QuPath and ImageJ annotations, and interface with scanpy, squidpy and Seurat through either SpatialData or generation of MetaData.
\n
Documentation is available here : https://histomaptx.readthedocs.io/en/latest/

## Features

- **Flexible File Support**: Read annotations from various formats (GeoJSON, gzipped, or zipped files)
- **Comprehensive Metrics**: Calculate area, perimeter, circularity, solidity, and other geometric properties
- **Interactive Visualization**: Generate 2D and 3D visualizations of tissue annotations
- **Annotation Ordering**: Control the rendering order of annotations for clearer visualization
- **Spatial Analysis**: Compute overlaps between annotations and Visium spots
- **Summary Statistics**: Generate detailed morphological summaries for each annotation

## Installation

```bash
pip install histomaptx
```

## Dependencies

HistoMap requires:
- geopandas
- pandas
- numpy
- matplotlib
- plotly
- scipy
- shapely

## Quick Start

```python
import histomap as hm
import squidpy as sq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load Visium data (with squidpy for example)
adata = sq.datasets.visium_fluo_moran_test()
spatial_data = adata.uns['spatial']

# Load annotations from a GeoJSON file
histo = hm.HistoMap("annotations.geojson", spatial_data)

# Display a summary of the annotations
print(histo.generate_summary())

# Plot the annotations
histo.plot_annotations()

# Change the plotting order
histo.change_plot_order(["Tumor", "Stroma", "Immune cells"])

# Compute overlaps between annotations and Visium spots
histo.compute_overlap_annotation()

# Plot spots colored by their overlap with a specific annotation
histo.plot_annotation_overlay("Tumor")
```

## Core Functionality

### Loading and Processing Annotations

HistoMap automatically processes annotation files in various formats and extracts key information:

```python
# Initialize with a GeoJSON file
histo = HistoMap("annotations.geojson", visium_spatial_data, '/path/to/image.tiff')

```

### Generating Statistical Summaries

Get comprehensive morphological statistics about each annotation:

```python
# Generate a summary DataFrame with key metrics
summary = histo.generate_summary()
print(summary)
```

The summary includes metrics such as:
- Total area and perimeter
- Mean aspect ratio, circularity, and compactness
- Centroid coordinates
- Solidity and extent
- Polygon counts

### Visualization

HistoMap offers multiple visualization options:

#### Basic Annotation Plot

```python
# Plot annotations with default settings
histo.plot_annotations()

# Customize fill and contour colors
histo.plot_annotations(fill=True, contour="black")

# Specify custom colors for each annotation
histo.plot_annotations(fill=["red", "blue", "green"], contour=["black", "black", "black"])
```

#### 3D Visualization by Annotation Order

```python
# Create a 3D plot with annotations at different z-levels
histo.plot_annotation_order(fill=True, elevation_factor=1.5)

# Create an interactive 3D plot using Plotly
histo.plot_annotation_order_interactive(fill=True, elevation_factor=1.5)
```

### Controlling Annotation Order

The order in which annotations are plotted can be crucial for visualization:

```python
# Display current plot order
histo.display_plot_order()

# Change plot order (annotations listed first will be on top)
histo.change_plot_order(["Tumor", "Stroma", "Immune cells"])
```

### Spatial Analysis with Visium Spots

Analyze the overlap between histological annotations and Visium spots:

```python
# Compute overlap between annotations and spots
histo.compute_overlap_annotation()

# Visualize spots colored by their overlap with a specific annotation
histo.plot_annotation_overlay("Tumor")

# Find spots that overlap with two different annotations
histo.plot_combined_annotation_overlap("Tumor", "Immune cells")
```

## Advanced Usage

### Custom Annotation Colors

You can customize the colors used for annotations to make your visualizations match your publication style:

```python
# Define custom colors
annotation_colors = {
    "Tumor": "#E41A1C",
    "Stroma": "#377EB8",
    "Immune cells": "#4DAF4A"
}

# Use a color list matching the annotation order
annotations = histo.data_exploded['Annotation'].unique()
color_list = [annotation_colors[ann] for ann in annotations]

# Plot with custom colors
histo.plot_annotations(fill=color_list, contour="black")
```

### Exporting Results

```python
# Generate and save a summary to CSV
summary = histo.generate_summary()
summary.to_csv("annotation_summary.csv", index=False)

# Save the figure
fig, ax = plt.subplots(figsize=(10, 10))
histo.plot_annotations()
plt.savefig("annotations.png", dpi=300, bbox_inches="tight")
```

### Working with Spot-Level Data

After computing overlaps, you can extract spots that overlap with specific annotations:

```python
# Compute overlaps
histo.compute_overlap_annotation()

# Get spots that overlap with both tumor and immune cells
overlaps = histo.plot_combined_annotation_overlap("Tumor", "Immune cells")
overlapping_spots = overlaps[2]  # Third return value contains the overlapping spots

# Get spot IDs for further analysis
overlapping_spot_ids = overlapping_spots.index.tolist()
```

## API Reference

### HistoMap Class

The main class for working with histological annotations and spatial data.

#### Methods

- `__init__(file_name, visium_spatialdata)`: Initialize with a GeoJSON file and Visium spatial data
- `generate_spot_geodata()`: Generate circular geometries for Visium spots
- `read_geojson_based_on_type(file_name)`: Read GeoJSON data based on file format
- `_extract_annotations()`: Process annotation data from GeoJSON
- `add_area_column()`: Add area calculations to annotations
- `generate_summary()`: Create a comprehensive summary of annotation metrics
- `display_plot_order()`: Show the current plot order of annotations
- `change_plot_order(order_list)`: Modify the order in which annotations are plotted
- `plot_annotations(fill, contour)`: Create a 2D plot of annotations
- `compute_overlap_annotation()`: Calculate overlap between spots and annotations
- `plot_annotation_overlay(annotation)`: Visualize overlap between spots and a specific annotation
- `plot_annotation_order(fill, contour, elevation_factor)`: Create a 3D plot with annotations at different z-levels
- `plot_annotation_order_interactive(fill, contour, elevation_factor)`: Create an interactive 3D plot using Plotly
- `plot_combined_annotation_overlap(annotation1, annotation2)`: Find spots overlapping with two annotations

## FAQ

### How does HistoMap handle large annotation files?

HistoMap efficiently processes GeoJSON files using GeoPandas' spatial indexing capabilities, which helps manage large datasets. For very large files, you may need to increase your system's memory allocation.

### Can I use HistoMap with non-Visium spatial data?

While HistoMap is optimized for Visium data, you can adapt it for other spatial transcriptomics platforms by constructing a compatible spatial data object.

### How can I integrate HistoMap with other spatial analysis tools?

HistoMap works well with the broader spatial transcriptomics ecosystem, including:
- Squidpy for additional spatial statistics
- Scanpy for cell-type identification
- Seaborn for advanced visualizations of results


## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use HistoMap in your research, please cite:

```
Unpublished
```

## Contact

For questions and feedback, please open an issue on the GitHub repository