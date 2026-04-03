# HistoMap
<img src="figures/logo.svg" align="right" width="150px" />
HistoMap is a Python library for analyzing and visualizing histological annotations alongside spatially resolved transcriptomics data (Visium, VisiumHD and Xenium). It provides tools for processing, analyzing, and visualizing GeoJSON-based tissue annotations with spatial transcriptomics spot data.  
It is integrated with QuPath and ImageJ annotations, and interface with scanpy, squidpy and Seurat through either SpatialData or generation of MetaData.  
  
  
**Documentation is available here :** https://histomaptx.readthedocs.io/en/latest/

## Installation

```bash
pip install histomaptx
```


## Core Functionality

### Loading and Processing Annotations

HistoMap automatically processes annotation files in geojson format from QuPath and extracts key information:

```python
import histomaptx as hm
# Initialize with a GeoJSON file
histo = hm.HistoMap("annotations.geojson", visium_spatial_data, '/path/to/image.tiff')

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
histo.plot_annotation_order()
```

### Controlling Annotation Order

The order in which annotations are plotted can be crucial for generating the final map:

```python
# Display current plot order
histo.display_plot_order()

# Change plot order (annotations listed first will be on top)
histo.change_plot_order(["Tumor", "Stroma", "Immune cells"])
```

### Compute overlap of annotation with spatial units

Compute the overlap between histological annotations and Visium spots:

```python
# Compute overlap between annotations and spots
histo.compute_overlap_annotation()

# Visualize spots colored by their overlap with a specific annotation
histo.plot_annotation_overlay("Tumor")

# Find spots that overlap with two different annotations
histo.plot_combined_annotation_overlap("Tumor", "Immune cells")
```

### Generate the annotation map 

Once overlaps are computed and positivity threshold set, we can generate the final annotation map

```python
histomap.generate_annotation_map(annotate_all=True)  

hm.plot_annotation_map(histomap, resolution='lowres') 
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

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use HistoMap in your research, please cite:

```
Unpublished
```

## Contact

For questions and feedback, please open an issue on the GitHub repository