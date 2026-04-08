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
import geopandas as gpd
from shapely.geometry import Polygon, Point
import numpy as np
from . import histomap_utils
import rasterio
import tifffile
from rasterio.errors import RasterioIOError

class HistoMap:
    def __init__(self, file_name, visium_spatialdata, full_res_path, visium_type=None, bin_size=None):
        self.file_name = file_name
        self.visium_spatialdata = visium_spatialdata

        self.visium_type = self.detect_spatialdata_type(visium_type)
        self.data = self.read_geojson_based_on_type(file_name)
        self.data_exploded = self.data.explode().set_crs(None, allow_override=True)

        self._extract_annotations()
        self.add_area_column()
        self.bin_size = bin_size
        self.spot_geodata = self.generate_spot_geodata()


        self.overlay_computed = False
        self.segmentation_dataframe = False
        self.plotting_image = self._retrieve_hires_image_from_spatialdata()
        self.full_res_path = full_res_path
        self.full_res_height, self.full_res_width = self._get_full_res_dimensions()
        self.tissue_detection = False
        self.activated_annotations = self.data_exploded['Annotation'].unique()
        self.disabled_annotations = []
        self.positive_threshold = self._setup_positive_threshold()
        self._setup_plot_order()
        self.annotation_colors = self._setup_annotation_color()
        
            
    def __str__(self):
        try:
            unique_annotations = self.data_exploded['Annotation'].unique()
            num_annotations = len(unique_annotations)
        except KeyError:
            unique_annotations = []
            num_annotations = 0

        statement = (
            f"HistoMap Object\n\n"
            f"--------Annotation--------\n"
            f"📂 Annotation file used: {self.file_name}\n"
            f"📂 Full resolution image used: {self.full_res_path}\n"
            f"🔢 Number of annotations present: {num_annotations}\n"
            f"🏷 Annotations: {', '.join(map(str, unique_annotations)) if num_annotations > 0 else 'N/A'}"
            f"\n🏷 Activated annotations: {', '.join(map(str, self.activated_annotations)) if num_annotations > 0 else 'N/A'}"
            f"\n🏷 Disabled annotations: {', '.join(map(str, self.disabled_annotations)) if num_annotations > 0 else 'None'}"
            f"\n\n-------Spatial Data--------\n"
            f"Type is {self.visium_type}\n\n"
            f"{self.visium_spatialdata}\n\n"
            f"-----HistoMap Workflow------\n"

        )
        
        return statement
    
    def _setup_annotation_color(self):
        """
        Creates a dataframe with columns 'annotation' and 'color' at object initialization.
        Maps each unique annotation to a distinct color using a color cycle.
        
        Returns:
        pandas.DataFrame: A dataframe containing annotation names and their assigned colors
        """
        import matplotlib.colors as mcolors
        
        # Get unique annotations
        unique_annotations = self.data_exploded['Annotation'].unique()
        
        # Create a color cycle using default matplotlib color cycle
        default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        
        # If we have more annotations than default colors, extend with more colors
        if len(unique_annotations) > len(default_colors):
            # Add more colors from matplotlib's tab20 colormap
            additional_colors = list(mcolors.TABLEAU_COLORS.values())
            available_colors = default_colors + additional_colors
        else:
            available_colors = default_colors
        
        # Create the dataframe
        color_data = []
        for i, annotation in enumerate(unique_annotations):
            color_idx = i % len(available_colors)  # Cycle through colors if needed
            color_data.append({
                'annotation': annotation,
                'color': available_colors[color_idx]
            })
        
        annotation_colors = pd.DataFrame(color_data)
        
        # Store as instance attribute
        self.annotation_colors = annotation_colors
        
        return annotation_colors
    
    def change_annotation_color(self, annotation_color_dict):
        """
        Updates the color for specified annotations in the class.
        
        Parameters:
        annotation_color_dict (dict): Dictionary mapping annotation names to colors
                                    e.g. {"Tumor": "red", "Stroma": "blue"}
        
        Returns:
        bool: True if all updates were successful, False otherwise
        """
        import matplotlib.colors as mcolors
        
        all_updates_successful = True
        
        for annotation, color in annotation_color_dict.items():
            # Check if annotation exists in the data
            if annotation not in self.data_exploded['Annotation'].unique():
                print(f"Warning: Annotation '{annotation}' does not exist in the data.")
                all_updates_successful = False
                continue
            
            # Validate color
            try:
                # This will raise ValueError if color is not recognized
                mcolors.to_rgba(color)
                
                # Update color in the annotation dataframe
                mask = self.annotation_colors['annotation'] == annotation
                if any(mask):
                    self.annotation_colors.loc[mask, 'color'] = color
                    print(f"Updated color for '{annotation}' to '{color}'")
                else:
                    # If annotation exists in data but not in color mapping (unlikely but possible)
                    new_row = pd.DataFrame({'annotation': [annotation], 'color': [color]})
                    self.annotation_colors = pd.concat([self.annotation_colors, new_row], ignore_index=True)
                    print(f"Added color mapping for '{annotation}': '{color}'")
            except ValueError:
                print(f"Warning: Color '{color}' is not recognized by matplotlib.")
                all_updates_successful = False
        
        return all_updates_successful

    def display_annotation_color(self):
        """
        Prints the dataframe containing annotations and their corresponding colors.
        
        Returns:
        pandas.DataFrame: A dataframe with annotation and color columns
        """
        if hasattr(self, 'annotation_colors') and not self.annotation_colors.empty:
            print("Annotation colors:")
            display_df = self.annotation_colors.copy()
            
            # Add a status column to show if annotation is active or disabled
            display_df['status'] = display_df['annotation'].apply(
                lambda x: 'Active' if x in self.activated_annotations else 'Disabled'
            )
            
            return display_df
        else:
            print("No annotation colors have been set up.")
            return None

    def _setup_positive_threshold(self):
        """Get the annotation and set the positive threshold to 50% by default"""
        return  pd.DataFrame({"annotation": self.activated_annotations, "threshold": 50, "Overlay":False})

    def detect_spatialdata_type(self, visium_type) -> str:
        """
        Detect whether a SpatialData object was generated from
        Visium or Visium HD data.

        Returns: "visium", "visium_hd", or "unknown"
        """
        allowed_types = {None, "visium", "visium_hd", "xenium"}
        if visium_type not in allowed_types:
            raise ValueError(
                f"Invalid visium_type '{visium_type}'. Must be one of: None, 'visium', 'visium_hd'."
            )
        if visium_type == None:
            shapes = list(self.visium_spatialdata.shapes.keys())
            tables = list(self.visium_spatialdata.tables.keys())

            # Heuristic 1: multiple resolutions (HD)
            if any("_002um" in k or "_008um" in k or "_016um" in k for k in shapes + tables):
                return "visium_hd"

            # Heuristic 2: geometry type (HD = Polygon, Visium = Point)
            try:
                import shapely
                geom_types = {g.geom_type for name, df in self.visium_spatialdata.shapes.items() for g in df.geometry.head(100)}
                if "Polygon" in geom_types or "MultiPolygon" in geom_types:
                    return "visium_hd"
            except Exception:
                pass

            # Heuristic 3: number of shapes
            total_shapes = sum(len(df) for df in self.visium_spatialdata.shapes.values())
            if total_shapes > 1e5:
                return "visium_hd"
            elif total_shapes < 1e4:
                return "visium"

            # Fallback
            return "visium"
        else:
            return visium_type

    
    def change_positive_threshold(self, new_threshold):
        """Update threshold values based on the provided dictionary, with validation"""
        for annotation, threshold in new_threshold.items():
            if annotation not in self.positive_threshold['annotation'].values:
                raise ValueError(f"Annotation '{annotation}' not found in the DataFrame.")
            if not (0 <= threshold <= 100):
                raise ValueError(f"Threshold for '{annotation}' must be between 0 and 100.")

            # Check if Overlay is True for the annotation
            overlay_status = self.positive_threshold.loc[
                self.positive_threshold['annotation'] == annotation, 'Overlay'
            ].values[0]  # Extract the value as a boolean

            if not overlay_status:
                raise ValueError(
                    f"Overlay for '{annotation}' should be computed first using compute_annotation_overlap."
                )

            # Apply updates
            self.positive_threshold.loc[self.positive_threshold['annotation'] == annotation, 'threshold'] = threshold
            self.set_positive(annotation, threshold)

    def set_positive(self, annotation, threshold):
        """Set the positive column for the given annotation and threshold."""
        if annotation+'_positive' not in self.spot_geodata.columns:
            self.spot_geodata[annotation+'_positive'] = False  # Initialize the column if not present

        # Update 'Annotation_positive' based on the given threshold
        self.spot_geodata[annotation+'_positive'] = self.spot_geodata[annotation+"_overlap"] >= threshold

        print(f"Updated 'Annotation_positive' for '{annotation}' with threshold {threshold}%.")

    def display_positive_threshold(self):
        """Return the threshold DataFrame"""
        return self.positive_threshold.copy()  # Return a copy to prevent accidental modifications




    def _get_full_res_dimensions(self):
        """Return (height, width) of full resolution image."""

        try:
            # Try rasterio first
            with rasterio.open(self.full_res_path) as image:
                return image.height, image.width

        except RasterioIOError as e:
            # Fallback for JP2000 / OME-TIFF issues
            try:
                with tifffile.TiffFile(self.full_res_path) as tif:
                    page = tif.pages[0]
                    return page.imagelength, page.imagewidth
            except Exception as tif_error:
                raise RuntimeError(
                    f"Failed to read image dimensions using both rasterio "
                    f"and tifffile.\nRasterio error: {e}\n"
                    f"Tifffile error: {tif_error}"
                )



    def _retrieve_hires_image_from_spatialdata(self, preferred_keys=None):
        """
        Retrieve highest resolution image from SpatialData object.

        Parameters
        ----------
        preferred_keys : list[str], optional
            List of substrings to prioritize when searching for image keys.
            Example: ["hires", "morphology", "focus"]

        Returns
        -------
        xarray.DataArray or None
        """

        images = self.visium_spatialdata.images

        if not images:
            return None

        preferred_keys = preferred_keys or ["hires", "morphology", "focus"]

        #  Try preferred matches first
        for key in images.keys():
            if any(p in key.lower() for p in preferred_keys):
                img = images[key]
                return self._get_highest_resolution(img)

        # 2Otherwisereturn first image
        first_key = next(iter(images))
        return self._get_highest_resolution(images[first_key])
        
    def _get_highest_resolution(self, image):
        # Extract DataArray from DataTree if needed
        if hasattr(image, "children"):  # DataTree
            scale0 = list(image.children.values())[0]
            image = scale0["image"]  # DataArray

        # Ensure channel-first format (c, y, x)
        if "c" not in image.dims:
            raise ValueError("Image must have channel dimension 'c'")

        n_channels = image.sizes["c"]

        # --- Normalize channel count ---
        if n_channels == 1:
            return image  # already grayscale

        elif n_channels == 3:
            return image  # already RGB-compatible

        elif n_channels > 3:
            # Default behavior: use first channel (e.g. DAPI)
            # Keep dimension so shape remains (1, y, x)
            return image.isel(c=0).expand_dims("c")

        else:
            raise ValueError("Unexpected number of channels")

    
    

    def generate_spot_geodata(self):
        """Function to generate a dataframe from the spatialdata object as circle geodata"""
        # Extract spatial coordinates from the GeoDataFrame in the Shapes section
        # Take the first Shapes key (could be improved for multiple resolutions)
        image_id = list(self.visium_spatialdata.shapes.keys())[0]  
        coords_gdf = self.visium_spatialdata.shapes[image_id]  # GeoDataFrame

        if self.visium_type == "visium":
            # Standard Visium: create circular polygons from spot centers
            spot_radius = float(coords_gdf['radius'].iloc[0])
            gdf = gpd.GeoDataFrame(
                {"spot_id": range(len(coords_gdf))},
                geometry=[Point(x, y).buffer(spot_radius) for x, y in zip(coords_gdf.geometry.x, coords_gdf.geometry.y)],
                crs=coords_gdf.crs
            )

        elif self.visium_type == "visium_hd":

            import re

            shape_keys = list(self.visium_spatialdata.shapes.keys())

            # Extract bin size from key pattern like 'L5E_square_002um'
            key_bin_map = {}

            for key in shape_keys:
                match = re.search(r"_(\d+)um$", key)
                if match:
                    bin_value = int(match.group(1))  # handles leading zeros automatically
                    key_bin_map[bin_value] = key

            if not key_bin_map:
                raise ValueError(
                    "Could not detect bin sizes from spatialdata shape keys. "
                    "Expected pattern like '*_016um'"
                )

            available_bins = sorted(key_bin_map.keys())

            # ---- Selection logic ----
            if self.bin_size is not None:
                if self.bin_size not in available_bins:
                    raise ValueError(
                        f"Requested bin_size {self.bin_size} not available. "
                        f"Available bin sizes: {available_bins}"
                    )
                selected_bin = self.bin_size
            else:
                # No bin size specified
                if len(available_bins) > 1:
                    selected_bin = 16
                else:
                    selected_bin = available_bins[0]

            selected_key = key_bin_map[selected_bin]
            coords_gdf = self.visium_spatialdata.shapes[selected_key]

            gdf = coords_gdf.copy().reset_index(drop=True)
            gdf["spot_id"] = range(len(gdf))

        elif self.visium_type == "xenium":
            pixel_size = 0.2125

            gdf = coords_gdf.copy()

            # Scale microns → pixels
            gdf["geometry"] = gdf.geometry.scale(
                xfact=1/pixel_size,
                yfact=1/pixel_size,
                origin=(0, 0)
            )

            # Fix invalid geometries
            from shapely.validation import make_valid
            import shapely

            gdf["geometry"] = gdf.geometry.apply(make_valid)

            # Snap to precision grid (prevents GEOS side location conflicts)
            gdf["geometry"] = shapely.set_precision(gdf.geometry, grid_size=1e-6)

            gdf = gdf.reset_index(drop=True)
            gdf["spot_id"] = range(len(gdf))



        else:
            raise ValueError(f"Unsupported visium_type '{self.visium_type}'")

        return gdf

    def read_geojson_based_on_type(self, file_name):
        """Function to detect file type and read the file accordingly"""
        if file_name.endswith('.gz'):
            with gzip.open(file_name, 'rb') as f:
                return gpd.read_file(io.BytesIO(f.read()))
        elif zipfile.is_zipfile(file_name):
            with zipfile.ZipFile(file_name, 'r') as zip_ref:
                geojson_file = zip_ref.namelist()[0]
                with zip_ref.open(geojson_file) as f:
                    return gpd.read_file(f)
        elif file_name.endswith('.geojson'):
            return gpd.read_file(file_name)
        else:
            raise ValueError("Unsupported file format")

    def _extract_annotations(self):
        """Extract annotations from 'classification' column and compute area."""
        cols = self.data_exploded.columns

        if 'classification' in cols and self.data_exploded['classification'].notna().any():
            self.data_exploded = self.data_exploded.dropna(subset=['classification'])
            self.data_exploded['classification'] = self.data_exploded['classification'].apply(
                lambda x: x if isinstance(x, dict) else ast.literal_eval(x)
            )
            self.data_exploded['Annotation'] = self.data_exploded['classification'].apply(
                lambda x: x.get('name', None).replace(' ', '_') if x.get('name') else None
            )
        elif 'name' in cols:
            # QuPath export without classification: annotation name is in the 'name' column directly
            self.data_exploded['Annotation'] = self.data_exploded['name'].apply(
                lambda x: x.replace(' ', '_') if isinstance(x, str) else None
            )
        else:
            raise ValueError(
                "Could not find annotation names in the GeoJSON file. "
                "Make sure your QuPath annotations have a name or classification set."
            )

    def _setup_plot_order(self):
        """Set up the plot_order based on the alphabetical order of annotations."""
        # Initialize plot_order if it doesn't exist
        if 'plot_order' not in self.data_exploded.columns:
            # Sort annotations alphabetically and assign a plot_order
            sorted_annotations = sorted(self.data_exploded['Annotation'].unique())
            
            # Create a mapping from annotation to plot_order based on alphabetical order
            annotation_to_order = {annotation: idx for idx, annotation in enumerate(sorted_annotations)}
            
            # Assign the plot_order based on the sorted order
            self.data_exploded['plot_order'] = self.data_exploded['Annotation'].map(annotation_to_order)

        # Sort the data by the plot_order
        self.data_exploded = self.data_exploded.sort_values(by='plot_order', ascending=True)

    def add_segmentation(self, segmentation_file, file_type="geojson"):
        """ read geojson segmentation from qupath and return a segmentation dataset """
        segmentation_df = self.read_geojson_based_on_type(segmentation_file)
        # only keep cell and not annotation
        segmentation_df = segmentation_df[segmentation_df['objectType'] == 'cell']
        self.spot_geodata['n_cell'] = 0
        # Loop over each spot and count overlapping cells
        for idx, spot in self.spot_geodata.iterrows():
            # Get the geometry of the current spot
            spot_geom = spot['geometry']
            
            # Find all cells that intersect with the spot
            overlapping_cells = segmentation_df[segmentation_df['geometry'].intersects(spot_geom)]
            
            # Count the number of overlapping cells
            self.spot_geodata.at[idx, 'n_cell'] = len(overlapping_cells)
        self.segmentation_dataframe = segmentation_df
        print('Segmentation added. Number of cells in the file :'+ str(len(segmentation_df.index))+ ". Average cell per spot : " + str(self.spot_geodata['n_cell'].mean()))
            
    def add_area_column(self):
        """Adds an area column to the dataframe (in square units)."""
        self.data_exploded['area'] = self.data_exploded['geometry'].area

    def generate_summary(self):
        """Generates a summary DataFrame for only activated annotations."""
        summary = pd.DataFrame()

        # Filter data to keep only activated annotations
        filtered_data = self.data_exploded[self.data_exploded['Annotation'].isin(self.activated_annotations)]

        # Compute metrics as before
        if not filtered_data.empty:
            summary['total_area'] = filtered_data.groupby('Annotation')['geometry'].apply(lambda x: x.area.sum())
            summary['total_perimeter'] = filtered_data.groupby('Annotation')['geometry'].apply(lambda x: x.length.sum())
            summary['polygon_count'] = filtered_data.groupby('Annotation').size()

            total_area = filtered_data.geometry.area.sum()
            summary['area_fraction'] = summary['total_area'] / total_area if total_area > 0 else 0

            summary.reset_index(inplace=True)

        return summary

    def display_plot_order(self):
        """Displays the plot order for each annotation in the DataFrame.
        If plot_order doesn't exist, it assigns a default plot order based on alphabetical annotation."""
        
        print("Current Plot Order:")
        
        # Check if 'plot_order' exists in the dataframe
        if 'plot_order' in self.data_exploded.columns:
            # Sort by 'plot_order' if it exists
            sorted_data = self.data_exploded[['Annotation', 'plot_order']].sort_values(by='plot_order', ascending=False)
        else:
            # If 'plot_order' doesn't exist, create a default plot order based on alphabetical sorting of 'Annotation'
            # Assign a rank based on the alphabetical order of 'Annotation'
            self.data_exploded['plot_order'] = self.data_exploded['Annotation'].astype('category').cat.codes
            
            # Sort by the default 'plot_order'
            sorted_data = self.data_exploded[['Annotation', 'plot_order']].sort_values(by='plot_order', ascending=True)

        # Filter out disabled annotations from the plot order display
        active_annotations = self.data_exploded[self.data_exploded['Annotation'].isin(self.activated_annotations)]
        
        unique_annotations = active_annotations[['Annotation', 'plot_order']].drop_duplicates()
        # Sort by plot_order
        unique_annotations = unique_annotations.sort_values(by='plot_order', ascending=True)
        return unique_annotations


    def change_plot_order(self, order_list:list):
        """
        Changes the plot order in the DataFrame based on the provided list of annotations.
        The list must contain only activated annotations and must have the same length as self.activated_annotations.
        No annotations can be missing from the order_list.

        Parameters:
        - order_list: A list of annotations in the desired plot order (top first).
        """
        if not isinstance(order_list, list):
            raise TypeError('Order_list should be a list of all activated annotations')
        # Ensure all annotations in the order_list are activated
        invalid_activated_annotations = set(order_list) - set(self.activated_annotations)
        
        if invalid_activated_annotations:
            raise ValueError(f"The following annotations are not activated and cannot be part of the custom order: {', '.join(invalid_activated_annotations)}. "
                            f"Activated annotations are: {', '.join(self.activated_annotations)}.")
        
        # Check if order_list has the same length as activated_annotations
        if len(order_list) != len(self.activated_annotations):
            raise ValueError(f"order_list must have the same length as activated_annotations. "
                            f"Expected length: {len(self.activated_annotations)}, but got length: {len(order_list)}. "
                            f"Activated annotations are: {', '.join(self.activated_annotations)}.")
        
        # Ensure all annotations in the order_list exist in the DataFrame
        valid_annotations = set(self.data_exploded['Annotation'].unique())
        invalid_order_annotations = set(order_list) - valid_annotations

        if invalid_order_annotations:
            raise ValueError(f"The provided order list contains invalid annotations: {', '.join(invalid_order_annotations)}. Available annotations are: {', '.join(valid_annotations)}.")

        # Ensure no annotations are missing in the order_list
        missing_annotations = set(self.activated_annotations) - set(order_list)
        
        if missing_annotations:
            raise ValueError(f"The following activated annotations are missing from the order list: {', '.join(missing_annotations)}. "
                            f"Please include all activated annotations. Activated annotations are: {', '.join(self.activated_annotations)}.")

        # Create a dictionary mapping each annotation to its new plot order
        plot_order_dict = {annotation: idx for idx, annotation in enumerate(order_list)}

        # Map the annotations in the DataFrame to the new plot order
        self.data_exploded['plot_order'] = self.data_exploded['Annotation'].map(plot_order_dict)

        # Sort by 'plot_order' in ascending order to plot top annotations first
        self.data_exploded = self.data_exploded.sort_values('plot_order', ascending=True)

        # Display the updated plot order
        self.display_plot_order()


        
    def compute_annotation_overlap(self, annotation_names):
        """Compute overlap, considering specific annotations or all activated annotations if 'all' is passed."""
        
        # If 'all' is passed, use all activated annotations
        if annotation_names == 'all':
            annotation_names = self.activated_annotations

        # Ensure annotation_names is a list
        if isinstance(annotation_names, str):
            annotation_names = [annotation_names]

        # Warn if any annotations are disabled
        disabled_annotations = [ann for ann in annotation_names if ann not in self.activated_annotations]
        if disabled_annotations:
            print(f"Some annotations are disabled and will be ignored: {', '.join(disabled_annotations)}")

        # Remove disabled annotations
        annotation_names = [ann for ann in annotation_names if ann in self.activated_annotations]

        if not annotation_names:
            raise ValueError("No valid annotations to compute overlap. They may be disabled.")

        # Check if all annotations exist in data_exploded
        missing_annotations = [ann for ann in annotation_names if ann not in self.data_exploded['Annotation'].unique()]
        if missing_annotations:
            raise ValueError(f"Annotations not found: {', '.join(missing_annotations)}")

        # Filter the data to include only the selected annotations
        annotation_subset = self.data_exploded[self.data_exploded['Annotation'].isin(annotation_names)]
        print('Starting calculating overlap')
        # Compute overlap for the selected annotations
        self.spot_geodata = histomap_utils.calculate_annotation_overlap_fast(self.spot_geodata, annotation_subset)
        for ann in annotation_names:
            overlap_col = ann + "_overlap"
            if self.visium_type == "visium_hd":

                sdata = self.visium_spatialdata
                bin_str = f"{int(self.bin_size):03d}um"
                table_key = next(k for k in sdata.tables.keys() if bin_str in k)

                if overlap_col not in sdata.tables[table_key].obs:
                    sdata.tables[table_key].obs[overlap_col] = (
                        self.spot_geodata[overlap_col].values
                    )

            elif self.visium_type == "visium":

                sdata = self.visium_spatialdata
                table_key = list(sdata.tables.keys())[0]

                if overlap_col not in sdata.tables[table_key].obs:
                    sdata.tables[table_key].obs[overlap_col] = (
                        self.spot_geodata[overlap_col].values
                    )


        # Call set_positive for each annotation in the list
        for ann in annotation_names:
            threshold = self.positive_threshold.loc[self.positive_threshold['annotation'] == ann, 'threshold'].values[0]
            self.set_positive(ann, threshold)

        # Update Overlay status for the computed annotations
        self.positive_threshold.loc[
            self.positive_threshold['annotation'].isin(annotation_names), 'Overlay'
        ] = True

        print(f"Computed annotation overlap for: {', '.join(annotation_names)}. Updated Overlay status.")




    def compute_tissue_overlap(self, positive, negative=None):
        """
        Compute tissue overlap using positive and optional negative annotations.
        Adds a 'tissue_detection' column to spot_geodata.
        
        - `positive`: List or string of annotation names considered as tissue.
        - `negative`: List or string of annotation names to exclude (optional).
        """

        # Ensure positive is a list
        if isinstance(positive, str):
            positive = [positive]
        
        # Validate annotations
        missing_positive = [ann for ann in positive if ann not in self.data_exploded['Annotation'].unique()]
        if missing_positive:
            raise ValueError(f"Positive annotations not found: {', '.join(missing_positive)}")

        if negative:
            if isinstance(negative, str):
                negative = [negative]
            missing_negative = [ann for ann in negative if ann not in self.data_exploded['Annotation'].unique()]
            if missing_negative:
                raise ValueError(f"Negative annotations not found: {', '.join(missing_negative)}")

        # Compute overlap for positive annotations
        self.compute_annotation_overlap(positive)
        
        # Create the tissue_detection column
        self.spot_geodata["tissue_detection"] = self.spot_geodata[[f"{ann}_overlap" for ann in positive]].sum(axis=1)

        # If negative annotations are provided, compute and subtract their overlap
        if negative:
            self.compute_annotation_overlap(negative)
            self.spot_geodata["tissue_detection"] -= self.spot_geodata[[f"{ann}_overlap" for ann in negative]].sum(axis=1)

        # Ensure no negative values after subtraction
        self.spot_geodata["tissue_detection"] = self.spot_geodata["tissue_detection"].clip(lower=0)
        self.tissue_detection = True

    def tissue_detection_summary(self):
        """
        Generate a summary table with statistics for the tissue detection values.

        Returns:
        - DataFrame with count, mean, std, min, 25%, 50%, 75%, and max values.
        """
        if "tissue_detection" not in self.spot_geodata.columns:
            raise ValueError("Column 'tissue_detection' not found in spot_geodata. Run compute_tissue_overlap first.")

        # Compute summary statistics
        summary = self.spot_geodata["tissue_detection"].describe().to_frame().T

        # Rename index for clarity
        summary.index = ["Tissue Detection Statistics"]

        return summary
    
    def filter_tissue_overlap(self, threshold=0):
        """
        Filter the HistoMap object based on a tissue detection threshold.

        Parameters:
        - threshold (float): Minimum tissue detection percentage to keep a spot.

        Returns:
        - The HistoMap object with filtered spot_geodata and visium_spatialdata.
        """
        import copy
        if "tissue_detection" not in self.spot_geodata.columns:
            raise ValueError("Column 'tissue_detection' not found in spot_geodata. Run compute_tissue_overlap first.")

        # Subset spot_geodata
        filtered_spots = self.spot_geodata[self.spot_geodata["tissue_detection"] >= threshold].copy()

        histomap_filtered = copy.deepcopy(self)
        histomap_filtered.spot_geodata = filtered_spots
        # Subset visium_spatialdata based on the filtered spot indices
        histomap_utils.filter_spatialdata(histomap_filtered, filtered_spots.index)
        return histomap_filtered
    
    def disable_annotation(self, annotations):
        """
        Move the given annotation(s) to the disabled annotations list.
        If annotations are present in activated annotations, they are removed from it.
        Set their plot_order to None.

        Parameters:
        - annotations: A string or list of annotations to disable.
        """
        # Ensure annotations is a list
        if isinstance(annotations, str):
            annotations = [annotations]

        # Ensure disabled_annotations exists
        if not hasattr(self, 'disabled_annotations'):
            self.disabled_annotations = []

        # Find already disabled annotations
        already_disabled = [ann for ann in annotations if ann in self.disabled_annotations]
        if already_disabled:
            print(f"Annotations already disabled: {', '.join(already_disabled)}")

        # Find annotations that can be disabled
        to_disable = [ann for ann in annotations if ann in self.activated_annotations]

        if not to_disable:
            print("No new annotations to disable.")
            return

        # Remove from activated annotations and add to disabled_annotations
        self.activated_annotations = [ann for ann in self.activated_annotations if ann not in to_disable]
        self.disabled_annotations = list(set(self.disabled_annotations) | set(to_disable))

        # Set the plot order of the disabled annotations to None
        self.data_exploded.loc[self.data_exploded['Annotation'].isin(to_disable), 'plot_order'] = None

        # Re-adjust the plot order of the remaining activated annotations (renumber them)
        remaining_activated = self.data_exploded[self.data_exploded['Annotation'].isin(self.activated_annotations)]

        # Group by 'Annotation' and assign a unique plot_order to each group
        remaining_activated['plot_order'] = remaining_activated.groupby('Annotation').ngroup()

        # Update plot_order for matching annotations directly
        for _, row in remaining_activated.iterrows():
            self.data_exploded.loc[self.data_exploded['Annotation'] == row['Annotation'], 'plot_order'] = row['plot_order']

        # Print updated lists for feedback
        print(f"Updated activated annotations: {self.activated_annotations}")
        print(f"Updated disabled annotations: {self.disabled_annotations}")




    def activate_annotation(self, annotations):
        """
        Move the given annotation(s) from the disabled annotations list back to the activated annotations list.
        If annotations are present in disabled annotations, they are removed from it.
        Set their plot order at the end.

        Parameters:
        - annotations: A string or list of annotations to activate.
        """
        # Ensure annotations is a list
        if isinstance(annotations, str):
            annotations = [annotations]

        # Ensure activated_annotations exists
        if not hasattr(self, 'activated_annotations'):
            self.activated_annotations = []

        # Find already activated annotations
        already_activated = [ann for ann in annotations if ann in self.activated_annotations]
        if already_activated:
            print(f"Annotations already activated: {', '.join(already_activated)}")

        # Find annotations that can be activated (those in disabled_annotations)
        to_activate = [ann for ann in annotations if ann in self.disabled_annotations]

        if not to_activate:
            print("No new annotations to activate.")
            return

        # Remove from disabled_annotations and add to activated_annotations
        self.disabled_annotations = [ann for ann in self.disabled_annotations if ann not in to_activate]
        self.activated_annotations.extend(to_activate)

        # Set the plot order for the activated annotations at the end
        current_max_plot_order = self.data_exploded['plot_order'].max() if not self.data_exploded['plot_order'].isnull().all() else -1
        for ann in to_activate:
            self.data_exploded.loc[self.data_exploded['Annotation'] == ann, 'plot_order'] = current_max_plot_order + 1
            current_max_plot_order += 1

        # Print updated lists for feedback
        print(f"Updated activated annotations: {self.activated_annotations}")
        print(f"Updated disabled annotations: {self.disabled_annotations}")

    
    def to_anndata(self):
        """Return a anndata object with metadata from histomap"""
        adata = self.visium_spatialdata.tables['table'].copy()
        spot_geodata = self.spot_geodata.copy()
        # transfer the metadata 
        df1 = spot_geodata.drop(columns='geometry', inplace=False)
        df2 = adata.obs
        # Perform the merge, keeping only the columns from df1 (spot_geodata)
        df_merged = pd.merge(df2, df1, on='spot_id', how='left', suffixes=('', '_drop'))
        # Drop the unwanted columns (from df2) 
        df_merged = df_merged.loc[:, ~df_merged.columns.str.endswith('_drop')]
        df_merged.index = df2.index
        adata.obs = df_merged
        return adata
    
    def to_spatialdata(self):
        """Return a SpatialData object with metadata from HistoMap"""
        spatialdata_obj = self.visium_spatialdata
        adata = spatialdata_obj.tables['table'].copy()
        # transfer the metadata 
        df1 = self.spot_geodata.drop(columns='geometry')
        df2 = adata.obs
        df_merged = pd.merge(df2, df1, on='spot_id', how='left')
        df_merged.index = df2.index
        adata.obs = df_merged
        spatialdata_obj.tables['table'] = adata
        return spatialdata_obj
    
    def spot_metadata_to_df(self):
        """Return a csv with the metadata"""
        adata = self.visium_spatialdata.tables['table']
        # transfer the metadata 
        df1 = self.spot_geodata.drop(columns='geometry', inplace=False)
        df2 = self.visium_spatialdata.tables['table'].obs
        df_merged = pd.merge(df2, df1, on='spot_id', how='left')
        df_merged.index = df2.index
        return df_merged
    
    def generate_annotation_map(self, annotate_all=True):
        """Generate annotation map based on plot_order and positivity"""
        # Create 'Annotation_map' column if it doesn't exist or reset it
        self.spot_geodata['Annotation_map'] = "None"
        for annotation in self.activated_annotations:
            # Check if the annotation exists in 'positive_threshold' with Overlay set to True
            positive_threshold_row = self.positive_threshold[self.positive_threshold['annotation'] == annotation]
            
            if positive_threshold_row.empty:
                raise ValueError(f"'{annotation}' not found in the object.")
            
            if not positive_threshold_row['Overlay'].iloc[0]:
                raise ValueError(f"Overlap for '{annotation}' was not computed. Run histomap.compute_overlap() first or disable the annotation layer with histomap.disable_annotation().")
            
            # Gather the corresponding 'Annotation_positive' column for the annotation
            positive_column = f"{annotation}_positive"

            # Check if the column exists in self.spot_geodata
            if positive_column not in self.spot_geodata.columns:
                raise ValueError(f"Column '{positive_column}' does not exist in 'spot_geodata'.")

            

            # Get the 'plot_order' for the current annotation
            plot_order = self.data_exploded.loc[self.data_exploded['Annotation'] == annotation, 'plot_order'].iloc[0]

            # Assign annotations to spots based on positive status and plot_order
            for idx, row in self.spot_geodata.iterrows():
                if row[positive_column]:  # If the spot is positive for this annotation
                    # If the spot is already assigned a positive annotation, check if it should be replaced
                    if self.spot_geodata.at[idx, 'Annotation_map'] == "None":
                        self.spot_geodata.at[idx, 'Annotation_map'] = annotation
                    else:
                        # If multiple annotations are positive, keep the one with the lowest plot_order
                        current_annotation = self.spot_geodata.at[idx, 'Annotation_map']
                        current_plot_order = self.data_exploded.loc[self.data_exploded['Annotation'] == current_annotation, 'plot_order'].iloc[0]

                        if plot_order < current_plot_order:
                            self.spot_geodata.at[idx, 'Annotation_map'] = annotation

        # After all positive annotations have been processed, if annotate_all is True, handle spots that are still "None"
        if annotate_all:
            print("Annotating remaining spots with the most overlapping annotation...")

            for idx, row in self.spot_geodata.iterrows():
                if row['Annotation_map'] == "None":
                    # Find the most overlapping annotation for the spot
                    max_overlap = 0
                    best_annotation = "None"

                    # Loop over activated annotations and calculate the overlap
                    for annotation in self.activated_annotations:
                        overlap_column = f"{annotation}_overlap"
                        if overlap_column in self.spot_geodata.columns:
                            overlap = row[overlap_column]
                            if overlap > max_overlap:
                                max_overlap = overlap
                                best_annotation = annotation

                    # If the spot has no overlap with any annotation (max_overlap == 0), it stays "None"
                    if max_overlap > 0:
                        self.spot_geodata.at[idx, 'Annotation_map'] = best_annotation

        print('Annotation map successfully generated')


    def save(self, filename):
        """Save the current object to a pickle file."""
        try:
            import pickle
            with open(filename, 'wb') as f:
                pickle.dump(self, f)
            print(f"Object saved to {filename}")
        except Exception as e:
            print(f"Error saving object: {e}")
