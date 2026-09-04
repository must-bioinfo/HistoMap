"""Tests for histomaptx.distances.compute_nearest_annotation_distance."""

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import Point

from histomaptx.distances import compute_nearest_annotation_distance


class _FakeHistoMap:
    """Minimal stand-in for a HistoMap object.

    compute_nearest_annotation_distance only touches `self.spot_geodata`, so a real HistoMap
    (which needs a full SpatialData object, an image, a GeoJSON file, ...) isn't necessary to
    exercise it in isolation.
    """

    def __init__(self, spot_geodata):
        self.spot_geodata = spot_geodata


def _make_spots(annotations, geometries):
    return gpd.GeoDataFrame({"Annotation_map": annotations}, geometry=geometries)


def test_edge_is_the_default_method():
    # Centers 10 apart, radius 3 each -> edge-to-edge gap is 10 - 3 - 3 = 4.
    spot_a = Point(0, 0).buffer(3)
    spot_b = Point(10, 0).buffer(3)
    histomap = _FakeHistoMap(_make_spots(["Tumor", "Stroma"], [spot_a, spot_b]))

    compute_nearest_annotation_distance(histomap, "Tumor")

    assert histomap.spot_geodata.loc[1, "distance_to_Tumor"] == pytest.approx(4.0)


def test_centroid_method_ignores_spot_radius():
    spot_a = Point(0, 0).buffer(3)
    spot_b = Point(10, 0).buffer(3)
    histomap = _FakeHistoMap(_make_spots(["Tumor", "Stroma"], [spot_a, spot_b]))

    compute_nearest_annotation_distance(histomap, "Tumor", method="centroid")

    assert histomap.spot_geodata.loc[1, "distance_to_Tumor"] == pytest.approx(10.0)


def test_edge_method_accounts_for_spot_radius():
    # Centers 10 apart, radius 3 each -> edge-to-edge gap is 10 - 3 - 3 = 4.
    spot_a = Point(0, 0).buffer(3)
    spot_b = Point(10, 0).buffer(3)
    histomap = _FakeHistoMap(_make_spots(["Tumor", "Stroma"], [spot_a, spot_b]))

    compute_nearest_annotation_distance(histomap, "Tumor", method="edge")

    assert histomap.spot_geodata.loc[1, "distance_to_Tumor"] == pytest.approx(4.0)


def test_edge_method_is_zero_for_touching_or_overlapping_spots():
    spot_a = Point(0, 0).buffer(3)
    spot_b = Point(5, 0).buffer(3)  # overlaps spot_a
    histomap = _FakeHistoMap(_make_spots(["Tumor", "Stroma"], [spot_a, spot_b]))

    compute_nearest_annotation_distance(histomap, "Tumor", method="edge")

    assert histomap.spot_geodata.loc[1, "distance_to_Tumor"] == pytest.approx(0.0)


def test_invalid_method_raises():
    spot_a = Point(0, 0).buffer(3)
    spot_b = Point(10, 0).buffer(3)
    histomap = _FakeHistoMap(_make_spots(["Tumor", "Stroma"], [spot_a, spot_b]))

    with pytest.raises(ValueError, match="method must be 'centroid' or 'edge'"):
        compute_nearest_annotation_distance(histomap, "Tumor", method="bogus")


def test_raises_when_no_spot_has_the_annotation():
    spot_a = Point(0, 0).buffer(3)
    spot_b = Point(10, 0).buffer(3)
    histomap = _FakeHistoMap(_make_spots(["Stroma", "Stroma"], [spot_a, spot_b]))

    with pytest.raises(ValueError, match="No spots found with annotation 'Tumor'"):
        compute_nearest_annotation_distance(histomap, "Tumor")


def test_picks_the_nearest_of_several_annotated_spots():
    # method='centroid' here purely for round, easy-to-check numbers - this test is about the
    # "nearest of several" selection logic, not about which distance method is the default
    # (see test_edge_is_the_default_method for that).
    tumor_near = Point(0, 0).buffer(1)
    tumor_far = Point(100, 0).buffer(1)
    other = Point(5, 0).buffer(1)
    histomap = _FakeHistoMap(
        _make_spots(["Tumor", "Tumor", "Stroma"], [tumor_near, tumor_far, other])
    )

    compute_nearest_annotation_distance(histomap, "Tumor", method="centroid")

    assert histomap.spot_geodata.loc[2, "distance_to_Tumor"] == pytest.approx(5.0)


def test_annotated_spots_are_left_unset():
    spot_a = Point(0, 0).buffer(3)
    spot_b = Point(10, 0).buffer(3)
    histomap = _FakeHistoMap(_make_spots(["Tumor", "Stroma"], [spot_a, spot_b]))

    compute_nearest_annotation_distance(histomap, "Tumor")

    # The reference annotation's own spots never get a distance value written for them.
    assert pd.isna(histomap.spot_geodata.loc[0, "distance_to_Tumor"])


def test_composite_labels_do_not_exact_match_a_bare_annotation_name():
    # Annotation_map can hold ':'-joined composites (generate_annotation_map with
    # conserve_hierarchy=True); the lookup here is an exact match, so a spot labeled
    # "Tumor:Stroma" does not count as a 'Tumor' reference spot.
    spot_a = Point(0, 0).buffer(3)
    spot_b = Point(10, 0).buffer(3)
    histomap = _FakeHistoMap(_make_spots(["Tumor:Stroma", "Stroma"], [spot_a, spot_b]))

    with pytest.raises(ValueError, match="No spots found with annotation 'Tumor'"):
        compute_nearest_annotation_distance(histomap, "Tumor")
