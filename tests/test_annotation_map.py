"""Tests for HistoMap.set_positive / change_positive_threshold (overlap thresholding),
plot_order precedence, and generate_annotation_map (priority vs conserve_hierarchy).
"""

import pandas as pd
import pytest

from histomaptx.histomap_object import HistoMap


class _FakeHistoMap:
    """A lightweight stand-in for HistoMap, bound to the real methods under test.

    These methods only touch a handful of plain attributes (spot_geodata, positive_threshold,
    activated_annotations, data_exploded), so a fully constructed HistoMap is not needed.
    """

    set_positive = HistoMap.set_positive
    change_positive_threshold = HistoMap.change_positive_threshold
    generate_annotation_map = HistoMap.generate_annotation_map

    def __init__(self, spot_geodata, positive_threshold=None, activated_annotations=None,
                 data_exploded=None):
        self.spot_geodata = spot_geodata
        self.positive_threshold = positive_threshold
        self.activated_annotations = activated_annotations
        self.data_exploded = data_exploded


def _threshold_df(annotations, thresholds=50, overlay=True):
    n = len(annotations)
    thresholds = [thresholds] * n if isinstance(thresholds, (int, float)) else thresholds
    overlay = [overlay] * n if isinstance(overlay, bool) else overlay
    return pd.DataFrame({"annotation": annotations, "threshold": thresholds, "Overlay": overlay})


def _plot_order_df(annotation_to_order):
    return pd.DataFrame({
        "Annotation": list(annotation_to_order.keys()),
        "plot_order": list(annotation_to_order.values()),
    })


# ---------------------------------------------------------------------------
# set_positive
# ---------------------------------------------------------------------------

def test_set_positive_thresholds_overlap_column():
    spot_geodata = pd.DataFrame({"artery_overlap": [10, 50, 60, 90]})
    histomap = _FakeHistoMap(spot_geodata)

    histomap.set_positive("artery", 50)

    # >= threshold counts as positive, including the exact boundary value
    assert list(spot_geodata["artery_positive"]) == [False, True, True, True]


def test_set_positive_creates_column_if_missing():
    spot_geodata = pd.DataFrame({"Tumor_overlap": [0, 100]})
    histomap = _FakeHistoMap(spot_geodata)

    assert "Tumor_positive" not in spot_geodata.columns
    histomap.set_positive("Tumor", 50)
    assert "Tumor_positive" in spot_geodata.columns
    assert list(spot_geodata["Tumor_positive"]) == [False, True]


# ---------------------------------------------------------------------------
# change_positive_threshold
# ---------------------------------------------------------------------------

def test_change_positive_threshold_updates_threshold_and_recomputes_positive():
    spot_geodata = pd.DataFrame({"artery_overlap": [10, 40, 60]})
    positive_threshold = _threshold_df(["artery"], thresholds=50, overlay=True)
    histomap = _FakeHistoMap(spot_geodata, positive_threshold=positive_threshold)

    histomap.change_positive_threshold({"artery": 30})

    assert histomap.positive_threshold.loc[0, "threshold"] == 30
    assert list(spot_geodata["artery_positive"]) == [False, True, True]


def test_change_positive_threshold_updates_several_annotations_at_once():
    spot_geodata = pd.DataFrame({
        "artery_overlap": [10, 60],
        "Tumor_overlap": [20, 80],
    })
    positive_threshold = _threshold_df(["artery", "Tumor"], thresholds=50, overlay=True)
    histomap = _FakeHistoMap(spot_geodata, positive_threshold=positive_threshold)

    histomap.change_positive_threshold({"artery": 5, "Tumor": 90})

    assert list(spot_geodata["artery_positive"]) == [True, True]
    assert list(spot_geodata["Tumor_positive"]) == [False, False]


def test_change_positive_threshold_raises_for_unknown_annotation():
    positive_threshold = _threshold_df(["artery"], overlay=True)
    histomap = _FakeHistoMap(pd.DataFrame({"artery_overlap": [10]}), positive_threshold=positive_threshold)

    with pytest.raises(ValueError, match="Annotation 'Tumor' not found"):
        histomap.change_positive_threshold({"Tumor": 50})


@pytest.mark.parametrize("bad_threshold", [-1, 101, 150])
def test_change_positive_threshold_raises_for_out_of_range_threshold(bad_threshold):
    positive_threshold = _threshold_df(["artery"], overlay=True)
    histomap = _FakeHistoMap(pd.DataFrame({"artery_overlap": [10]}), positive_threshold=positive_threshold)

    with pytest.raises(ValueError, match="must be between 0 and 100"):
        histomap.change_positive_threshold({"artery": bad_threshold})


def test_change_positive_threshold_raises_when_overlay_not_computed():
    positive_threshold = _threshold_df(["artery"], overlay=False)
    histomap = _FakeHistoMap(pd.DataFrame({"artery_overlap": [10]}), positive_threshold=positive_threshold)

    with pytest.raises(ValueError, match="should be computed first"):
        histomap.change_positive_threshold({"artery": 50})


# ---------------------------------------------------------------------------
# generate_annotation_map: plot_order precedence (conserve_hierarchy=False, the default)
# ---------------------------------------------------------------------------

def test_priority_lowest_plot_order_wins():
    # artery has the lowest plot_order (0) -> highest precedence, should win over Stroma (1)
    # and Tumor (2) wherever it's positive too.
    spot_geodata = pd.DataFrame({
        "artery_positive": [True, False, True],
        "Stroma_positive": [True, True, False],
        "Tumor_positive": [False, True, False],
    })
    histomap = _FakeHistoMap(
        spot_geodata,
        positive_threshold=_threshold_df(["artery", "Stroma", "Tumor"], overlay=True),
        activated_annotations=["artery", "Stroma", "Tumor"],
        data_exploded=_plot_order_df({"artery": 0, "Stroma": 1, "Tumor": 2}),
    )

    histomap.generate_annotation_map(annotate_all=False, conserve_hierarchy=False)

    # row 0: artery only -> artery. row 1: Stroma(1) and Tumor(2) -> Stroma wins. row 2: artery only -> artery.
    assert list(spot_geodata["Annotation_map"]) == ["artery", "Stroma", "artery"]


def test_priority_is_independent_of_activated_annotations_order():
    # Same data as above, but activated_annotations is processed in a different order; the
    # winner should still be determined purely by plot_order, not iteration order.
    spot_geodata = pd.DataFrame({
        "artery_positive": [True, False, True],
        "Stroma_positive": [True, True, False],
        "Tumor_positive": [False, True, False],
    })
    histomap = _FakeHistoMap(
        spot_geodata,
        positive_threshold=_threshold_df(["artery", "Stroma", "Tumor"], overlay=True),
        activated_annotations=["Tumor", "artery", "Stroma"],  # shuffled on purpose
        data_exploded=_plot_order_df({"artery": 0, "Stroma": 1, "Tumor": 2}),
    )

    histomap.generate_annotation_map(annotate_all=False, conserve_hierarchy=False)

    assert spot_geodata.loc[0, "Annotation_map"] == "artery"
    assert spot_geodata.loc[2, "Annotation_map"] == "artery"


def test_priority_single_positive_annotation_wins_by_default():
    spot_geodata = pd.DataFrame({
        "artery_positive": [False],
        "Stroma_positive": [True],
    })
    histomap = _FakeHistoMap(
        spot_geodata,
        positive_threshold=_threshold_df(["artery", "Stroma"], overlay=True),
        activated_annotations=["artery", "Stroma"],
        data_exploded=_plot_order_df({"artery": 0, "Stroma": 1}),
    )

    histomap.generate_annotation_map(annotate_all=False, conserve_hierarchy=False)

    assert spot_geodata.loc[0, "Annotation_map"] == "Stroma"


def test_generate_annotation_map_raises_for_unregistered_annotation():
    histomap = _FakeHistoMap(
        pd.DataFrame({"artery_positive": [True]}),
        positive_threshold=_threshold_df([], overlay=True),  # empty: 'artery' missing
        activated_annotations=["artery"],
        data_exploded=_plot_order_df({"artery": 0}),
    )

    with pytest.raises(ValueError, match="'artery' not found in the object"):
        histomap.generate_annotation_map()


def test_generate_annotation_map_raises_when_overlay_not_computed():
    histomap = _FakeHistoMap(
        pd.DataFrame({"artery_positive": [True]}),
        positive_threshold=_threshold_df(["artery"], overlay=False),
        activated_annotations=["artery"],
        data_exploded=_plot_order_df({"artery": 0}),
    )

    with pytest.raises(ValueError, match="Overlap for 'artery' was not computed"):
        histomap.generate_annotation_map()


def test_generate_annotation_map_raises_when_positive_column_missing():
    histomap = _FakeHistoMap(
        pd.DataFrame({"unrelated_column": [1]}),
        positive_threshold=_threshold_df(["artery"], overlay=True),
        activated_annotations=["artery"],
        data_exploded=_plot_order_df({"artery": 0}),
    )

    with pytest.raises(ValueError, match="Column 'artery_positive' does not exist"):
        histomap.generate_annotation_map()


# ---------------------------------------------------------------------------
# generate_annotation_map: conserve_hierarchy=True
# ---------------------------------------------------------------------------

def test_conserve_hierarchy_combines_positive_annotations_ordered_by_plot_order():
    spot_geodata = pd.DataFrame({
        "artery_positive": [True],
        "Stroma_positive": [True],
        "Tumor_positive": [False],
    })
    histomap = _FakeHistoMap(
        spot_geodata,
        positive_threshold=_threshold_df(["artery", "Stroma", "Tumor"], overlay=True),
        activated_annotations=["artery", "Stroma", "Tumor"],
        data_exploded=_plot_order_df({"artery": 0, "Stroma": 1, "Tumor": 2}),
    )

    histomap.generate_annotation_map(annotate_all=False, conserve_hierarchy=True)

    assert spot_geodata.loc[0, "Annotation_map"] == "artery:Stroma"


def test_conserve_hierarchy_orders_composite_by_plot_order_not_activation_order():
    spot_geodata = pd.DataFrame({
        "artery_positive": [True],
        "Stroma_positive": [True],
    })
    histomap = _FakeHistoMap(
        spot_geodata,
        positive_threshold=_threshold_df(["artery", "Stroma"], overlay=True),
        activated_annotations=["Stroma", "artery"],  # listed in reverse of plot_order
        data_exploded=_plot_order_df({"artery": 0, "Stroma": 1}),
    )

    histomap.generate_annotation_map(annotate_all=False, conserve_hierarchy=True)

    # plot_order (artery=0 before Stroma=1) wins over activated_annotations iteration order
    assert spot_geodata.loc[0, "Annotation_map"] == "artery:Stroma"


def test_conserve_hierarchy_single_positive_has_no_colon():
    spot_geodata = pd.DataFrame({
        "artery_positive": [True],
        "Stroma_positive": [False],
    })
    histomap = _FakeHistoMap(
        spot_geodata,
        positive_threshold=_threshold_df(["artery", "Stroma"], overlay=True),
        activated_annotations=["artery", "Stroma"],
        data_exploded=_plot_order_df({"artery": 0, "Stroma": 1}),
    )

    histomap.generate_annotation_map(annotate_all=False, conserve_hierarchy=True)

    assert spot_geodata.loc[0, "Annotation_map"] == "artery"


def test_conserve_hierarchy_false_vs_true_differ_only_on_multi_positive_spots():
    spot_geodata_a = pd.DataFrame({
        "artery_positive": [True, False],
        "Stroma_positive": [True, True],
    })
    spot_geodata_b = spot_geodata_a.copy()

    common_kwargs = dict(
        positive_threshold=_threshold_df(["artery", "Stroma"], overlay=True),
        activated_annotations=["artery", "Stroma"],
        data_exploded=_plot_order_df({"artery": 0, "Stroma": 1}),
    )

    winner_take_all = _FakeHistoMap(spot_geodata_a, **common_kwargs)
    winner_take_all.generate_annotation_map(annotate_all=False, conserve_hierarchy=False)

    combined = _FakeHistoMap(spot_geodata_b, **common_kwargs)
    combined.generate_annotation_map(annotate_all=False, conserve_hierarchy=True)

    # spot 0 is positive for both -> differs between modes
    assert spot_geodata_a.loc[0, "Annotation_map"] == "artery"
    assert spot_geodata_b.loc[0, "Annotation_map"] == "artery:Stroma"
    # spot 1 is positive for only Stroma -> identical in both modes
    assert spot_geodata_a.loc[1, "Annotation_map"] == "Stroma"
    assert spot_geodata_b.loc[1, "Annotation_map"] == "Stroma"


# ---------------------------------------------------------------------------
# generate_annotation_map: annotate_all fallback (applies regardless of conserve_hierarchy)
# ---------------------------------------------------------------------------

def test_annotate_all_true_assigns_most_overlapping_annotation_to_unmatched_spots():
    spot_geodata = pd.DataFrame({
        "artery_positive": [False],
        "Stroma_positive": [False],
        "artery_overlap": [10],
        "Stroma_overlap": [30],
    })
    histomap = _FakeHistoMap(
        spot_geodata,
        positive_threshold=_threshold_df(["artery", "Stroma"], overlay=True),
        activated_annotations=["artery", "Stroma"],
        data_exploded=_plot_order_df({"artery": 0, "Stroma": 1}),
    )

    histomap.generate_annotation_map(annotate_all=True, conserve_hierarchy=False)

    assert spot_geodata.loc[0, "Annotation_map"] == "Stroma"  # 30 > 10


def test_annotate_all_true_leaves_zero_overlap_spots_as_none():
    spot_geodata = pd.DataFrame({
        "artery_positive": [False],
        "artery_overlap": [0],
    })
    histomap = _FakeHistoMap(
        spot_geodata,
        positive_threshold=_threshold_df(["artery"], overlay=True),
        activated_annotations=["artery"],
        data_exploded=_plot_order_df({"artery": 0}),
    )

    histomap.generate_annotation_map(annotate_all=True, conserve_hierarchy=False)

    assert spot_geodata.loc[0, "Annotation_map"] == "None"


def test_annotate_all_false_leaves_unmatched_spots_as_none():
    spot_geodata = pd.DataFrame({
        "artery_positive": [False],
        "artery_overlap": [80],  # would win under annotate_all=True, but it's disabled here
    })
    histomap = _FakeHistoMap(
        spot_geodata,
        positive_threshold=_threshold_df(["artery"], overlay=True),
        activated_annotations=["artery"],
        data_exploded=_plot_order_df({"artery": 0}),
    )

    histomap.generate_annotation_map(annotate_all=False, conserve_hierarchy=False)

    assert spot_geodata.loc[0, "Annotation_map"] == "None"
