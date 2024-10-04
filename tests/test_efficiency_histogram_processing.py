import pytest
import numpy as np
import boost_histogram as bh
from hist import Hist
from ddmisid.utils.histogram import (
    NullEffProcessor,
    NegativeEffProcessor,
    EfficiencyHistogramProcessor,
    _epsilon,
)


# Fixture to create a mock histogram
@pytest.fixture
def mock_hist():
    hist = bh.Histogram(bh.axis.Regular(10, 0, 10), storage=bh.storage.Weight())
    for i in range(6):
        hist[i] = (i, i**2)  # Positive values
    hist[6] = (
        -0.001,
        0.0000000001,
    )  # Very small negative value, incompatible with zero, but within tolerance
    hist[7] = (0.0, 1.0)  # A null value
    hist[8] = (-1.0, 4.0)  # Negative value, large variance -> compatible with zero
    hist[9] = (-1.0, 0.1)  # Negative value, small variance -> not compatible with zero
    return hist


def test_null_eff_processor(mock_hist):
    """Test the NullEffProcessor to handle zero values."""
    null_processor = NullEffProcessor()
    processed_hist = null_processor.process(mock_hist)

    view = processed_hist.view()

    # Check that the null central value is changed to _epsilon
    assert view[7].value == _epsilon, "Null value not set to epsilon"

    # Check that variance remains unchanged unless it is zero
    assert (
        view[7].variance == 1.0
    ), "Variance should remain unchanged for non-zero values"


def test_negative_eff_processor(mock_hist):
    """Test the NegativeEffProcessor to handle negative values and raise errors for values outside tolerance."""
    negative_processor = NegativeEffProcessor()

    # Expect a ValueError to be raised due to the large negative value at index 9
    with pytest.raises(
        ValueError, match="PID efficiency value below 0.0 detected at index"
    ):
        negative_processor.process(mock_hist)

    # Process histogram again without the problematic index
    view = mock_hist.view()

    # Check that the negative value with large variance (index 8) is squashed to epsilon
    assert view[8].value == _epsilon, "Negative value not squashed to epsilon"
    # Check that the variance remains unchanged for compatible-with-zero cases
    assert (
        view[8].variance == 4.0
    ), "Variance incorrectly modified for large negative value"

    # Check that the small negative value is handled properly (index 6)
    assert view[6].value == _epsilon, "Small negative value not squashed to epsilon"
    assert (
        view[6].variance == 0.001
    ), "Variance incorrectly modified for small negative value"


def test_efficiency_hist_processor(mock_hist):
    """Test the combined EfficiencyHistogramProcessor."""
    combined_processor = EfficiencyHistogramProcessor()

    # Expect a ValueError to be raised due to the large negative value at index 9
    with pytest.raises(
        ValueError, match="PID efficiency value below 0.0 detected at index"
    ):
        combined_processor.process(mock_hist)

    # Check null value handling
    view = mock_hist.view()

    # Ensure that the small negative value (index 6) is processed correctly
    assert (
        view[6].value == _epsilon
    ), "Small negative value not squashed to epsilon by combined processor"
    assert (
        view[6].variance == 0.001
    ), "Variance not set correctly for small negative value by combined processor"

    # Ensure the negative value with large variance (index 8) is handled properly
    assert view[8].value == _epsilon, "Negative value not set to epsilon"
    assert view[8].variance == 4.0, "Variance incorrectly modified for negative value"
