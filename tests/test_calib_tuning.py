import pytest
from config.calib_tuning import CalibSamples, MCTunings, InvalidYearError


def test_calibration_samples():
    """Test the calibration samples."""
    calib_samples = CalibSamples()

    # hadrons
    assert calib_samples.fetch("2018", "hadron") == "Turbo18"
    assert calib_samples.fetch("2017", "hadron") == "Turbo17"
    assert calib_samples.fetch("2016", "hadron") == "Turbo16"
    assert calib_samples.fetch("2015", "hadron") == "Turbo15"
    assert calib_samples.fetch("2012", "hadron") == "21"
    assert calib_samples.fetch("2011", "hadron") == "21r1"

    # electrons
    assert calib_samples.fetch("2018", "e") == "Electron18"
    assert calib_samples.fetch("2017", "e") == "Electron17"
    assert calib_samples.fetch("2016", "e") == "Electron16"
    assert calib_samples.fetch("2015", "e") == "Electron15"
    assert calib_samples.fetch("2012", "e") == "21"
    assert calib_samples.fetch("2011", "e") == "21r1"

    # invalid year
    with pytest.raises(InvalidYearError):
        calib_samples.fetch("2020", "hadron")

    # invalid species
    with pytest.raises(ValueError):
        calib_samples.fetch("2018", "invalid")


def test_mc_tunings():
    """Test the MC tuning mapping."""

    # test the MC tunings
    mc_tunings = MCTunings()

    # test the MC tunings
    assert mc_tunings.fetch("2018") == "MC15TuneV1"
    assert mc_tunings.fetch("2017") == "MC15TuneV1"
    assert mc_tunings.fetch("2016") == "MC15TuneV1"
    assert mc_tunings.fetch("2015") == "MC15TuneV1"
    assert mc_tunings.fetch("2012") == "MC12TuneV4"
    assert mc_tunings.fetch("2011") == "MC12TuneV4"

    # invalid year
    with pytest.raises(InvalidYearError):
        mc_tunings.fetch("2020")
