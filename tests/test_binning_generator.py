import pytest
from pathlib import Path
import json
from ddmisid.utils import DefaultBinningGenerator, PIDCalibAliasFactory


@pytest.fixture
def setup_binning():
    """Test binning and species dictionaries."""
    example_species = {
        "kaon": "K",
    }

    example_binning = {
        "P": [10000.0, 20000.0, 30000.0, 40000.0, 100000.0],
        "PT": [500.0, 1000.0, 1500.0, 2000.0, 4000.0],
        "ETA": [2.0, 2.5, 3.0, 3.5, 4.0],
        "nTracks": [0, 100, 200, 300, 400],
    }

    binning_generator = DefaultBinningGenerator(
        species="kaon",
        species_alias="K",
        binning=example_binning,
        binning_alias="TestBinning",
    )

    yield binning_generator

    # Clean up after test
    output_path = Path("tests/data/binning_2016/TestBinning.json")
    if output_path.exists():
        output_path.unlink()


def test_process_variable():
    """Verify PIDCalib2 alias compliance."""
    assert PIDCalibAliasFactory.process_variable("P", "2011", "K") == "P"
    assert PIDCalibAliasFactory.process_variable("PT", "2018", "K") == "Brunel_PT"
    assert (
        PIDCalibAliasFactory.process_variable("nTracks", "2016", "K")
        == "nTracks_Brunel"
    )
    assert (
        PIDCalibAliasFactory.process_variable("ETA", "2016", "e_B_Jpsi") == "Brunel_ETA"
    )
    assert PIDCalibAliasFactory.process_variable("ETA", "2012", "e_B_Jpsi") == "ETA"
    assert (
        PIDCalibAliasFactory.process_variable("nTracks", "2011", "e_B_Jpsi")
        == "nTracks"
    )


def test_binning_build(setup_binning):
    """Test binning build method."""
    expected_binning = {
        "K": {
            "Brunel_P": [10000.0, 20000.0, 30000.0, 40000.0, 100000.0],
            "Brunel_PT": [500.0, 1000.0, 1500.0, 2000.0, 4000.0],
            "Brunel_ETA": [2.0, 2.5, 3.0, 3.5, 4.0],
            "nTracks_Brunel": [0, 100, 200, 300, 400],
        }
    }

    # Generate binning and write to JSON file
    setup_binning.build(year="2016", outdir="tests/data", verbose=False)

    # Read JSON file back in, and compare to expected binning
    with open(Path("tests/data/binning_2016/kaon_TestBinning.json")) as f:
        binning = json.load(f)

    assert binning == expected_binning


def test_binning_variable_compliance(setup_binning):
    """Ensure the binning variables comply with the PIDCalib2 calibration-sample guidelines."""
    # run 2 aliases
    binvars = setup_binning.get_binning_variables(year="2016")
    assert binvars == [
        "Brunel_P",
        "Brunel_PT",
        "Brunel_ETA",
        "nTracks_Brunel",
    ]

    # run 1 aliases
    binvars = setup_binning.get_binning_variables(year="2012")
    assert binvars == ["P", "PT", "ETA", "nTracks"]
