import pytest
from ddmisid.bml_fitter import BMLFitter
import numpy as np
import json
import pickle
import yaml
from pathlib import Path
import hist


# Mock configuration and data for use in the tests
def mock_config():
    h1 = hist.Hist(hist.axis.Regular(3, 0, 3, name="x"))
    h1.fill(x=[0.5, 1.5, 2.5], weight=[1, 2, 3])

    h2 = hist.Hist(hist.axis.Regular(3, 0, 3, name="x"))
    h2.fill(x=[0.5, 1.5, 2.5], weight=[4, 5, 6])

    return {
        "template1": h1,
        "template2": h2,
    }


def mock_data():
    return np.array([10, 20, 30])


# A pytest fixture for creating a BMLFitter instance
@pytest.fixture
def mock_fitter(tmp_path="test_data"):
    # Create temporary files for data and config
    Path(tmp_path).mkdir(exist_ok=True, parents=True)
    data_file = f"{tmp_path}/data.pkl"
    config_file = f"{tmp_path}/config.yml"

    # Save the mock data and config to the temporary files
    with open(data_file, "wb") as file:
        pickle.dump(mock_data(), file)
    with open(config_file, "w") as file:
        yaml.dump(mock_config(), file)

    # Create and return the BMLFitter instance
    return BMLFitter(str(data_file), str(config_file))


# Test the initialization of the BMLFitter class
def test_bml_fitter_init(mock_fitter):
    assert isinstance(mock_fitter.config, dict)
    assert isinstance(mock_fitter.data, np.ndarray)


# Test the _normalise_templates method
def test_normalise_templates(mock_fitter):
    mock_fitter._normalise_templates()
    for template_name in mock_fitter.config.keys():
        template = getattr(mock_fitter, template_name)
        assert np.isclose(np.sum(template), 1.0)


# Test the build_workspace method
def test_build_workspace(mock_fitter):
    spec = mock_fitter.build_workspace()
    assert isinstance(spec, dict)
    assert "channels" in spec


# Test the export_workspace method
def test_export_workspace(mock_fitter, tmp_path):
    filename = tmp_path / "workspace.json"
    mock_fitter.build_workspace()
    mock_fitter.export_workspace(str(filename))
    with open(filename) as file:
        exported_spec = json.load(file)
    assert isinstance(exported_spec, dict)
    assert "channels" in exported_spec


# Add additional test functions here to cover other methods and functionality
# Test the fit method
def test_fit(mock_fitter):
    result = mock_fitter.fit()
    assert isinstance(result, np.ndarray)
    assert "fit_result" in mock_fitter.result_dict


# Test the visualize method (without actual plotting)
def test_visualize(mock_fitter, monkeypatch):
    # We're not testing the actual plot, so we'll patch plt.show to avoid displaying it
    monkeypatch.setattr(plt, "show", lambda: None)
    mock_fitter.fit()
    # Just call the method to ensure it does not raise an exception
    mock_fitter.visualize()


# Test the save_results method (without saving to a file)
def test_save_results_no_save(mock_fitter):
    mock_fitter.fit()
    result_dict = mock_fitter.save_results(save=False)
    assert "fit_result" in result_dict


# Test the save_results method (with saving to a file)
def test_save_results_save(mock_fitter, tmp_path):
    filename = tmp_path / "results.json"
    mock_fitter.fit()
    result_dict = mock_fitter.save_results(save=True, filename=str(filename))
    with open(filename) as file:
        saved_result = json.load(file)
    assert "fit_result" in saved_result


def test_incorrect_template_type(tmp_path):
    # Create temporary files for data and config
    data_file = tmp_path / "data.pkl"
    config_file = tmp_path / "config.yml"

    incorrect_config = {
        "template1": "incorrect_type",
    }

    # Save the incorrect config and mock data to temporary files
    with open(data_file, "wb") as file:
        pickle.dump(mock_data(), file)
    with open(config_file, "w") as file:
        yaml.dump(incorrect_config(), file)

    # The initialization should raise a ValueError due to the incorrect template type
    with pytest.raises(ValueError):
        BMLFitter(str(data_file), str(config_file))
