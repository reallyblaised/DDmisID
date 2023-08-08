"""
Unit testing for make_templates.py
"""

import pytest
import ddmisid.make_templates
import numpy as np
import pickle
import os
import shutil
import hist


np.random.seed(42)


def create_template_maker():
    """
    Create instance of ParticleTemplateMaker class without eff hists.
    """
    template_maker = ddmisid.make_templates.HadronTemplateMaker("none", "0", [])
    return template_maker


def populate_template_maker(template_maker, reco_categories, eff_hists):
    template_maker.template_bins = [ f"{reco}" for reco in reco_categories ]
    template_maker.eff_hists = eff_hists


def get_test_data(dims, num_reco_categories, get_bins=False):
    def make_eff_hist(dims):
        num_points = 1000
        data = [ np.random.normal(loc=10, scale=2, size=num_points) for _ in range(dims) ]
        hist_axes = [ hist.axis.Variable(edges=np.linspace(0, 20, 21)) for _ in range(dims) ]
        h = hist.Hist(*hist_axes, storage=hist.storage.Weight())
        h.fill(*data, weight=0.5)
        return h
    reco_categories = [ f"0_to_{reco}_like" for reco in range(num_reco_categories) ]
    eff_hists = {}
    for reco in reco_categories:
        eff_hists[reco] = make_eff_hist(dims)
    if not get_bins:
        return create_template_maker(), reco_categories, eff_hists
    else:
        return create_template_maker(), reco_categories, eff_hists, {
            reco: [ num for num in range(21) ] for reco in reco_categories
        }


def open_hist(hist_path):
    with open(f"{hist_path}", "rb") as file:
        eff_hist = pickle.load(file)
    return eff_hist


def test_template_maker_creation():
    template_maker = create_template_maker()
    assert isinstance(template_maker, ddmisid.make_templates.HadronTemplateMaker), "incorrect instantiation"


@pytest.mark.parametrize(
    "template_maker,reco_categories,eff_hists",
    [ get_test_data(2, 4), get_test_data(3, 5), get_test_data(3, 2), get_test_data(2, 7) ]
)
def test_eff_hist_to_hist(template_maker, reco_categories, eff_hists):
    populate_template_maker(template_maker, reco_categories, eff_hists)
    res = template_maker.make_hist()
    assert isinstance(res, hist.Hist), "template maker is not creating or returning a histogram"

    for reco in reco_categories:
        assert eff_hists[reco] == res[..., reco], "self.hist histogram not filled correctly"


@pytest.mark.parametrize(
    "template_maker,reco_categories,eff_hists", 
    [ get_test_data(4, 4), get_test_data(1, 5), get_test_data(4, 2), get_test_data(5, 3) ]
)
def test_illegal_input(template_maker, reco_categories, eff_hists):
    populate_template_maker(template_maker, reco_categories, eff_hists)
    try:
        template_maker.make_hist()
    except ValueError:
        pass
    except AssertionError:
        pass
    else:
        raise AttributeError("illegal input is not caught")


# @pytest.mark.parametrize(
#     "template_maker,reco_categories,eff_hists,binning",
#     [ get_test_data(2, 4, True), get_test_data(3, 5, True), get_test_data(3, 2, True), get_test_data(2, 7, True) ]
# )
# def test_template_saving(template_maker, reco_categories, eff_hists, binning):
#     populate_template_maker(template_maker, reco_categories, eff_hists)
#     template_maker.make_hist()
#     template_maker.save_templates(binning, path="tests/test_templates")

#     template_file_paths = []
#     for (_, _, file_names) in os.walk("tests/test_templates"):
#         template_file_paths.extend(file_names)

#     num_bins = 1
#     for _, binning_axis in binning.items():
#         num_bins *= len(binning_axis)-1
#     assert len(template_file_paths) == num_bins, "incorrect number of files generated"
#     shutil.rmtree("tests/test_templates")
    


