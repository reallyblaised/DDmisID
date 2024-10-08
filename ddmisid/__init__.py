__version__ = "0.1.0"

from .utils.io import (
    read_config,
    load_root,
    simple_load,
    write_df,
    update_write_df,
    extract_sel_dict_branches,
)
from .utils.plot import (
    data_pull_plot,
    make_legend,
    save_to,
)
from .auth import kinit
from .engine import _load_validated_config, config
from .utils.histogram import load_hist, EfficiencyHistogramProcessor

# Attempt to load the pre-validated config when importing the package
try:
    _load_validated_config()
except FileNotFoundError:
    config = None
