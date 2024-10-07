__version__ = "0.1.0"

from .utils.io import read_config

from .auth import kinit
from .engine import _load_validated_config, config

# Attempt to load the pre-validated config when importing the package
try:
    _load_validated_config()
except FileNotFoundError:
    config = None
