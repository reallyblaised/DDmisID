"""
Utility functions for plotting and path manipulation.
"""
__author__ = "Blaise Delaney"
__email__ = "Blaise Delaney at cern.ch"

from collections.abc import Callable
from functools import wraps
from pathlib import Path
from typing import TypeVar
from typing import Union, Any
import functools
import yaml
from typing_extensions import ParamSpec
import time
from termcolor2 import c as tc
import pandas as pd
import uproot
from tqdm import tqdm

P = ParamSpec("P")
R = TypeVar("R")
T = TypeVar("T")


def debug(func: Callable[P, T]) -> Callable[P, T]:
    """Print the function signature"""

    @functools.wraps(func)
    def wrapper_debug(*args: P.args, **kwargs: P.kwargs) -> T:
        args_repr = [repr(a) for a in args]
        kwargs_repr = [f"{k}={v!r}" for k, v in kwargs.items()]
        signature = ", ".join(args_repr + kwargs_repr)
        print(f"Calling {func.__name__}({signature})")
        return func(*args, **kwargs)

    return wrapper_debug


def timing(func: Callable[P, T]) -> Callable[P, T]:
    """Print the runtime of the decorated function"""

    @functools.wraps(func)
    def wrapper_timer(*args: P.args, **kwargs: P.kwargs) -> T:
        start_time = time.perf_counter()  # 1
        value = func(*args, **kwargs)
        end_time = time.perf_counter()  # 2
        run_time = end_time - start_time  # 3
        print(tc(f"Success: finished {func.__name__!r} in {run_time:.4f} secs\n").green)
        return value

    return wrapper_timer


def check_argpath(func: Callable[P, R]) -> Callable[P, R]:
    """Ascertain correct path of binning config yml file"""

    @wraps(func)
    def inner(path: str, **kwargs: P.kwargs) -> R:
        try:
            Path(path)
        except IOError:
            print("Incorrect input path")
        features = func(path, **kwargs)
        return features

    return inner


@check_argpath
def read_config(
    path: str,
    key: Union[None, str] = None,
) -> Any:
    """Read in the feature from config yml file after checking it exists"""

    with open(path, "r") as stream:
        in_dict = yaml.load(stream, Loader=yaml.FullLoader)
    if key is not None:
        try:
            key in in_dict
            feature_list = in_dict[key]
        except ValueError:
            print(f"'{key}' key not in dict")

    return feature_list


@timing
def load_ntuple(
    file_path: str,
    key: str,
    tree_name: str = "DecayTree",
    branches: list[str] | None = None,
    library: str = "ak",  # default to awkward
    cut: list[str] | str | None = None,
    name: str | None = None,
    max_entries: int | None = None,
    batch_size: str | None = "50 KB",
    **kwargs,
) -> Any:
    """Load file using pkl or uproot, depending on file extension"""
    ext = Path(file_path).suffix
    if ext == ".pkl":
        df = pd.read_pickle(file_path)
    elif ext == ".root":
        file_tree_path = f"{file_path}:{key}/{tree_name}"
        df = load_root(
            file=file_tree_path,
            library=library,
            branches=branches,
            cut=cut,
            name=name,
            max_entries=max_entries,
            batch_size=batch_size,
            **kwargs,
        )
    else:
        raise ValueError("File extension not recognised")

    return df


def load_root(
    file: str,
    library: str,
    branches: list[str] | None = None,
    cut: list[str] | None = None,
    name: str | None = None,
    max_entries: int | None = None,
    batch_size: str | None = "200 MB",
    **kwargs,
) -> pd.DataFrame:
    """Wrapper for uproot.iterate() to load ROOT files into a pandas DataFrame"""

    events = uproot.open(f"{file}")

    # if pandas, batch and concatenate
    if library == "pd":
        bevs = events.num_entries_for(batch_size, branches, entry_stop=max_entries)
        tevs = events.num_entries
        nits = round(tevs / bevs + 0.5)
        aggr = []
        for batch in tqdm(
            events.iterate(
                expressions=branches,
                cut=cut,
                library=library,
                entry_stop=max_entries,
                step_size=batch_size,
                **kwargs,
            ),
            total=nits,
            ascii=True,
            desc=f"Batches loaded",
        ):
            aggr.append(batch)
        # concatenate batches into one dataframe
        df = pd.concat(aggr)

        # assign TeX label to df for plotting
        if name is not None:
            df.name = name

    # else, load into awkward or numpy objects
    else:
        df = events.arrays(
            expressions=branches,
            cut=cut,
            library=library,
            entry_stop=max_entries,
            **kwargs,
        )

    print(f"\nSUCCESS: loaded with {len(df)} entries")
    return df
