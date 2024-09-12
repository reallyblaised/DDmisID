"""
Utility functions for plotting and path manipulation.
"""

__author__ = "Blaise Delaney"
__email__ = "Blaise Delaney at cern.ch"

from collections.abc import Callable
from functools import wraps
from pathlib import Path
from typing import TypeVar
from typing import Union, Any, List, Optional
import functools
import yaml
from typing_extensions import ParamSpec
import time
from termcolor2 import c as tc
import pandas as pd
import uproot
from tqdm import tqdm
import awkward as ak
from numpy.typing import ArrayLike
import re
import pickle

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
        print((f"{func.__name__!r} execution completed in {run_time:.4f} secs\n"))
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
    key: Optional[str] = None,
    tree_name: Optional[str] = None,
    branches: Optional[List[str]] = None,
    library: str = "ak",  # default to awkward
    cut: Optional[List[str]] = None,
    name: Optional[str] = None,
    max_entries: Optional[int] = None,
    batch_size: Optional[str] = "50 KB",
    **kwargs,
) -> Any:
    """Load file using pkl or uproot, depending on file extension"""
    ext = Path(file_path).suffix
    if ext == ".pkl":
        df = pd.read_pickle(file_path)
    elif ext == ".root":
        df = load_root(
            file_path=file_path,
            key=key,
            tree_name=tree_name,
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
    file_path: str,
    library: str,
    key: Optional[str] = None,
    tree_name: Optional[str] = None,
    branches: Optional[List[str]] = None,
    cut: Optional[Union[List[str], str]] = None,
    name: Optional[str] = None,
    max_entries: Optional[int] = None,
    batch_size: Optional[str] = "200 MB",
    **kwargs,
) -> Any:
    """Wrapper for uproot.iterate() to load ROOT files into a pandas DataFrame"""

    if key is not None:
        events = uproot.open(f"{file_path}:{key}/{tree_name}")
    else:
        events = uproot.open(f"{file_path}:{tree_name}")

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


def open_count(
    file_path: str,
    key: str | None = None,
    tree_name: str | None = None,
) -> int:
    """Open ROOT file and simply count the number of entries"""
    if key is not None:
        events = uproot.open(f"{file_path}:{key}/{tree_name}")
    else:
        events = uproot.open(f"{file_path}:{tree_name}")

    return events.num_entries


def simple_load(
    path: str,
    key: str | None = None,
    tree: str | None = "DecayTree",
    branches: list[str] | None = None,
    max_events: int | None = None,
    library: str | None = "pd",
    cut: str | None = None,
    **kwargs: Any,
) -> Any:
    """Load a pandas DataFrame from a ROOT file.

    Parameters
    ----------
    path : str
        Path to the ROOT file.
    key : str | None
        Key to target directory containing the tree.
    tree : str | None
        Name of the tree to load.
    branches : list[str] | None
        List of branches to load. If None (default), load all branches.
    max_events : int | None
        Maximum number of events to load. If None (default), load all events.
    library : str | None
        Library to use for loading the data. If None (default), use akward arrays.
    cut : str | None
        Cut to apply to the data. If None (default), no cut is applied.

    Returns
    -------
    Any
        Loaded data. If library is None, return awkward arrays. Default is pandas DataFrame.

    """
    if key is not None:
        events = uproot.open(f"{path}:{key}/{tree}")
    else:
        events = uproot.open(f"{path}:{tree}")

    # load into pandas DataFrame
    return events.arrays(
        library=library, expressions=branches, entry_stop=max_events, cut=cut, **kwargs
    )


def book_output_file(
    path: str,
    compression: uproot.compression.Compression = uproot.ZLIB(4),
) -> uproot.WritableDirectory:
    """Book the output file with default compression ZLIB(4)

    Parameters
    ----------
    path : str
        path to output file
    compression : uproot.compression.Compression
        compression algorithm to use (default ZLIB(4))

    Returns
    -------
    outfile : uproot.rootio.ROOTDirectory
        Path to output file
    """
    outpath = Path(path)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    outfile = uproot.recreate(outpath.absolute(), compression=compression)

    return outfile


def update_write_df(
    outfile: uproot.WritableDirectory,
    data: pd.DataFrame | ak.Array,
    key: str | None = None,
    treename: str = "DecayTree",
) -> None:
    """Update a ROOT file to include the new data

    Parameters
    ----------
    data : pd.DataFrame | ak.Array
        data to write out

    outfile : uproot.rootio.ROOTDirectory
        Booked outfile

    key : str
        name of the output directory

    treename : str
        name of the tree to write out
    """
    # handling the directory structure, if any
    if key is not None:
        outfile[f"{key}/{treename}"] = data
        # display what we have written out
        print(
            f"Entries written to file {key}/{treename}: ",
            outfile[f"{key}/{treename}"].num_entries,
        )

    else:
        outfile[f"{treename}"] = data
        # display what we have written out
        print(
            f"Entries written to file DecayTree: ",
            outfile[f"{treename}"].num_entries,
        )


def write_df(
    data: pd.DataFrame | ak.Array,
    path: str,
    key: str | None = None,
    treename: str = "DecayTree",
) -> None:
    """Write anew (recreate) the dataframe to a ROOT file

    Parameters
    ----------
    data : pd.DataFrame | ak.Array
        data to write out

    path : str
        path to output file

    key : str
        name of the output directory

    treename : str
        name of the tree to write out
    """
    # book output file
    outfile = book_output_file(path)
    # write anew
    update_write_df(outfile, data, key, treename)


def extract_sel_dict_branches(selection_dict: dict) -> list:
    """
    Extracts unique branches from the selection conditions in the provided dictionary,
    combining all branches found across all keys into a single list without duplicates.

    Parameters
    ----------
        selection_dict (dict): A dictionary with values containing selection strings.

    Returns
    -------
        list: A list of unique branch names used across all selections.
    """
    branch_pattern = re.compile(r"\b[A-Za-z_]+\b(?=\s*[!<>=])")
    all_branches = set()

    for selection in selection_dict.values():
        branches = branch_pattern.findall(selection)
        all_branches.update(branches)

    return list(all_branches)


def load_hist(f):
    """Load boost_histogram/Hist from pickle file"""
    with open(f, "rb") as f_in:
        return pickle.load(f_in)
