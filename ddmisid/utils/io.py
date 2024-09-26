"""
Utility functions for plotting and path manipulation.
"""

from collections.abc import Callable
from functools import wraps
from pathlib import Path
from typing import TypeVar
from typing import Union, Any, List, Optional
import functools
import yaml
from typing_extensions import ParamSpec
import time
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
    """Ascertain correct path of binning config yml file."""

    @wraps(func)
    def inner(path: str, **kwargs: P.kwargs) -> R:
        if not Path(path).is_file():
            raise FileNotFoundError(f"Config file '{path}' not found.")
        return func(path, **kwargs)

    return inner


@check_argpath
def read_config(
    path: str,
    key: Union[None, str] = None,
) -> dict[str, Any] | Any:
    """Read in the feature from config yml file after checking it exists."""
    with open(path, "r") as stream:
        in_dict = yaml.load(stream, Loader=yaml.FullLoader)

    if key is not None:
        if key in in_dict:
            return in_dict[key]
        else:
            raise KeyError(f"'{key}' key not found in config file.")

    return in_dict


def load_root(
    path: str,
    key: Optional[str] = None,
    tree: Optional[str] = None,
    branches: Optional[List[str]] = None,
    cut: Optional[Union[List[str], str]] = None,
    name: Optional[str] = None,
    max_events: Optional[int] = None,
    batch_size: Optional[str] = "2 MB",
    library: str = "ak",  # 'pd' for pandas, 'ak' for awkward arrays
    **kwargs,
) -> Any:
    """Wrapper for uproot.iterate() to load ROOT files into pandas DataFrame or awkward Array."""

    # Open the ROOT file and tree
    if key is not None:
        events = uproot.open(f"{path}:{key}/{tree}")
    else:
        events = uproot.open(f"{path}:{tree}")

    # batched loa
    bevs = events.num_entries_for(batch_size, branches, entry_stop=max_events)
    tevs = events.num_entries
    nits = round(tevs / bevs + 0.5)

    # differentiate the library used
    if library == "pd":
        aggr = []
        for batch in tqdm(
            events.iterate(
                expressions=branches,
                cut=cut,
                library=library,
                entry_stop=max_events,
                step_size=batch_size,
                **kwargs,
            ),
            total=nits,
            ascii=True,
            desc=f"Batches loaded",
        ):
            aggr.append(batch)
        # Concatenate batches into one dataframe
        df = pd.concat(aggr)

        # Assign TeX label to DataFrame for plotting
        if name is not None:
            df.name = name

    # If awkward, batch and concatenate into awkward array
    elif library == "ak":
        aggr = []
        for batch in tqdm(
            events.iterate(
                expressions=branches,
                cut=cut,
                library=library,
                entry_stop=max_events,
                step_size=batch_size,
                **kwargs,
            ),
            total=nits,
            ascii=True,
            desc=f"Batches loaded",
        ):
            aggr.append(batch)

        # Concatenate batches into one awkward array along axis 0 (default)
        df = ak.concatenate(aggr, axis=0)

        # Assign TeX label to awkward array for plotting
        if name is not None:
            df.name = name

    # Raise an error if an unsupported library is passed
    else:
        raise ValueError(
            "Unsupported library. Use 'pd' for pandas or 'ak' for awkward arrays."
        )

    print(f"\nSUCCESS: loaded with {len(df)} entries")
    return df


@timing
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


@timing
def simple_load(
    path: str,
    key: str | None = None,
    tree: str | None = "DecayTree",
    branches: list[str] | None = None,
    max_events: int | None = None,
    library: str | None = "pd",
    cut: str | None = None,
    timeout: int | None = 500,
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
        events = uproot.open(f"{path}:{key}/{tree}", timeout=timeout)
    else:
        events = uproot.open(f"{path}:{tree}", timeout=timeout)

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