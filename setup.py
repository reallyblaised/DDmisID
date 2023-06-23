from setuptools import setup, find_packages

setup(
    name="ddmisid",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
        "boost-histogram",
        "hist",
        "uproot",
        "iminuit<2.21",
        "numba",
        "pyhf",
        "tqdm",
        "colorama",
        "termcolor",
        "tabulate",
        "pyyaml",
        "snakemake",
        "termcolor2",
    ],
)
