from setuptools import setup, find_packages

setup(
    name="ddmisid",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pidcalib2",
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
        "xrootd @ https://test-files.pythonhosted.org/packages/f9/0b/db7c22839324172286cddb2c5dbd5cc3aaf51c78e455d3d57d23744013f3/xrootd-20180823104.tar.gz",
        "colorama",
        "termcolor",
        "tabulate",
        "pyyaml",
        "xxhash",
    ],
)
