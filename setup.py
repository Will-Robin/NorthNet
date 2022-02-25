from setuptools import setup, version

setup(
    name="NorthNet",
    version="0.2",
    author="William E. Robinson",
    packages = ["NorthNet"],
    install_requires=[
        "networkx >= 2.6.3",
        "numpy >= 1.21.5",
        "rdkit >= 2021.09.4",
        "scipy >= 1.7.3"
    ],
)
