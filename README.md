# NorthNet

A Python package which aims to provide tools in understanding chemical reactions at the network level.

Have a look in the scripts folder for examples of how to use NorthNet (there is a walkthrough in `NetworkGenerationNotebook.ipynb`).

## Core functionality

- Reaction-rule based network generation
- Network manipulation and editing
- Conversion of networks into mathematical models
- Experimental reaction data integration with reaction networks
- *Future* data manipulation guided by reaction networks.

## Installation

The programs contained in this repository require Python 3.9.2 and the rdkit, numpy, scipy, networkx, pandas, matplotlib and sklearn libraries (see `environment.yml`).

This software should work on all systems capable of installing the dependencies described above. It has been run successfully on MacOS (10.15) and Windows (Windows 10) machines.

A typical installation time should take less than 1 hour. The software can be installed as follows:

### 1. Clone

Clone the repository from GitHub. If you're having trouble with cloning, just download the zip file and store it on your computer wherever you would like.

### 2. Create a virtual environment

It is possible to create a virtual environment Anaconda, Miniconda (conda) or pip

#### Using conda:

Create a virtual environment with conda:

`conda create --name northnet-env`

Activate the virtual environment:

`conda activate northnet-env`

**Go to Install dependencies.**

#### Using pip:

Create a virtual environment:

`pip install virtualenv`

`virtualenv northnet-env`

Activate the virtual environment:

Mac:

`source northnet-env/bin/activate`

Windows:

`northnet-env\Scripts\activate`

**Go to Install dependencies.**

### Install dependencies

Use `pip` or `conda` to install the following dependencies.

e.g. using conda:
- `conda install -c anaconda scipy`
- `pip install networkx`
- `conda install -c rdkit`

### Install NorthNet

In command line/terminal, navigate to the folder containing the NorthNet code, then type:

conda:
  - `conda develop NorthNet`

  (you may have to run `conda install conda-build` first)

pip:
  - Install from the repository root using pip: `pip install .`,
  - Or in editable mode (so edits are immediately reflected): `pip install -e .`
