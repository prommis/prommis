# Getting started

This page gets you started with the PrOMMiS Python package.

It depends on the [IDAES](https://idaes-pse.readthedocs.io/en/stable/) and [WaterTAP](https://watertap.readthedocs.io/en/stable/) Python packages, which means that these will be automatically installed as part of its installation.

## Install PrOMMiS

PrOMMiS is distributed as a Python package, so can be installed through Python's standard tool, `pip`. These instructions assume you can run commands from a Mac OSX, Linux, or Windows Powershell terminal.

### Create Python environment (optional)

Before installing the PrOMMiS software, we recommend you set up an "environment" that will let you install Python packages without affecting your default setup.
Instructions are given below for doing this using Conda and Miniforge, but there are other methods such as [`pyenv`](https://github.com/pyenv/pyenv) and [`virtualenv`](https://virtualenv.pypa.io/en/latest/) that have similar effects.

### Create Conda / Miniforge environment

Conda is the name for the packaging environment created by the Anaconda corporation.
Miniforge is software that installs Python packages, and sets up
Python environments, without commercial and licensing concerns associated
with the Anaconda software and package libraries.

Long story short, install miniforge from these links:

- [Home page](https://github.com/conda-forge/miniforge)
- [Releases](https://github.com/conda-forge/miniforge/releases)

This will enable the `conda` command you can use to install a new Python
environment:

```
conda create -n prommis python=3.11 -y
```

Before you do any work you need to "activate" the environment:

```
conda activate prommis
```

### Install PrOMMiS

We will now install the PrOMMiS Python package into this environment.

```
pip install prommis
```

This will take a couple of minutes.
Once it is done, you can continue to try an example.

### Install IDAES “extensions”

After installing PrOMMiS, you will need to install binary packages for the solvers.
These “extensions”, as they are called, can be installed, after installing IDAES, using this command:

```
idaes get-extensions
```

:::{note}
Depending on your operating system, additional steps might be needed. For more information, refer to the [IDAES installation guide](https://idaes-pse.readthedocs.io/en/stable/tutorials/getting_started/index.html).
:::

## Try an example

This section fetches and runs a Jupyter Notebook containing a PrOMMiS flowsheet.

We start by creating a working directory, then download the notebook from GitHub,
then -- after making sure the Jupyter Notebook code is installed -- open that
notebook in Jupyter.

### Create subdirectory to work in

```sh
mkdir prommis-work
cd prommis-work
```

### Download the UKy flowsheet Jupyter Notebook

::::{tab-set}

:::{tab-item} Console (curl)
Run the following commands using the [curl](https://curl.se/) command-line tool.

```sh
curl -o uky_flowsheet-solution.ipynb https://raw.githubusercontent.com/prommis/prommis/refs/heads/main/docs/tutorials/uky_flowsheet-solution.ipynb
curl -o uky_flowsheet.png https://raw.githubusercontent.com/prommis/prommis/refs/heads/main/docs/tutorials/uky_flowsheet.png
```

:::

:::{tab-item} Web browser

Navigate to [the Jupyter notebook](https://github.com/prommis/prommis/blob/main/docs/tutorials/uky_flowsheet-solution.ipynb) and hit Download ({octicon}`download`) link in upper right of the file preview box.

Do the same with [the flowsheet diagram](https://github.com/prommis/prommis/blob/main/docs/tutorials/uky_flowsheet.png)

:::

::::

### Install Jupyter Lab

In order to run the Jupyter Notebook you'll need the current web interface for Jupyter Notebooks,
called Jupyter Lab, to be installed.

```sh
pip install jupyterlab
```

### Open the notebook

From here, all you need to do is run Jupyter Lab using the notebook you downloaded as input:

```sh
 jupyter lab uky_flowsheet-solution.ipynb
 ```

 This should open a new browser window or tab which shows the Jupyter Lab interface and the UKy notebook.

 :::{tip}
 To hide the sidebar, and view just the notebook, hit "Control-B".
 :::

### Run the notebook

For those unfamiliar with Jupyter Notebooks, the most common operation is to press Ctrl + Enter to execute the code in the current cell and advance to the next one. In this way you can execute a notebook step by step, pausing as needed to read and examine output. For more information see the [JupyterLab documentation](https://jupyterlab.readthedocs.io/en/latest/).

## Next steps

To use PrOMMiS libraries in your own Python code, you will need to create Python modules that create and build flowsheets, such as the module for the [UKy flowsheet](https://github.com/prommis/prommis/blob/main/src/prommis/uky/uky_flowsheet.py).

To understand how this and other examples work and start to create your own, please refer to documentation of the [IDAES](https://idaes-pse.readthedocs.io/en/stable/) process systems engineering toolkit software.

## For PrOMMiS developers

Unless otherwise noted, these commands assume that the working directory is the root of the local clone of this repository (i.e., the directory containing the [README file](../README.md)).

### Installation

```sh
conda create --name prommis-dev --yes python=3.11
conda activate prommis-dev
git clone https://github.com/prommis/prommis && cd prommis
pip install -r requirements-dev.txt
```

### Before committing

Before running any of these commands, ensure the `prommis-dev` Conda environment has been activated:

```sh
conda activate prommis-dev
```

#### Sort import statements

```sh
isort src/prommis
```

#### Formatting code

```sh
black .
```

#### Running linter (Pylint)

```sh
pylint prommis
```

#### Running spell checker (Typos)

```sh
typos
```

Note: if the `typos` executable is not found, it can be installed by running `conda install --yes -c conda-forge typos` after activating the `prommis-dev` Conda environment.

#### Running tests

```sh
pytest                          # run the complete test suite
pytest -k test_my_flowsheet.py  # run only test defined in the file named test_my_flowsheet.py
```

#### Building documentation

From the `docs/` subdirectory:

```sh
jupyter-book build .
```

#### Testing (executing) notebooks

From the `docs/` subdirectory:

```sh
pytest --nbmake -m "solution" .
```
