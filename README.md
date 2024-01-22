# prommis
Process Optimization and Modeling for Minerals Sustainability

## Getting started (for PrOMMiS contributors)

Unless otherwise noted, these commands assume that the working directory is the root of the local clone of this repository (i.e. the directory containing this README file).

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
