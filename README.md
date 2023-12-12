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

#### Formatting code

```sh
black .
```

#### Running linter (Pylint)

```sh
pylint prommis
```

#### Running tests

```sh
pytest                          # run the complete test suite
pytest -k test_my_flowsheet.py  # run only test defined in the file named test_my_flowsheet.py
```
