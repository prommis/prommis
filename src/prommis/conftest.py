#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
import pytest

try:
    # NOTE attempting to import watertap here ensures that watertap has been imported
    # (if available) before any test file is run
    # this is necessary since watertap installs its own Ipopt wrapper on import
    # which can affect solves even for models that don't use WaterTAP
    # see prommis/prommis#52 for one such example
    import watertap  # pylint: disable=unused-import
except ModuleNotFoundError:
    pass


IDAES_MARKERS = {
    "build": "Test of model build methods",
    "unit": "Quick tests that do not require a solver, must run in < 2 s",
    "component": "Quick tests that may require a solver",
    "integration": "Long duration tests",
    "solver": "Test requires a solver",
}


def pytest_configure(config: pytest.Config):
    for spec, descr in IDAES_MARKERS.items():
        config.addinivalue_line("markers", f"{spec}: {descr}")
