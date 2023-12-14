import pytest

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
