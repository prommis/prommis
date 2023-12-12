from pytest import Config


IDAES_MARKERS = {
    "unit": "Quick tests that do not require a solver, must run in < 2 s",
    "component": "Quick tests that may require a solver",
    "integration": "Long duration tests",
}


def pytest_configure(config: Config):
    for spec, descr in IDAES_MARKERS.items():
        config.addinivalue_line("markers", f"{spec}: {descr}")
