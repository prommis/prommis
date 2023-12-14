from typing import List

import pytest


NOTEBOOK_MARKERS = {
    "tutorial": "Notebook implementing a tutorial",
    "solution": "Notebook for the solution part of a tutorial",
}


def pytest_configure(config: pytest.Config):
    config.pluginmanager.set_blocked("python")

    for name, descr in NOTEBOOK_MARKERS.items():
        config.addinivalue_line("markers", f"{name}: {descr}")


def pytest_collection_modifyitems(items: List[pytest.Item]):
    for item in items:
        if "tutorials" in item.path.parts:
            item.add_marker("tutorial")
            if item.path.stem.endswith("-solution"):
                item.add_marker("solution")
