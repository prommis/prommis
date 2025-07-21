#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
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
