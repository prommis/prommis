#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from typing import List

import pytest


def pytest_addoption(parser: pytest.Parser):
    parser.addoption(
        "--xfail-known-issues",
        action="store_true",
        default=False,
    )


def pytest_configure(config: pytest.Config):
    config.addinivalue_line(
        "markers", "known_issue(issue_number): Test affected by known issue(s)"
    )


def _get_known_issues(
    item: pytest.Item, repository: str = "prommis/prommis"
) -> List[str]:
    links = []
    for marker in item.iter_markers(name="known_issue"):
        for number in marker.args:
            links.append(f"{repository}#{number}")
    return links


def _xfail_known_issues(item: pytest.Item) -> None:
    known_issues = _get_known_issues(item)
    if known_issues:
        issues_to_display = ", ".join(known_issues)
        marker = pytest.mark.xfail(
            strict=False,
            reason=f"This test is known to be affected by issue(s): {issues_to_display}",
        )
        item.add_marker(marker)


def pytest_collection_modifyitems(config: pytest.Config, items: List[pytest.Item]):
    _handle_known_issues = lambda *args: None
    if config.getoption("--xfail-known-issues"):
        _handle_known_issues = _xfail_known_issues

    for item in items:
        _handle_known_issues(item)
