import pytest

import prommis


def test_version():
    assert bool(prommis.__version__)
