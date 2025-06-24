#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

import warnings

from pyomo.environ import SolverFactory

from idaes.core.solvers import get_solver

import pytest

from prommis.superstructure.check_superstructure_inputs import check_plant_lifetime_params

### Test plant lifetime parameters
def test_plant_start_type():
    """Test that TypeError is raised if plant_start is not an int."""
    with pytest.raises(TypeError, match="plant_start is not of type int."):
        check_plant_lifetime_params(plant_start="2020", plant_lifetime=10)


def test_plant_lifetime_type():
    """Test that TypeError is raised if plant_lifetime is not an int."""
    with pytest.raises(TypeError, match="plant_lifetime is not of type int."):
        check_plant_lifetime_params(plant_start=2020, plant_lifetime="10")


def test_plant_lifetime_min():
    """Test that ValueError is raised if plant_lifetime is less than 3."""
    with pytest.raises(
        ValueError, match="Plant lifetime must be a minimum of three years."
    ):
        check_plant_lifetime_params(plant_start=2020, plant_lifetime=2)


def test_valid_inputs():
    """Test that no error is raised for valid inputs."""
    check_plant_lifetime_params(plant_start=2020, plant_lifetime=5)