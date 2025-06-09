#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import ConcreteModel, Var

from idaes.core import FlowsheetBlock

import pytest

from prommis.cmi_precipitator import (
    precipitate_properties as precipitate_thermo_prop_pack,
)


@pytest.mark.unit
def test_build():
    # Define precipitate species present in system
    precipitate_comp_list = ["FeOH3"]

    # define equilibrium constant for precipitation/dissolution reactions
    precipitate_log_keq_dict = {
        "E3": -33.498,
    }

    # define reaction stoichiometry for precipitates
    precipitate_stoich_dict = {
        "E1": {"FeOH3": 0},
        "E2": {"FeOH3": 0},
        "E3": {"FeOH3": -1},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.precipitate_properties = precipitate_thermo_prop_pack.PrecipitateParameter(
        precipitate_comp_list=precipitate_comp_list,
        logkeq_dict=precipitate_log_keq_dict,
        stoich_dict=precipitate_stoich_dict,
    )

    m.fs.state = m.fs.precipitate_properties.build_state_block(m.fs.time)

    assert len(m.fs.state) == 1

    assert isinstance(m.fs.state[0].moles_precipitate_comp, Var)

    for i in m.fs.precipitate_properties.component_list:
        m.fs.state[0].moles_precipitate_comp[i].set_value(0.5)

    m.fs.state.fix_initialization_states()

    for i in m.fs.precipitate_properties.component_list:
        assert m.fs.state[0].moles_precipitate_comp[i].fixed
