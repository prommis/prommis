####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import ConcreteModel, Var

from idaes.core import FlowsheetBlock

import pytest

from prommis.cmi_custom_precipitator import aqueous_properties as aq_thermo_prop_pack


@pytest.mark.unit
def test_build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Define aqueous species present in system
    aq_comp_list = ["HNO3", "H^+", "OH^-", "NO3^-", "Fe^3+"]

    # define aqueous equilibrium constants
    aq_log_keq_dict = {
        "E1": 1.084,
        "E2": -14,
    }

    # define reaction stoichiometry for aqueous components
    eq_stoich_dict = {
        "E1": {"HNO3": -1, "H^+": 1, "NO3^-": 1, "OH^-": 0, "Fe^3+": 0},
        "E2": {"H^+": 1, "OH^-": 1, "HNO3": 0, "NO3^-": 0, "Fe^3+": 0},
        "E3": {"OH^-": 3, "Fe^3+": 1, "HNO3": 0, "NO3^-": 0, "H^+": 0},
    }

    m.fs.aq_properties = aq_thermo_prop_pack.AqueousParameter(
        aq_comp_list=aq_comp_list,
        eq_rxn_logkeq_dict=aq_log_keq_dict,
        eq_rxn_stoich_dict=eq_stoich_dict,
    )

    m.fs.state = m.fs.aq_properties.build_state_block(m.fs.time)
    assert len(m.fs.state) == 1

    assert isinstance(m.fs.state[0].molality_aq_comp, Var)
    assert isinstance(m.fs.state[0].flow_vol, Var)

    m.fs.state[0].flow_vol.set_value(1)
    for i in m.fs.aq_properties.component_list:
        m.fs.state[0].molality_aq_comp[i].set_value(0.5)

    m.fs.state.fix_initialization_states()

    for i in m.fs.aq_properties.component_list:
        assert m.fs.state[0].molality_aq_comp[i].fixed
