#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import ConcreteModel, Var

from idaes.core import FlowsheetBlock

import pytest

from prommis.precipitate.precipitate_liquid_properties import AqueousParameter
from prommis.precipitate.precipitate_solids_properties import PrecipitateParameters
from prommis.properties.mixed_acid_properties import (
    MixedAcidParameterBlock,
    get_aliases,
)


@pytest.mark.parametrize("use_mixed_acid", [False, True])
@pytest.mark.unit
def test_build(use_mixed_acid):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    if use_mixed_acid:
        m.fs.prec_solid = PrecipitateParameters(
            aqueous_aliases=get_aliases(include_sulfates=True)
        )
    else:
        m.fs.prec_solid = PrecipitateParameters()

    m.fs.state = m.fs.prec_solid.build_state_block(m.fs.time)

    assert len(m.fs.state) == 1

    assert isinstance(m.fs.state[0].flow_mol_comp, Var)

    m.fs.state[0].flow_mol_comp[:].set_value(0.5)

    m.fs.state.fix_initialization_states()

    for j in m.fs.prec_solid.component_list:
        assert m.fs.state[0].flow_mol_comp[j].fixed


@pytest.mark.unit
def test_react_default():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.prec_solid = PrecipitateParameters()

    assert m.fs.prec_solid.react == {
        "Sc2(C2O4)3(s)": "Sc",
        "Y2(C2O4)3(s)": "Y",
        "La2(C2O4)3(s)": "La",
        "Ce2(C2O4)3(s)": "Ce",
        "Pr2(C2O4)3(s)": "Pr",
        "Nd2(C2O4)3(s)": "Nd",
        "Sm2(C2O4)3(s)": "Sm",
        "Gd2(C2O4)3(s)": "Gd",
        "Dy2(C2O4)3(s)": "Dy",
        "Al2(C2O4)3(s)": "Al",
        "Ca(C2O4)(s)": "Ca",
        "Fe2(C2O4)3(s)": "Fe",
    }


@pytest.mark.unit
def test_react_ionated():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.prec_solid = PrecipitateParameters(
        aqueous_aliases=get_aliases(include_sulfates=True)
    )

    assert m.fs.prec_solid.react == {
        "Sc2(C2O4)3(s)": "Sc_3+",
        "Y2(C2O4)3(s)": "Y_3+",
        "La2(C2O4)3(s)": "La_3+",
        "Ce2(C2O4)3(s)": "Ce_3+",
        "Pr2(C2O4)3(s)": "Pr_3+",
        "Nd2(C2O4)3(s)": "Nd_3+",
        "Sm2(C2O4)3(s)": "Sm_3+",
        "Gd2(C2O4)3(s)": "Gd_3+",
        "Dy2(C2O4)3(s)": "Dy_3+",
        "Al2(C2O4)3(s)": "Al_3+",
        "Ca(C2O4)(s)": "Ca_2+",
        "Fe2(C2O4)3(s)": "Fe_3+",
    }

    # Solid components themselves are not renamed
    assert "Fe2(C2O4)3(s)" in m.fs.prec_solid.component_list


@pytest.mark.parametrize("use_mixed_acid", [False, True])
@pytest.mark.unit
def test_react_species_in_liquid(use_mixed_acid):
    # The species the solids react with must exist in the aqueous
    # property package when both use the same aliases
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    if use_mixed_acid:
        aliases = get_aliases(include_sulfates=True)
        m.fs.prec_solid = PrecipitateParameters(aqueous_aliases=aliases)
        m.fs.prec_sol = AqueousParameter(aqueous_aliases=aliases)
        m.fs.leach_soln = MixedAcidParameterBlock(include_sulfates=True)
        liquid_species = [m.fs.prec_sol, m.fs.leach_soln]
    else:
        m.fs.prec_solid = PrecipitateParameters()
        m.fs.prec_sol = AqueousParameter()
        liquid_species = [m.fs.prec_sol]

    for liquid in liquid_species:
        for j in m.fs.prec_solid.component_list:
            assert m.fs.prec_solid.react[j] in liquid.component_list
