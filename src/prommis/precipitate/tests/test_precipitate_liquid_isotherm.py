#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
import io
import logging

from pyomo.environ import ConcreteModel, Constraint, Var, value
from pyomo.common.log import LoggingIntercept

from idaes.core import FlowsheetBlock

import pytest

from prommis.precipitate.precipitate_liquid_properties_isotherm import AqueousParameter
from prommis.properties.mixed_acid_properties import (
    MixedAcidParameterBlock,
    get_aliases,
)


@pytest.mark.unit
def test_deprecation_warning_on_parameter_block():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    stream = io.StringIO()
    with LoggingIntercept(stream, "pyomo", logging.WARNING):
        m.fs.prec_sol = AqueousParameter()

    output = stream.getvalue()
    assert "DEPRECATED" in output
    assert "PrOMMiS 1.1" in output


@pytest.mark.unit
def test_deprecation_warning_on_state_block():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    stream = io.StringIO()
    with LoggingIntercept(None, "pyomo", 0):
        m.fs.prec_sol = AqueousParameter()

    stream = io.StringIO()
    with LoggingIntercept(stream, "pyomo", logging.WARNING):
        m.fs.state = m.fs.prec_sol.build_state_block(m.fs.time)

    output = stream.getvalue()
    assert "DEPRECATED" in output
    assert "PrOMMiS 1.1" in output


@pytest.mark.parametrize("use_mixed_acid", [False, True])
@pytest.mark.unit
def test_build(use_mixed_acid):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    if use_mixed_acid:
        m.fs.prec_sol = AqueousParameter(
            aqueous_aliases=get_aliases(include_sulfates=True, include_oxalates=True)
        )
    else:
        m.fs.prec_sol = AqueousParameter()

    m.fs.state = m.fs.prec_sol.build_state_block(m.fs.time)

    assert len(m.fs.state) == 1

    assert isinstance(m.fs.state[0].flow_vol, Var)
    assert isinstance(m.fs.state[0].conc_mass_comp, Var)
    assert isinstance(m.fs.state[0].flow_mol_comp, Var)
    assert isinstance(m.fs.state[0].conc_mol_comp, Var)
    assert isinstance(m.fs.state[0].temperature, Var)
    assert isinstance(m.fs.state[0].pressure, Var)

    assert isinstance(m.fs.state[0].flow_mol_constraint, Constraint)
    assert isinstance(m.fs.state[0].molar_concentration_constraint, Constraint)

    m.fs.state[0].flow_vol.set_value(10)
    for i in m.fs.prec_sol.dissolved_elements:
        m.fs.state[0].conc_mass_comp[i].set_value(0.1)
        m.fs.state[0].flow_mol_comp[i].set_value(0.5)

    m.fs.state.fix_initialization_states()

    assert m.fs.state[0].flow_vol.fixed
    for j in m.fs.prec_sol.dissolved_elements:
        assert m.fs.state[0].conc_mass_comp[j].fixed
        assert not m.fs.state[0].flow_mol_comp[j].fixed

        assert m.fs.state[0].flow_mol_constraint[j].active
        assert m.fs.state[0].molar_concentration_constraint[j].active


@pytest.mark.unit
def test_default_names():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.prec_sol = AqueousParameter()

    assert list(m.fs.prec_sol.component_list) == [
        "H2O",
        "Al",
        "Ca",
        "Fe",
        "H",
        "Cl",
        "HSO4",
        "SO4",
        "H2C2O4",
        "Sc",
        "Y",
        "La",
        "Ce",
        "Pr",
        "Nd",
        "Sm",
        "Gd",
        "Dy",
    ]
    assert set(m.fs.prec_sol.dissolved_elements) == set(m.fs.prec_sol.component_list)


@pytest.mark.unit
def test_ionated_names():
    aliases = get_aliases(include_sulfates=True, include_oxalates=True)

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.default = AqueousParameter()
    m.fs.ionated = AqueousParameter(aqueous_aliases=aliases)

    assert len(m.fs.ionated.component_list) == len(m.fs.default.component_list)
    assert set(m.fs.ionated.dissolved_elements) == {
        aliases[j] for j in m.fs.default.dissolved_elements
    }
    assert "Fe_3+" in m.fs.ionated.component_list
    assert "Fe" not in m.fs.ionated.component_list
    assert "H2C2O4" in m.fs.ionated.component_list

    # Parameters are unchanged apart from the names they are indexed by
    for j in m.fs.default.component_list:
        assert value(m.fs.ionated.mw[aliases[j]]) == pytest.approx(
            value(m.fs.default.mw[j]), rel=1e-12
        )


@pytest.mark.unit
def test_ionated_names_match_mixed_acid():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.prec_sol = AqueousParameter(
        aqueous_aliases=get_aliases(include_sulfates=True, include_oxalates=True)
    )
    m.fs.leach_soln = MixedAcidParameterBlock(
        include_sulfates=True, include_oxalates=True
    )

    assert set(m.fs.prec_sol.component_list) <= set(m.fs.leach_soln.component_list)
