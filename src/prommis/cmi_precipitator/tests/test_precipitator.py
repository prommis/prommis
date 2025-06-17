#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import ConcreteModel, Constraint, assert_optimal_termination, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_total_constraints,
    number_unused_variables,
    number_variables,
)
from idaes.core.util.scaling import (
    constraint_scaling_transform,
    get_scaling_factor,
    set_scaling_factor,
)

import pytest

from prommis.cmi_precipitator import aqueous_properties as aqueous_thermo_prop_pack
from prommis.cmi_precipitator import (
    precipitate_properties as precipitate_thermo_prop_pack,
)
from prommis.cmi_precipitator.opt_based_precipitator import Precipitator

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
class TestPrec(object):
    @pytest.fixture(scope="class")
    def prec(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        # Define aqueous species present in system
        aqueous_comp_list = ["HNO3", "H^+", "OH^-", "NO3^-", "Fe^3+"]
        # Define precipitate species present in system
        precipitate_comp_list = ["FeOH3"]

        # define aqueous equilibrium constants
        aqueous_log_keq_dict = {
            "E1": 1.084,
            "E2": -14,
        }

        # define equilibrium constant for precipitation/dissolution reactions
        precipitate_log_keq_dict = {
            "E3": -33.498,
        }

        # define reaction stoichiometry for aqueous components
        aqueous_stoich_dict = {
            "E1": {"HNO3": -1, "H^+": 1, "NO3^-": 1, "OH^-": 0, "Fe^3+": 0},
            "E2": {"H^+": 1, "OH^-": 1, "HNO3": 0, "NO3^-": 0, "Fe^3+": 0},
            "E3": {"OH^-": 3, "Fe^3+": 1, "HNO3": 0, "NO3^-": 0, "H^+": 0},
        }

        # define reaction stoichiometry for precipitates
        precipitate_stoich_dict = {
            "E1": {"FeOH3": 0},
            "E2": {"FeOH3": 0},
            "E3": {"FeOH3": -1},
        }

        m.fs.aqueous_properties = aqueous_thermo_prop_pack.AqueousParameter(
            aqueous_comp_list=aqueous_comp_list,
            logkeq_dict=aqueous_log_keq_dict,
            stoich_dict=aqueous_stoich_dict,
        )
        m.fs.precipitate_properties = precipitate_thermo_prop_pack.PrecipitateParameter(
            precipitate_comp_list=precipitate_comp_list,
            logkeq_dict=precipitate_log_keq_dict,
            stoich_dict=precipitate_stoich_dict,
        )

        m.fs.unit = Precipitator(
            property_package_aqueous=m.fs.aqueous_properties,
            property_package_precipitate=m.fs.precipitate_properties,
        )

        # initial concentrations
        pH = 7
        Hp_init = 10 ** (-pH)
        OHm_init = 10 ** (-14 + pH)

        Feppp_init = 1.05e-1
        NO3m_init = 1.8e-1
        HNO3m_init = 1e-10
        FeOH3_init = 1e-10
        Fe2O3_init = 1e-10

        m.fs.unit.aqueous_inlet.flow_vol[0].fix(1)
        m.fs.unit.aqueous_inlet.molality_aqueous_comp[0, "HNO3"].fix(HNO3m_init)
        m.fs.unit.aqueous_inlet.molality_aqueous_comp[0, "H^+"].fix(Hp_init)
        m.fs.unit.aqueous_inlet.molality_aqueous_comp[0, "OH^-"].fix(OHm_init)
        m.fs.unit.aqueous_inlet.molality_aqueous_comp[0, "NO3^-"].fix(NO3m_init)
        m.fs.unit.aqueous_inlet.molality_aqueous_comp[0, "Fe^3+"].fix(Feppp_init)
        m.fs.unit.precipitate_inlet.moles_precipitate_comp[0, "FeOH3"].fix(FeOH3_init)

        return m

    @pytest.mark.unit
    def test_config(self, prec):

        assert len(prec.fs.unit.config()) == 9

        assert not prec.fs.unit.config.dynamic
        assert not prec.fs.unit.config.has_holdup
        assert prec.fs.unit.config.has_equilibrium_reactions
        assert not prec.fs.unit.config.has_phase_equilibrium
        assert not prec.fs.unit.config.has_heat_of_reaction
        assert (
            prec.fs.unit.config.property_package_aqueous is prec.fs.aqueous_properties
        )
        assert (
            prec.fs.unit.config.property_package_precipitate
            is prec.fs.precipitate_properties
        )

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, prec):
        assert hasattr(prec.fs.unit, "aqueous_inlet")
        assert len(prec.fs.unit.aqueous_inlet.vars) == 2
        assert hasattr(prec.fs.unit.aqueous_inlet, "flow_vol")
        assert hasattr(prec.fs.unit.aqueous_inlet, "molality_aqueous_comp")

        assert hasattr(prec.fs.unit, "aqueous_outlet")
        assert len(prec.fs.unit.aqueous_outlet.vars) == 2
        assert hasattr(prec.fs.unit.aqueous_outlet, "flow_vol")
        assert hasattr(prec.fs.unit.aqueous_outlet, "molality_aqueous_comp")

        assert hasattr(prec.fs.unit, "precipitate_inlet")
        assert len(prec.fs.unit.precipitate_inlet.vars) == 1
        assert hasattr(prec.fs.unit.precipitate_inlet, "moles_precipitate_comp")

        assert hasattr(prec.fs.unit, "precipitate_outlet")
        assert len(prec.fs.unit.precipitate_outlet.vars) == 1
        assert hasattr(prec.fs.unit.precipitate_outlet, "moles_precipitate_comp")

        assert hasattr(prec.fs.unit, "log_q_precipitate_equilibrium_rxn_eqns")
        assert hasattr(prec.fs.unit, "precip_rxns_log_cons")
        assert hasattr(prec.fs.unit, "aqueous_mole_balance_eqns")
        assert hasattr(prec.fs.unit, "precipitate_mole_balance_eqns")
        assert hasattr(prec.fs.unit, "min_logs")
        assert hasattr(prec.fs.unit, "vol_balance")

        assert number_variables(prec.fs.unit) == 18
        assert number_total_constraints(prec.fs.unit) == 11
        assert number_unused_variables(prec.fs.unit) == 0

    @pytest.mark.component
    def test_units(self, prec):
        assert_units_consistent(prec.fs.unit)

    @pytest.mark.component
    def test_dof(self, prec):
        assert degrees_of_freedom(prec) == 1

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, prec):
        # scale model
        set_scaling_factor(prec.fs.unit.log_k_equilibrium_rxn_eqns[0.0, "E1"], 1e-4)
        set_scaling_factor(prec.fs.unit.log_k_equilibrium_rxn_eqns[0.0, "E2"], 1e-10)
        set_scaling_factor(
            prec.fs.unit.log_q_precipitate_equilibrium_rxn_eqns[0.0, "E3"], 1e-10
        )
        for con in prec.fs.component_data_objects(Constraint, active=True):
            scaling_factor = get_scaling_factor(con, default=None)
            if scaling_factor is not None:
                constraint_scaling_transform(con, scaling_factor)

        solver = get_solver()
        solver.options["nlp_scaling_method"] = "user-scaling"
        results = solver.solve(prec)

        assert_optimal_termination(results)

        dt = DiagnosticsToolbox(prec)
        dt.assert_no_numerical_warnings()

        # fixing final concentration of precipitate to optimal value to avoid non-zero dof warnings
        optimal_feoh3 = value(
            prec.fs.unit.precipitate_outlet.moles_precipitate_comp[0, "FeOH3"]
        )
        prec.fs.unit.precipitate_outlet.moles_precipitate_comp[0, "FeOH3"].fix(
            optimal_feoh3
        )
        dt.assert_no_structural_warnings()

    @pytest.mark.solver
    @pytest.mark.component
    def test_solution(self, prec):
        # aqeous species final concentrations
        assert pytest.approx(0.104766, abs=1e-5) == value(
            prec.fs.unit.aqueous_outlet.molality_aqueous_comp[0, "Fe^3+"]
        )
        assert pytest.approx(1.025171e-5, abs=1e-5) == value(
            prec.fs.unit.aqueous_outlet.molality_aqueous_comp[0, "HNO3"]
        )
        assert pytest.approx(0.000691, abs=1e-5) == value(
            prec.fs.unit.aqueous_outlet.molality_aqueous_comp[0, "H^+"]
        )
        assert pytest.approx(0.1799897, abs=1e-5) == value(
            prec.fs.unit.aqueous_outlet.molality_aqueous_comp[0, "NO3^-"]
        )
        assert pytest.approx(1.446944e-11, abs=1e-5) == value(
            prec.fs.unit.aqueous_outlet.molality_aqueous_comp[0, "OH^-"]
        )

        # precipitate species final amounts
        assert pytest.approx(0.00023372, abs=1e-5) == value(
            prec.fs.unit.precipitate_outlet.moles_precipitate_comp[0, "FeOH3"]
        )
