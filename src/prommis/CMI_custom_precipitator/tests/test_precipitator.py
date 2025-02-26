#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
from pyomo.environ import ConcreteModel, assert_optimal_termination, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    InitializationStatus,
)
from idaes.core.solvers import get_solver
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_total_constraints,
    number_unused_variables,
    number_variables,
)
from idaes.core.util.scaling import set_scaling_factor, unscaled_variables_generator

import pytest

from prommis.CMI_custom_precipitator import AqueousProperties as aq_thermo_prop_pack
from prommis.CMI_custom_precipitator import (
    PrecipitateProperties as precip_thermo_prop_pack,
)
from prommis.CMI_custom_precipitator.opt_based_precipitator import Precipitator


class TestPrec(object):
    @pytest.fixture(scope="class")
    def prec(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        # Define aqueous species present in system
        aq_comp_list = ["HNO3", "H^+", "OH^-", "NO3^-", "Fe^3+"]
        # Define precipitate species present in system
        precip_comp_list = ["FeOH3"]

        # define aqueous equilibrium constants
        aq_log_keq_dict = {
            "E1": 1.084,
            "E2": -14,
        }

        # define equilibrium constant for precipitation/dissolution reactions
        precip_log_keq_dict = {
            "E3": -33.498,
        }

        # define reaction stoichiometry for aqueous components
        eq_stoich_dict = {
            "E1": {"HNO3": -1, "H^+": 1, "NO3^-": 1, "OH^-": 0, "Fe^3+": 0},
            "E2": {"H^+": 1, "OH^-": 1, "HNO3": 0, "NO3^-": 0, "Fe^3+": 0},
            "E3": {"OH^-": 3, "Fe^3+": 1, "HNO3": 0, "NO3^-": 0, "H^+": 0},
        }

        # define reaction stoichiometry for precipitates
        precip_eq_stoich_dict = {
            "E1": {"FeOH3": 0},
            "E2": {"FeOH3": 0},
            "E3": {"FeOH3": -1},
        }

        m.fs.aq_properties = aq_thermo_prop_pack.AqueousParameter(
            aq_comp_list=aq_comp_list,
            eq_rxn_logkeq_dict=aq_log_keq_dict,
            eq_rxn_stoich_dict=eq_stoich_dict,
        )
        m.fs.precip_properties = precip_thermo_prop_pack.PrecipitateParameter(
            precip_comp_list=precip_comp_list,
            precip_eq_rxn_logkeq_dict=precip_log_keq_dict,
            precip_eq_rxn_stoich_dict=precip_eq_stoich_dict,
        )

        m.fs.unit = Precipitator(
            property_package_aqueous=m.fs.aq_properties,
            property_package_precipitate=m.fs.precip_properties,
        )

        # initial concentrations
        pH = 7
        Hp_init = 10 ** (-pH)
        OHm_init = 10 ** (-14 + pH)

        Feppp_init = 1.05e-1
        NO3m_init = 1.8e-1
        HNO3m_init = 1e-20
        FeOH3_init = 1e-20
        Fe2O3_init = 1e-20

        m.fs.unit.aqueous_inlet.flow_vol[0].fix(1)
        m.fs.unit.aqueous_inlet.molality_aq_comp[0, "HNO3"].fix(HNO3m_init)
        m.fs.unit.aqueous_inlet.molality_aq_comp[0, "H^+"].fix(Hp_init)
        m.fs.unit.aqueous_inlet.molality_aq_comp[0, "OH^-"].fix(OHm_init)
        m.fs.unit.aqueous_inlet.molality_aq_comp[0, "NO3^-"].fix(NO3m_init)
        m.fs.unit.aqueous_inlet.molality_aq_comp[0, "Fe^3+"].fix(Feppp_init)
        m.fs.unit.precipitate_inlet.moles_precip_comp[0, "FeOH3"].fix(FeOH3_init)

        return m

    @pytest.mark.unit
    def test_config(self, prec):

        assert len(prec.fs.unit.config()) == 9

        assert not prec.fs.unit.config.dynamic
        assert not prec.fs.unit.config.has_holdup
        assert prec.fs.unit.config.has_equilibrium_reactions
        assert not prec.fs.unit.config.has_phase_equilibrium
        assert not prec.fs.unit.config.has_heat_of_reaction
        assert prec.fs.unit.config.property_package_aqueous is prec.fs.aq_properties
        assert (
            prec.fs.unit.config.property_package_precipitate
            is prec.fs.precip_properties
        )

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, prec):
        assert hasattr(prec.fs.unit, "aqueous_inlet")
        assert len(prec.fs.unit.aqueous_inlet.vars) == 2
        assert hasattr(prec.fs.unit.aqueous_inlet, "flow_vol")
        assert hasattr(prec.fs.unit.aqueous_inlet, "molality_aq_comp")

        assert hasattr(prec.fs.unit, "aqueous_outlet")
        assert len(prec.fs.unit.aqueous_outlet.vars) == 2
        assert hasattr(prec.fs.unit.aqueous_outlet, "flow_vol")
        assert hasattr(prec.fs.unit.aqueous_outlet, "molality_aq_comp")

        assert hasattr(prec.fs.unit, "precipitate_inlet")
        assert len(prec.fs.unit.precipitate_inlet.vars) == 1
        assert hasattr(prec.fs.unit.precipitate_inlet, "moles_precip_comp")

        assert hasattr(prec.fs.unit, "precipitate_outlet")
        assert len(prec.fs.unit.precipitate_outlet.vars) == 1
        assert hasattr(prec.fs.unit.precipitate_outlet, "moles_precip_comp")

        assert hasattr(prec.fs.unit, "log_q_precip_equil_rxn_eqns")
        assert hasattr(prec.fs.unit, "precip_rxns_log_cons")
        assert hasattr(prec.fs.unit, "aq_mole_balance_eqns")
        assert hasattr(prec.fs.unit, "precip_mole_balance_eqns")
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
    @pytest.mark.component
    def test_solve(self, prec):
        # scale model
        set_scaling_factor(
            prec.fs.unit.cv_aqueous.properties_out[0.0].molality_aq_comp["H^+"], 1e4
        )
        set_scaling_factor(
            prec.fs.unit.cv_aqueous.properties_out[0.0].molality_aq_comp["OH^-"], 1e11
        )
        set_scaling_factor(
            prec.fs.unit.cv_aqueous.properties_out[0.0].molality_aq_comp["NO3^-"], 1e2
        )
        set_scaling_factor(
            prec.fs.unit.cv_aqueous.properties_out[0.0].molality_aq_comp["HNO3"], 1e5
        )
        set_scaling_factor(
            prec.fs.unit.cv_aqueous.properties_out[0.0].molality_aq_comp["Fe^3+"], 1e2
        )

        # set scaling for precipitate final amount
        set_scaling_factor(
            prec.fs.unit.cv_precipitate.properties_out[0.0].moles_precip_comp["FeOH3"],
            1e4,
        )

        results = prec.fs.unit.solve_unit()
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.component
    def test_solution(self, prec):
        # aqeous species final concentrations
        assert pytest.approx(0.104766, abs=1e-5) == value(
            prec.fs.unit.aqueous_outlet.molality_aq_comp[0, "Fe^3+"]
        )
        assert pytest.approx(1.025171e-5, abs=1e-5) == value(
            prec.fs.unit.aqueous_outlet.molality_aq_comp[0, "HNO3"]
        )
        assert pytest.approx(0.000691, abs=1e-5) == value(
            prec.fs.unit.aqueous_outlet.molality_aq_comp[0, "H^+"]
        )
        assert pytest.approx(0.1799897, abs=1e-5) == value(
            prec.fs.unit.aqueous_outlet.molality_aq_comp[0, "NO3^-"]
        )
        assert pytest.approx(1.446944e-11, abs=1e-5) == value(
            prec.fs.unit.aqueous_outlet.molality_aq_comp[0, "OH^-"]
        )

        # precipitate species final amounts
        assert pytest.approx(0.00023372, abs=1e-5) == value(
            prec.fs.unit.precipitate_outlet.moles_precip_comp[0, "FeOH3"]
        )
