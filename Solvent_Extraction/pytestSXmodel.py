import pytest
from pyomo.environ import (
    ConcreteModel,
    value
)
from idaes.core import FlowsheetBlock
from REESXmodel import REESX
from REEAqdistribution import REESolExAqParameters
from REEOgdistribution import REESolExOgParameters

from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.initialization.block_triangularization import BlockTriangularizationInitializer
from idaes.core.initialization import InitializationStatus

from pyomo.util.check_units import assert_units_consistent

solver = get_solver()

class TestSXmodel:
    @pytest.fixture(scope="class")
    def SolEx_frame(self):

        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.prop_a = REESolExAqParameters()
        m.fs.prop_o = REESolExOgParameters()

        m.fs.solex = REESX(number_of_finite_elements=3,
                            aqueous_streams = {"Acidsoln":{"property_package":m.fs.prop_a, "flow_direction":1}},
                            organic_streams = {"Orgacid":{"property_package":m.fs.prop_o, "flow_direction":2}})


        #Aqueous feed fixing
        m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Al"].fix(820)
        m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Ca"].fix(5230)
        m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Fe"].fix(270)
        m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Si"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Sc"].fix(209.31)
        m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Y"].fix(637.74)
        m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["La"].fix(2032.77)
        m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Ce"].fix(4516.13)
        m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Pr"].fix(756.64)
        m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Nd"].fix(2047.85)
        m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Pm"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Sm"].fix(369.1)
        m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Eu"].fix(25.81)
        m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Gd"].fix(174.38)
        m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Tb"].fix(75.28)
        m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Dy"].fix(101.12)
        m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Ho"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Er"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Tm"].fix(41.60)
        m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Yb"].fix(65.65)
        m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Lu"].fix(31.71)
        m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["Th"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].conc_mass_comp["U"].fix(0.01)
        m.fs.solex.Acidsoln_inlet_state[0].flow_vol.fix(4.4)

        #Organic feed fixing
        m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Al"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Ca"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Fe"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Si"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Sc"].fix(321.34)
        m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Y"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["La"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Ce"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Pr"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Nd"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Pm"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Sm"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Eu"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Gd"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Tb"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Dy"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Ho"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Er"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Tm"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Yb"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Lu"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["Th"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].conc_mass_comp["U"].fix(0)

        m.fs.solex.Orgacid_inlet_state[0].flow_vol.fix(62.01)

        return m
    
    @pytest.mark.unit
    def test_dof(self, SolEx_frame):
        m = SolEx_frame
        assert degrees_of_freedom(m) == 0
    
    @pytest.mark.unit
    def test_units(self, SolEx_frame):
        assert_units_consistent(SolEx_frame)
    
    @pytest.mark.component
    def test_solution(self, SolEx_frame):
        m = SolEx_frame
        assert value(m.fs.solex.Acidsoln[0,3].conc_mass_comp["Al"]) == pytest.approx(
            730, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].conc_mass_comp["Ca"]) == pytest.approx(
            4680, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].conc_mass_comp["Fe"]) == pytest.approx(
            250, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].conc_mass_comp["Si"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].conc_mass_comp["Sc"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].conc_mass_comp["Y"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].conc_mass_comp["La"]) == pytest.approx(
            30.84, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].conc_mass_comp["Ce"]) == pytest.approx(
            0.36, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].conc_mass_comp["Pr"]) == pytest.approx(
            0.0312, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].conc_mass_comp["Nd"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].conc_mass_comp["Pm"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].conc_mass_comp["Sm"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].conc_mass_comp["Eu"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].conc_mass_comp["Gd"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].conc_mass_comp["Tb"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].conc_mass_comp["Dy"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].conc_mass_comp["Ho"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].conc_mass_comp["Er"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].conc_mass_comp["Tm"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].conc_mass_comp["Yb"]) == pytest.approx(
            0.47, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].conc_mass_comp["Lu"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].conc_mass_comp["Th"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].conc_mass_comp["U"]) == pytest.approx(
            0, rel=1e-2
        )

        assert value(m.fs.solex.Orgacid[0,1].conc_mass_comp["Al"]) == pytest.approx(
            6.03, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].conc_mass_comp["Ca"]) == pytest.approx(
            39.64, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].conc_mass_comp["Fe"]) == pytest.approx(
            1.19, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].conc_mass_comp["Si"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].conc_mass_comp["Sc"]) == pytest.approx(
            336.25, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].conc_mass_comp["Y"]) == pytest.approx(
            45.41, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].conc_mass_comp["La"]) == pytest.approx(
            142.56, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].conc_mass_comp["Ce"]) == pytest.approx(
            321.59, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].conc_mass_comp["Pr"]) == pytest.approx(
            53.88, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].conc_mass_comp["Nd"]) == pytest.approx(
            145.83, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].conc_mass_comp["Pm"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].conc_mass_comp["Sm"]) == pytest.approx(
            26.28, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].conc_mass_comp["Eu"]) == pytest.approx(
            1.84, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].conc_mass_comp["Gd"]) == pytest.approx(
            12.42, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].conc_mass_comp["Tb"]) == pytest.approx(
            5.36, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].conc_mass_comp["Dy"]) == pytest.approx(
            7.20, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].conc_mass_comp["Ho"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].conc_mass_comp["Er"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].conc_mass_comp["Tm"]) == pytest.approx(
            2.96, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].conc_mass_comp["Yb"]) == pytest.approx(
            4.64, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].conc_mass_comp["Lu"]) == pytest.approx(
            2.26, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].conc_mass_comp["Th"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].conc_mass_comp["U"]) == pytest.approx(
            0, rel=1e-2
        )
    
    @pytest.mark.component
    def test_block_triangularization(self, model):
        initializer = BlockTriangularizationInitializer(constraint_tolerance=2e-5)
        initializer.initialize(model.fs.solex)

        assert initializer.summary[model.fs.solex]["status"] == InitializationStatus.Ok

