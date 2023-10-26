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
        m.fs.solex.Acidsoln_inlet_state[0].concentration["Al"].fix(820)
        m.fs.solex.Acidsoln_inlet_state[0].concentration["Ca"].fix(5230)
        m.fs.solex.Acidsoln_inlet_state[0].concentration["Fe"].fix(270)
        m.fs.solex.Acidsoln_inlet_state[0].concentration["Si"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].concentration["Sc"].fix(209.31)
        m.fs.solex.Acidsoln_inlet_state[0].concentration["Y"].fix(637.74)
        m.fs.solex.Acidsoln_inlet_state[0].concentration["La"].fix(2032.77)
        m.fs.solex.Acidsoln_inlet_state[0].concentration["Ce"].fix(4516.13)
        m.fs.solex.Acidsoln_inlet_state[0].concentration["Pr"].fix(756.64)
        m.fs.solex.Acidsoln_inlet_state[0].concentration["Nd"].fix(2047.85)
        m.fs.solex.Acidsoln_inlet_state[0].concentration["Pm"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].concentration["Sm"].fix(369.1)
        m.fs.solex.Acidsoln_inlet_state[0].concentration["Eu"].fix(25.81)
        m.fs.solex.Acidsoln_inlet_state[0].concentration["Gd"].fix(174.38)
        m.fs.solex.Acidsoln_inlet_state[0].concentration["Tb"].fix(75.28)
        m.fs.solex.Acidsoln_inlet_state[0].concentration["Dy"].fix(101.12)
        m.fs.solex.Acidsoln_inlet_state[0].concentration["Ho"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].concentration["Er"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].concentration["Tm"].fix(41.60)
        m.fs.solex.Acidsoln_inlet_state[0].concentration["Yb"].fix(65.65)
        m.fs.solex.Acidsoln_inlet_state[0].concentration["Lu"].fix(31.71)
        m.fs.solex.Acidsoln_inlet_state[0].concentration["Th"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].concentration["U"].fix(0.01)
        m.fs.solex.Acidsoln_inlet_state[0].flow_vol.fix(4.4)

        #Organic feed fixing
        m.fs.solex.Orgacid_inlet_state[0].concentration["Al"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].concentration["Ca"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].concentration["Fe"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].concentration["Si"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].concentration["Sc"].fix(321.34)
        m.fs.solex.Orgacid_inlet_state[0].concentration["Y"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].concentration["La"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].concentration["Ce"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].concentration["Pr"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].concentration["Nd"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].concentration["Pm"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].concentration["Sm"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].concentration["Eu"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].concentration["Gd"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].concentration["Tb"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].concentration["Dy"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].concentration["Ho"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].concentration["Er"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].concentration["Tm"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].concentration["Yb"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].concentration["Lu"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].concentration["Th"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].concentration["U"].fix(0)

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
        assert value(m.fs.solex.Acidsoln[0,3].concentration["Al"]) == pytest.approx(
            730, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].concentration["Ca"]) == pytest.approx(
            4680, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].concentration["Fe"]) == pytest.approx(
            250, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].concentration["Si"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].concentration["Sc"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].concentration["Y"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].concentration["La"]) == pytest.approx(
            30.84, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].concentration["Ce"]) == pytest.approx(
            0.36, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].concentration["Pr"]) == pytest.approx(
            0.0312, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].concentration["Nd"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].concentration["Pm"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].concentration["Sm"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].concentration["Eu"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].concentration["Gd"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].concentration["Tb"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].concentration["Dy"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].concentration["Ho"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].concentration["Er"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].concentration["Tm"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].concentration["Yb"]) == pytest.approx(
            0.47, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].concentration["Lu"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].concentration["Th"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].concentration["U"]) == pytest.approx(
            0, rel=1e-2
        )

        assert value(m.fs.solex.Orgacid[0,1].concentration["Al"]) == pytest.approx(
            6.03, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].concentration["Ca"]) == pytest.approx(
            39.64, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].concentration["Fe"]) == pytest.approx(
            1.19, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].concentration["Si"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].concentration["Sc"]) == pytest.approx(
            336.25, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].concentration["Y"]) == pytest.approx(
            45.41, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].concentration["La"]) == pytest.approx(
            142.56, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].concentration["Ce"]) == pytest.approx(
            321.59, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].concentration["Pr"]) == pytest.approx(
            53.88, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].concentration["Nd"]) == pytest.approx(
            145.83, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].concentration["Pm"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].concentration["Sm"]) == pytest.approx(
            26.28, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].concentration["Eu"]) == pytest.approx(
            1.84, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].concentration["Gd"]) == pytest.approx(
            12.42, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].concentration["Tb"]) == pytest.approx(
            5.36, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].concentration["Dy"]) == pytest.approx(
            7.20, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].concentration["Ho"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].concentration["Er"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].concentration["Tm"]) == pytest.approx(
            2.96, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].concentration["Yb"]) == pytest.approx(
            4.64, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].concentration["Lu"]) == pytest.approx(
            2.26, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].concentration["Th"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].concentration["U"]) == pytest.approx(
            0, rel=1e-2
        )
    
    @pytest.mark.component
    def test_block_triangularization(self, model):
        initializer = BlockTriangularizationInitializer(constraint_tolerance=2e-5)
        initializer.initialize(model.fs.solex)

        assert initializer.summary[model.fs.solex]["status"] == InitializationStatus.Ok

