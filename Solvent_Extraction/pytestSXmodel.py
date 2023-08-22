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
        m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Al"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Ca"].fix(0.02)
        m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Fe"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Si"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Sc"].fix(0.92)
        m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Y"].fix(2.82)
        m.fs.solex.Acidsoln_inlet_state[0].mass_flow["La"].fix(8.98)
        m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Ce"].fix(19.94)
        m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Pr"].fix(3.34)
        m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Nd"].fix(9.04)
        m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Pm"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Sm"].fix(1.63)
        m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Eu"].fix(0.11)
        m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Gd"].fix(0.77)
        m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Tb"].fix(0.33)
        m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Dy"].fix(0.45)
        m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Ho"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Er"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Tm"].fix(0.18)
        m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Yb"].fix(0.29)
        m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Lu"].fix(0.14)
        m.fs.solex.Acidsoln_inlet_state[0].mass_flow["Th"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].mass_flow["U"].fix(0)

        m.fs.solex.Acidsoln_inlet_state[0].volumetric_flow.fix(4.4)

        #Organic feed fixing
        m.fs.solex.Orgacid_inlet_state[0].mass_flow["Al"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].mass_flow["Ca"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].mass_flow["Fe"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].mass_flow["Si"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].mass_flow["Sc"].fix(19.93)
        m.fs.solex.Orgacid_inlet_state[0].mass_flow["Y"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].mass_flow["La"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].mass_flow["Ce"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].mass_flow["Pr"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].mass_flow["Nd"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].mass_flow["Pm"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].mass_flow["Sm"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].mass_flow["Eu"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].mass_flow["Gd"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].mass_flow["Tb"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].mass_flow["Dy"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].mass_flow["Ho"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].mass_flow["Er"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].mass_flow["Tm"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].mass_flow["Yb"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].mass_flow["Lu"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].mass_flow["Th"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].mass_flow["U"].fix(0)

        m.fs.solex.Orgacid_inlet_state[0].volumetric_flow.fix(62.01)

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
        assert value(m.fs.solex.Acidsoln[0,3].mass_flow["Al"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].mass_flow["Ca"]) == pytest.approx(
            0.02, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].mass_flow["Fe"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].mass_flow["Si"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].mass_flow["Sc"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].mass_flow["Y"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].mass_flow["La"]) == pytest.approx(
            0.14, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].mass_flow["Ce"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].mass_flow["Pr"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].mass_flow["Nd"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].mass_flow["Pm"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].mass_flow["Sm"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].mass_flow["Eu"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].mass_flow["Gd"]) == pytest.approx(
            0.01, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].mass_flow["Tb"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].mass_flow["Dy"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].mass_flow["Ho"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].mass_flow["Er"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].mass_flow["Tm"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].mass_flow["Yb"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].mass_flow["Lu"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].mass_flow["Th"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].mass_flow["U"]) == pytest.approx(
            0, rel=1e-2
        )

        assert value(m.fs.solex.Orgacid[0,1].mass_flow["Al"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].mass_flow["Ca"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].mass_flow["Fe"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].mass_flow["Si"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].mass_flow["Sc"]) == pytest.approx(
            20.85, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].mass_flow["Y"]) == pytest.approx(
            2.82, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].mass_flow["La"]) == pytest.approx(
            8.84, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].mass_flow["Ce"]) == pytest.approx(
            19.94, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].mass_flow["Pr"]) == pytest.approx(
            3.34, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].mass_flow["Nd"]) == pytest.approx(
            9.04, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].mass_flow["Pm"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].mass_flow["Sm"]) == pytest.approx(
            1.63, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].mass_flow["Eu"]) == pytest.approx(
            0.11, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].mass_flow["Gd"]) == pytest.approx(
            0.77, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].mass_flow["Tb"]) == pytest.approx(
            0.33, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].mass_flow["Dy"]) == pytest.approx(
            0.45, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].mass_flow["Ho"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].mass_flow["Er"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].mass_flow["Tm"]) == pytest.approx(
            0.18, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].mass_flow["Yb"]) == pytest.approx(
            0.29, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].mass_flow["Lu"]) == pytest.approx(
            0.14, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].mass_flow["Th"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].mass_flow["U"]) == pytest.approx(
            0, rel=1e-2
        )
