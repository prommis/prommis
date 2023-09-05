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
        m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Al"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Ca"].fix(0.02)
        m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Fe"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Si"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Sc"].fix(0.92)
        m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Y"].fix(2.82)
        m.fs.solex.Acidsoln_inlet_state[0].flow_mass["La"].fix(8.98)
        m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Ce"].fix(19.94)
        m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Pr"].fix(3.34)
        m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Nd"].fix(9.04)
        m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Pm"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Sm"].fix(1.63)
        m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Eu"].fix(0.11)
        m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Gd"].fix(0.77)
        m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Tb"].fix(0.33)
        m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Dy"].fix(0.45)
        m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Ho"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Er"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Tm"].fix(0.18)
        m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Yb"].fix(0.29)
        m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Lu"].fix(0.14)
        m.fs.solex.Acidsoln_inlet_state[0].flow_mass["Th"].fix(0)
        m.fs.solex.Acidsoln_inlet_state[0].flow_mass["U"].fix(0)

        m.fs.solex.Acidsoln_inlet_state[0].flow_vol.fix(4.4)

        #Organic feed fixing
        m.fs.solex.Orgacid_inlet_state[0].flow_mass["Al"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].flow_mass["Ca"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].flow_mass["Fe"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].flow_mass["Si"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].flow_mass["Sc"].fix(19.93)
        m.fs.solex.Orgacid_inlet_state[0].flow_mass["Y"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].flow_mass["La"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].flow_mass["Ce"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].flow_mass["Pr"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].flow_mass["Nd"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].flow_mass["Pm"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].flow_mass["Sm"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].flow_mass["Eu"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].flow_mass["Gd"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].flow_mass["Tb"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].flow_mass["Dy"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].flow_mass["Ho"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].flow_mass["Er"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].flow_mass["Tm"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].flow_mass["Yb"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].flow_mass["Lu"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].flow_mass["Th"].fix(0)
        m.fs.solex.Orgacid_inlet_state[0].flow_mass["U"].fix(0)

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
        assert value(m.fs.solex.Acidsoln[0,3].flow_mass["Al"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].flow_mass["Ca"]) == pytest.approx(
            0.02, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].flow_mass["Fe"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].flow_mass["Si"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].flow_mass["Sc"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].flow_mass["Y"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].flow_mass["La"]) == pytest.approx(
            0.14, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].flow_mass["Ce"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].flow_mass["Pr"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].flow_mass["Nd"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].flow_mass["Pm"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].flow_mass["Sm"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].flow_mass["Eu"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].flow_mass["Gd"]) == pytest.approx(
            0.01, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].flow_mass["Tb"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].flow_mass["Dy"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].flow_mass["Ho"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].flow_mass["Er"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].flow_mass["Tm"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].flow_mass["Yb"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].flow_mass["Lu"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].flow_mass["Th"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Acidsoln[0,3].flow_mass["U"]) == pytest.approx(
            0, rel=1e-2
        )

        assert value(m.fs.solex.Orgacid[0,1].flow_mass["Al"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].flow_mass["Ca"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].flow_mass["Fe"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].flow_mass["Si"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].flow_mass["Sc"]) == pytest.approx(
            20.85, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].flow_mass["Y"]) == pytest.approx(
            2.82, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].flow_mass["La"]) == pytest.approx(
            8.84, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].flow_mass["Ce"]) == pytest.approx(
            19.94, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].flow_mass["Pr"]) == pytest.approx(
            3.34, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].flow_mass["Nd"]) == pytest.approx(
            9.04, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].flow_mass["Pm"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].flow_mass["Sm"]) == pytest.approx(
            1.63, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].flow_mass["Eu"]) == pytest.approx(
            0.11, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].flow_mass["Gd"]) == pytest.approx(
            0.77, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].flow_mass["Tb"]) == pytest.approx(
            0.33, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].flow_mass["Dy"]) == pytest.approx(
            0.45, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].flow_mass["Ho"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].flow_mass["Er"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].flow_mass["Tm"]) == pytest.approx(
            0.18, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].flow_mass["Yb"]) == pytest.approx(
            0.29, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].flow_mass["Lu"]) == pytest.approx(
            0.14, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].flow_mass["Th"]) == pytest.approx(
            0, rel=1e-2
        )
        assert value(m.fs.solex.Orgacid[0,1].flow_mass["U"]) == pytest.approx(
            0, rel=1e-2
        )
