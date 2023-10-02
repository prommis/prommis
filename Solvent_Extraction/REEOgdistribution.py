# Imports 

from pyomo.environ import (Param, Set, units, Var)
from idaes.core import (declare_process_block_class, PhysicalParameterBlock, StateBlock,
                        StateBlockData, Component, Phase, MaterialFlowBasis)
from idaes.core.util.initialization import fix_state_vars
from idaes.core.util.exceptions import BurntToast


@declare_process_block_class("REESolExOgParameters")
class REESolExOgParameterData(PhysicalParameterBlock):
    def build(self):
        super().build()

        self.liquid = Phase()

        # Solvents
        self.DEHPA = Component()

        # Contaminants
        self.Al = Component()
        self.Ca = Component()
        self.Fe = Component()
        # self.Si = Component()

        # REEs
        self.Sc = Component()
        self.Y = Component()
        self.La = Component()
        self.Ce = Component()
        self.Pr = Component()
        self.Nd = Component()
        # self.Pm = Component()
        self.Sm = Component()
        # self.Eu = Component()
        self.Gd = Component()
        # self.Tb = Component()
        self.Dy = Component()
        # self.Ho = Component()
        # self.Er = Component()
        # self.Tm = Component()
        # self.Yb = Component()
        # self.Lu = Component()
        # self.Th = Component()
        # self.U = Component()

        # separate solutes

        self.dissolved_elements = Set(initialize=[
            "Al",
            "Ca",
            "Fe",
            # "Si",
            "Sc",
            "Y",
            "La",
            "Ce",
            "Pr",
            "Nd",
            # "Pm",
            "Sm",
            # "Eu",
            "Gd",
            # "Tb",
            "Dy",
            # "Ho",
            # "Er",
            # "Tm",
            # "Yb",
            # "Lu",
            # "Th",
            # "U"
        ])

        self.K_distribution = Param(self.dissolved_elements, initialize={
            "Al": 3.6 / 100,
            "Ca": 3.7 / 100,
            "Fe": 2.1 / 100,
            # "Si":0/100,
            "Sc": 100 / 100,
            "Y": 100 / 100,
            "La": 75.2 / 100,
            "Ce": 95.7 / 100,
            "Pr": 96.5 / 100,
            "Nd": 99.2 / 100,
            # "Pm":100/100,
            "Sm": 100 / 100,
            # "Eu":99.9/100,
            "Gd": 98.6 / 100,
            # "Tb":99.3/100,
            "Dy": 99.9 / 100,
            # "Ho":99.5/100,
            # "Er":99.5/100,
            # "Tm":98.6/100,
            # "Yb":80.7/100,
            # "Lu":99.5/100,
            # "Th":5/100,
            # "U":99.5/100
        }, mutable=True,
                                    units=units.dimensionless,
                                    doc="The fraction of component that goes from aqueous to organic phase")

        self._state_block_class = REESolExOgStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": units.hour,
                "mass": units.kg,
                "amount": units.mol,
                "length": units.m,
                "temperature": units.K
            }
        )


class _REESolExOgStateBlock(StateBlock):
    def fix_initialization_states(self):
        fix_state_vars(self)


@declare_process_block_class("REESolExOgStateBlock", block_class=_REESolExOgStateBlock)
class REESolExOgStateBlockData(StateBlockData):
    def build(self):
        super().build()

        self.flow_mass = Var(self.params.dissolved_elements, units=units.g / units.hour)

        self.flow_vol = Var(units=units.L / units.hour)

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mass

    def get_material_flow_terms(self, j):
        if j in self.params.dissolved_elements:
            return self.flow_mass[j]
        elif j == "DEHPA":
            return self.flow_vol * (975.8 * units.g / units.L)
        else:
            raise BurntToast()

    def define_state_vars(self):
        return {
            "flow_vol": self.flow_vol,
            "flow_mass": self.flow_mass
        }

