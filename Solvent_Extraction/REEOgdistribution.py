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
        self.Si = Component()

        # REEs
        self.Sc = Component()
        self.Y = Component()
        self.La = Component()
        self.Ce = Component()
        self.Pr = Component()
        self.Nd = Component()
        self.Pm = Component()
        self.Sm = Component()
        self.Eu = Component()
        self.Gd = Component()
        self.Tb = Component()
        self.Dy = Component()
        self.Ho = Component()
        self.Er = Component()
        self.Tm = Component()
        self.Yb = Component()
        self.Lu = Component()
        self.Th = Component()
        self.U = Component()

        # separate solutes

        self.dissolved_elements = Set(initialize = [
            "Al",
            "Ca",
            "Fe",
            "Si",
            "Sc",
            "Y",
            "La",
            "Ce",
            "Pr",
            "Nd",
            "Pm",
            "Sm",
            "Eu",
            "Gd",
            "Tb",
            "Dy",
            "Ho",
            "Er",
            "Tm",
            "Yb",
            "Lu",
            "Th",
            "U"
        ])

        
        self._state_block_class = REESolExOgStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": units.hour,
                "mass": units.g, 
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

        self.concentration = Var(self.params.dissolved_elements, units=units.mg/units.L, bounds=(0,None))      # Concentration added

        self.flow_vol = Var(units=units.L/units.hour)

        self.organic_vol = Var(units=units.L)

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mass
    
    def get_material_flow_terms(self, j):
        if j in self.params.dissolved_elements:  
            return self.flow_vol * self.concentration[j]/1000                 # Concentration added
        elif j=="DEHPA":
            return self.flow_vol * (975.8 * units.g/units.L) 
        else:
            raise BurntToast()
        
    def define_state_vars(self):
        return {
            "flow_vol": self.flow_vol,
            "concentration": self.concentration,
            "organic_vol":self.organic_vol
        }

