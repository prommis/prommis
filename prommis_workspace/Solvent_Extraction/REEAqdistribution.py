from pyomo.environ import Set, units, Var
from idaes.core import (
    declare_process_block_class,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    Component,
    Phase,
    MaterialFlowBasis,
)
from idaes.core.util.initialization import fix_state_vars
from idaes.core.util.exceptions import BurntToast


@declare_process_block_class("REESolExAqParameters")
class REESolExAqParameterData(PhysicalParameterBlock):
    def build(self):
        super().build()

        self.liquid = Phase()

        # Solvents
        self.H2SO4 = Component()

        # Contaminants
        self.Al = Component()
        self.Ca = Component()
        self.Fe = Component()

        # REEs
        self.Sc = Component()
        self.Y = Component()
        self.La = Component()
        self.Ce = Component()
        self.Pr = Component()
        self.Nd = Component()
        self.Sm = Component()
        self.Gd = Component()
        self.Dy = Component()

        self.dissolved_elements = Set(
            initialize=[
                "Al",
                "Ca",
                "Fe",
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
        )

        self._state_block_class = REESolExAqStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": units.hour,
                "mass": units.kg,
                "amount": units.mol,
                "length": units.m,
                "temperature": units.K,
            }
        )


class _REESolExAqStateBlock(StateBlock):
    def fix_initialization_states(self):
        fix_state_vars(self)


@declare_process_block_class("REESolExAqStateBlock", block_class=_REESolExAqStateBlock)
class REESolExAqStateBlockData(StateBlockData):
    def build(self):
        super().build()

        self.conc_mass_comp = Var(
            self.params.dissolved_elements,
            units=units.mg / units.L,
            bounds=(1e-20, None),
        )

        self.flow_vol = Var(units=units.L / units.hour, bounds=(1e-8, None))

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mass

    def get_material_flow_terms(self, p, j):
        if j in self.params.dissolved_elements:
            return self.flow_vol * units.convert(
                self.conc_mass_comp[j], to_units=units.kg / units.L
            )
        elif j == "H2SO4":
            return self.flow_vol * (1840 * units.g / units.L)
        else:
            raise BurntToast()

    def define_state_vars(self):
        return {"flow_vol": self.flow_vol, "conc_mass_comp": self.conc_mass_comp}
