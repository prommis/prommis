from pyomo.environ import Param, Set, Var, units

from idaes.core import (
    Component,
    MaterialFlowBasis,
    Phase,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    declare_process_block_class,
)
from idaes.core.util.initialization import fix_state_vars


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

        self.mw = Param(
            self.component_list,
            units=units.kg / units.mol,
            initialize={
                "DEHPA": 322.431e-3,
                "Sc": 44.946e-3,
                "Y": 88.905e-3,
                "La": 138.905e-3,
                "Ce": 140.116e-3,
                "Pr": 140.907e-3,
                "Nd": 144.242e-3,
                "Sm": 150.36e-3,
                "Gd": 157.25e-3,
                "Dy": 162.50e-3,
                "Al": 26.982e-3,
                "Ca": 40.078e-3,
                "Fe": 55.845e-3,
            },
        )

        # density of DEHPA
        self.dens_mol = Param(
            initialize=975.8e-3,
            units=units.kg / units.litre,
            mutable=True,
        )

        self._state_block_class = REESolExOgStateBlock

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


class _REESolExOgStateBlock(StateBlock):
    def fix_initialization_states(self):
        fix_state_vars(self)


@declare_process_block_class("REESolExOgStateBlock", block_class=_REESolExOgStateBlock)
class REESolExOgStateBlockData(StateBlockData):
    def build(self):
        super().build()

        self.conc_mass_comp = Var(
            self.params.dissolved_elements,
            units=units.mg / units.L,
            bounds=(1e-20, None),
        )

        self.flow_vol = Var(units=units.L / units.hour, bounds=(1e-8, None))

        self.conc_mol_comp = Var(
            self.params.dissolved_elements,
            units=units.mol / units.L,
            initialize=1e-5,
            bounds=(1e-20, None),
        )

        # Concentration conversion constraint
        @self.Constraint(self.params.dissolved_elements)
        def molar_concentration_constraint(b, j):
            return (
                units.convert(
                    b.conc_mol_comp[j] * b.params.mw[j], to_units=units.mg / units.litre
                )
                == b.conc_mass_comp[j]
            )

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def get_material_flow_terms(self, p, j):
        if j == "DEHPA":
            return self.flow_vol * self.params.dens_mol / self.params.mw[j]
        else:
            return units.convert(
                self.flow_vol * self.conc_mass_comp[j] / self.params.mw[j],
                to_units=units.mol / units.hour,
            )

    def define_state_vars(self):
        return {"flow_vol": self.flow_vol, "conc_mass_comp": self.conc_mass_comp}
