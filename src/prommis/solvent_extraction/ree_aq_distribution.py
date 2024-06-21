#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Initial property package for the aqueous phase solution of the solvent extraction
unit operation.

Authors: Arkoprabho Dasgupta

"""

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


@declare_process_block_class("REESolExAqParameters")
class REESolExAqParameterData(PhysicalParameterBlock):
    """
    This is a property package for the aqueous phase solution of the solvent extraction
    unit operation of the University of Kentucky pilot plant flowsheet.

    This  includes the following components:

    * Solvent: H2O
    * Acid components: H, SO4, HSO4, H2SO4
    * Rare Earths: Sc, Y, La, Ce, Pr, Nd, Sm, Gd, Dy
    * Impurities: Al, Ca, Fe

    Any kind of reactive interactions of the acid components are not considered.

    """

    def build(self):
        super().build()

        self.liquid = Phase()

        # Solvents
        self.H2O = Component()

        # Inerts
        self.H2SO4 = Component()
        self.H = Component()
        self.SO4 = Component()
        self.HSO4 = Component()

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
                "H2SO4",
                "H",
                "SO4",
                "HSO4",
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
                "H2SO4": 98.079e-3,
                "H2O": 18.015e-3,
                "H": 1.008e-3,
                "SO4": 96.06e-3,
                "HSO4": 97.068e-3,
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

        # density of H2O
        self.dens_mass = Param(
            initialize=1,
            units=units.kg / units.litre,
            mutable=True,
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
    """
    State block for aqueous phase solution of the solvent extraction process.

    """

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
        if j == "H2O":
            return self.flow_vol * self.params.dens_mass / self.params.mw[j]
        else:
            return units.convert(
                self.flow_vol * self.conc_mass_comp[j] / self.params.mw[j],
                to_units=units.mol / units.hour,
            )

    def get_material_density_terms(self, p, j):
        if j == "H2O":
            return units.convert(
                self.params.dens_mass / self.params.mw[j],
                to_units=units.mol / units.m**3,
            )
        else:
            return units.convert(
                self.conc_mass_comp[j] / self.params.mw[j],
                to_units=units.mol / units.m**3,
            )

    def define_state_vars(self):
        return {"flow_vol": self.flow_vol, "conc_mass_comp": self.conc_mass_comp}
