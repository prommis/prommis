#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Initial property package for the organic phase solution of the solvent extraction
unit operation.

Authors: Arkoprabho Dasgupta

"""

from pyomo.environ import Param, Var, units, Constraint, PositiveReals, Reals

from idaes.core import (
    Component,
    MaterialFlowBasis,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    declare_process_block_class,
    LiquidPhase,
)
from idaes.core.util.initialization import fix_state_vars
from idaes.core.util.misc import add_object_reference


@declare_process_block_class("REESolExOgParameters")
class REESolExOgParameterData(PhysicalParameterBlock):
    """
    This is a property package for the organic phase solution of the solvent extraction
    unit operation of the University of Kentucky pilot plant flowsheet.

    This  includes the following components:

    * Solvent: Kerosene
    * Extractant: DEHPA
    * Rare Earths: Sc, Y, La, Ce, Pr, Nd, Sm, Gd, Dy
    * Impurities: Al, Ca, Fe

    Kerosene is not considered to be involved in any reaction.

    """

    def build(self):
        super().build()

        self.organic = LiquidPhase()

        # Solvent
        self.Kerosene = Component()

        # Extractant
        self.DEHPA = Component()

        # Contaminants
        self.Al_o = Component()
        self.Ca_o = Component()
        self.Fe_o = Component()

        # REEs
        self.Sc_o = Component()
        self.Y_o = Component()
        self.La_o = Component()
        self.Ce_o = Component()
        self.Pr_o = Component()
        self.Nd_o = Component()
        self.Sm_o = Component()
        self.Gd_o = Component()
        self.Dy_o = Component()

        self.mw = Param(
            self.component_list,
            units=units.kg / units.mol,
            initialize={
                "Kerosene": 170e-3,
                "DEHPA": 322.431e-3,
                "Sc_o": 44.946e-3,
                "Y_o": 88.905e-3,
                "La_o": 138.905e-3,
                "Ce_o": 140.116e-3,
                "Pr_o": 140.907e-3,
                "Nd_o": 144.242e-3,
                "Sm_o": 150.36e-3,
                "Gd_o": 157.25e-3,
                "Dy_o": 162.50e-3,
                "Al_o": 26.982e-3,
                "Ca_o": 40.078e-3,
                "Fe_o": 55.845e-3,
            },
        )

        # density of Kerosene
        self.dens_mass = Param(
            initialize=0.82,
            units=units.kg / units.litre,
            mutable=True,
        )

        # Heat capacity of kerosene
        self.cp_mol = Param(
            mutable=True,
            initialize=341.7,
            doc="Molar heat capacity of kerosene [J/mol.K]",
            units=units.J / units.mol / units.K,
        )

        self.temperature_ref = Param(
            within=PositiveReals,
            mutable=True,
            default=298.15,
            doc="Reference temperature [K]",
            units=units.K,
        )

        self._state_block_class = REESolExOgStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "flow_vol": {"method": None},
                "conc_mass_comp": {"method": None},
                "conc_mol_comp": {"method": None},
                "dens_mass": {"method": "_dens_mass"},
            }
        )
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
    """
    State block for organic phase solution of the solvent extraction process.

    """

    def build(self):
        super().build()

        self.conc_mass_comp = Var(
            self.params.component_list,
            units=units.mg / units.L,
            initialize=1e-8,
            bounds=(1e-20, None),
        )

        self.flow_vol = Var(units=units.L / units.hour, bounds=(1e-8, None))

        self.conc_mol_comp = Var(
            self.params.component_list,
            units=units.mol / units.L,
            initialize=1e-5,
            bounds=(1e-20, None),
        )

        self.temperature = Var(
            domain=Reals,
            initialize=298.15,
            bounds=(298.1, None),
            doc="State temperature [K]",
            units=units.K,
        )

        self.pressure = Var(
            domain=Reals,
            initialize=101325.0,
            bounds=(1e3, 1e6),
            doc="State pressure [Pa]",
            units=units.Pa,
        )

        # Concentration conversion constraint
        @self.Constraint(self.params.component_list)
        def molar_concentration_constraint(b, j):
            return (
                units.convert(
                    b.conc_mol_comp[j] * b.params.mw[j], to_units=units.mg / units.litre
                )
                == b.conc_mass_comp[j]
            )

        if not self.config.defined_state:
            # Concentration of kerosene based on assumed density
            self.kerosene_concentration = Constraint(
                expr=self.conc_mass_comp["Kerosene"] == 8.2e5 * units.mg / units.L
            )

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def _dens_mass(self):
        add_object_reference(self, "dens_mass", self.params.dens_mass)

    def get_material_flow_terms(self, p, j):
        if j == "Kerosene":
            return self.flow_vol * self.params.dens_mass / self.params.mw[j]
        else:
            return units.convert(
                self.flow_vol * self.conc_mass_comp[j] / self.params.mw[j],
                to_units=units.mol / units.hour,
            )

    def get_material_density_terms(self, p, j):
        if j == "Kerosene":
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
        return {
            "flow_vol": self.flow_vol,
            "conc_mass_comp": self.conc_mass_comp,
            "temperature": self.temperature,
            "pressure": self.pressure,
        }
