#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Property package for precipitate of CMR wastewater treatment.

Authors: Chenyu Wang
"""

from pyomo.environ import Param, Set, Var, units

import idaes.core.util.scaling as iscale
from idaes.core import (
    Solvent,
    Solute,
    Component,
    MaterialFlowBasis,
    Phase,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    declare_process_block_class,
)
from idaes.core.util.initialization import fix_state_vars
from idaes.core.util.misc import add_object_reference


# -----------------------------------------------------------------------------
# Precipitate solution property package
@declare_process_block_class("AqueousParameter")
class AqueousParameterData(PhysicalParameterBlock):
    """
    Property package for aqueous solution generated in precipitator.

    Includes the following components:

    * Rare Earths: Sc, Y, La, Ce, Pr, Nd, Sm, Gd, Dy
    * Impurities: Al, Ca, Fe

    Assumes the equilibrium reaction has a fixed rate and fixed partition coefficients
    based on:

    Wang, Y., Ziemkiewicz, P., Noble, A., A Hybrid Experimental and Theoretical Approach
    to Optimize Recovery of Rare Earth Elements from Acid Mine Drainage,
    Minerals, 2022, 12. 236

    self.split can be substituted by surrogate model
    """

    def build(self):
        super().build()

        self.liquid = Phase()

        # Solvents
        self.H2O = Component()

        # Contaminants
        self.Ca = Component()
        self.Co = Component()
        self.Fe = Component()
        self.Mg = Component()
        self.Mn = Component()
        self.Na = Component()
        self.Ni = Component()
        self.Zn = Component()

        # REE
        self.Dy = Component()
        self.Gd = Component()
        self.Nd = Component()
        self.Pr = Component()

        # anions
        self.SO4 = Component()
        self.OH = Component()

        # parameter based on pH 7
        # TODO add surrogate model/equation
        self.split = Param(
            self.component_list,
            units=units.kg / units.kg,
            initialize={
                "H2O": 1e-20,
                "Ca": 1e-20,
                "Co": 99.9999,
                "Fe": 99.9999,
                "Mg": 1e-20,
                "Mn": 1e-20,
                "Na": 1e-20,
                "Ni": 1e-20,
                "Zn": 52.416,
                "Dy": 1e-20,
                "Gd": 1e-20,
                "Nd": 1e-20,
                "Pr": 1e-20,
                "SO4": 7.304,
                "OH": 1e-20,
            },
        )

        self.mw = Param(
            self.component_list,
            units=units.kg / units.mol,
            initialize={
                "H2O": 18e-3,
                "Ca": 40.08e-3,
                "Co": 58.93e-3,
                "Fe": 55.845e-3,
                "Mg": 24.305e-3,
                "Mn": 54.94e-3,
                "Na": 22.99e-3,
                "Ni": 58.69e-3,
                "Zn": 65.38e-3,
                "Dy": 162.5e-3,
                "Gd": 157.25e-3,
                "Nd": 144.24e-3,
                "Pr": 140.91e-3,
                "SO4": 96.06e-3,
                "OH": 17.01e-3,
            },
        )

        self.dissolved_elements = Set(
            initialize=[
                "Ca",
                "Co",
                "Fe",
                "Mg",
                "Mn",
                "Na",
                "Ni",
                "Zn",
                "Dy",
                "Gd",
                "Nd",
                "Pr",
                "SO4",
                "OH",
                "H2O",
            ]
        )

        # Assume dilute acid, density of pure water
        self.dens_mass = Param(
            initialize=1,
            units=units.kg / units.litre,
            mutable=True,
        )

        self._state_block_class = AqueousStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "flow_vol": {"method": None},
                "conc_mass_comp": {"method": None},
                "dens_mol": {"method": "_dens_mol"},
                "flow_mol_comp": {"method": None},
            }
        )
        obj.add_default_units(
            {
                "time": units.hour,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


class _AqueousStateBlock(StateBlock):
    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        # Fix state variables
        fix_state_vars(self)


@declare_process_block_class("AqueousStateBlock", block_class=_AqueousStateBlock)
class AqueousStateBlockkData(StateBlockData):
    """
    State block for aqueous solution generated in oxalate precipitator.

    """

    def build(self):
        super().build()

        self.flow_vol = Var(
            units=units.L / units.hour,
            initialize=10,
            bounds=(1e-8, None),
        )
        self.conc_mass_comp = Var(
            self.params.dissolved_elements,
            units=units.mg / units.L,
            initialize=1e-5,
            bounds=(1e-20, None),
        )
        self.flow_mol_comp = Var(
            self.params.dissolved_elements,
            units=units.mol / units.hour,
            initialize=1e-5,
            bounds=(1e-20, None),
        )

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

        # Concentration conversion constraint
        @self.Constraint(self.params.dissolved_elements)
        def flow_mol_constraint(b, j):
            return (
                units.convert(
                    b.flow_vol * b.conc_mass_comp[j] / b.params.mw[j],
                    to_units=units.mol / units.hour,
                )
                == b.flow_mol_comp[j]
            )

        iscale.set_scaling_factor(self.flow_vol, 1e1)
        iscale.set_scaling_factor(self.conc_mass_comp, 1e2)
        iscale.set_scaling_factor(self.flow_mol_comp, 1e3)
        iscale.set_scaling_factor(self.conc_mol_comp, 1e5)

    def _dens_mass(self):
        add_object_reference(self, "dens_mass", self.params.dens_mass)

    def get_material_flow_terms(self, p, j):
        # Note conversion to mol/hour
        if j == "H2O":
            # Assume constant density of 1 kg/L
            return self.flow_vol * self.params.dens_mass / self.params.mw[j]
        else:
            # Need to convert from moles to mass
            return units.convert(
                self.flow_vol * self.conc_mass_comp[j] / self.params.mw[j],
                to_units=units.mol / units.hour,
            )

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        return {
            "flow_vol": self.flow_vol,
            "conc_mass_comp": self.conc_mass_comp,
        }
