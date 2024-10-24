#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Initial property package for precipitate.

Authors: Alejandro Garciadiego
"""

from pyomo.environ import Param, Set, Var, units, exp

import idaes.core.util.scaling as iscale
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
from idaes.core.util.misc import add_object_reference


# -----------------------------------------------------------------------------
# Precipitate solution property package
@declare_process_block_class("AqueousParameter")
class AqueousParameterData(PhysicalParameterBlock):
    """
    Property package for aqueous solution generated in oxalate precipitator.

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
        self.Al = Component()
        self.Ca = Component()
        self.Fe = Component()
        self.H = Component()
        self.Cl = Component()
        self.HSO4 = Component()
        self.SO4 = Component()

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

        # parameter based on pH 1.5
        # TODO add surrogate model/equation
        self.split = Var(
            self.component_list,
            units=units.kg / units.kg,
            initialize={
                "H2O": 1e-20,
                "Sc": 0.3161,
                "Y": 0.7446,
                "La": 0.5151,
                "Ce": 0.6807,
                "Pr": 0.78,
                "Nd": 0.8155,
                "Sm": 0.8735,
                "Gd": 0.8801,
                "Dy": 0.8716,
                "Al": 0.009,
                "Ca": 0.0001,
                "Fe": 0.0244,
                "H": 1e-20,
                "Cl": 1e-20,
                "HSO4": 1e-20,
                "SO4": 1e-20,
            },
            bounds=(1e-30, 1),
        )

        self.acid_flow = Var(
            units=units.dimensionless,
            initialize=6.4,
            bounds=(1e-6, 100),
        )

        # parameter based on pH 1.5
        # TODO add surrogate model/equation
        self.E_D = Param(
            self.component_list,
            units=units.dimensionless,
            initialize={
                "H2O": 100,
                "Sc": 6.42030,
                "Y": 4.551786,
                "La": 4.3717,
                "Ce": 1.18848,
                "Pr": 2.09604,
                "Nd": 1.01030,
                "Sm": 2.296176,
                "Gd": 3.07276,
                "Dy": 4.8608,
                "Al": 50,
                "Ca": 14.49274,
                "Fe": 8.659561,
                "H": 100,
                "Cl": 100,
                "HSO4": 100,
                "SO4": 100,
            },
        )

        # parameter based on pH 1.5
        # TODO add surrogate model/equation
        self.N_D = Param(
            self.component_list,
            units=units.dimensionless,
            initialize={
                "H2O": 100,
                "Sc": 6.42030,
                "Y": 4.67403,
                "La": 4.6340,
                "Ce": 2.737238,
                "Pr": 3.44364,
                "Nd": 2.419137,
                "Sm": 3.7201,
                "Gd": 4.1995,
                "Dy": 4.73106,
                "Al": 0.9,
                "Ca": 4.45302,
                "Fe": 3.6495,
                "H": 100,
                "Cl": 100,
                "HSO4": 100,
                "SO4": 100,
            },
        )

        self.mw = Param(
            self.component_list,
            units=units.kg / units.mol,
            initialize={
                "H2O": 18e-3,
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
                "H": 1.008e-3,
                "Cl": 35.453e-3,
                "HSO4": 97.064e-3,
                "SO4": 96.056e-3,
            },
        )

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
                "H",
                "Cl",
                "HSO4",
                "SO4",
                "H2O",
            ]
        )

        self.split_elements = Set(
            initialize=[
                "Al",
                # "Ca",
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

        iscale.set_scaling_factor(self.flow_vol, 1e-2)
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
