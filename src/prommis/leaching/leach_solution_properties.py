#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Initial property package for REE leach solutions from coal refuse.

Authors: Andrew Lee
"""

from pyomo.environ import Constraint, Param, Set, Var, units

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
# Leach solution property package
@declare_process_block_class("LeachSolutionParameters")
class LeachSolutionParameterData(PhysicalParameterBlock):
    def build(self):
        super().build()

        self.liquid = Phase()

        # Solvent
        self.H2O = Component()

        # Acid related species
        self.H = Component()
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

        # Contaminants
        self.Al = Component()
        self.Ca = Component()
        self.Fe = Component()

        self.mw = Param(
            self.component_list,
            units=units.kg / units.mol,
            initialize={
                "H2O": 18e-3,
                "H": 1e-3,
                "HSO4": 97e-3,
                "SO4": 96e-3,
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

        # Inherent reaction for partial dissociation of HSO4
        self._has_inherent_reactions = True
        self.inherent_reaction_idx = Set(initialize=["Ka2"])
        self.inherent_reaction_stoichiometry = {
            ("Ka2", "liquid", "H"): 1,
            ("Ka2", "liquid", "HSO4"): -1,
            ("Ka2", "liquid", "SO4"): 1,
            ("Ka2", "liquid", "H2O"): 0,
            ("Ka2", "liquid", "Sc"): 0,
            ("Ka2", "liquid", "Y"): 0,
            ("Ka2", "liquid", "La"): 0,
            ("Ka2", "liquid", "Ce"): 0,
            ("Ka2", "liquid", "Pr"): 0,
            ("Ka2", "liquid", "Nd"): 0,
            ("Ka2", "liquid", "Sm"): 0,
            ("Ka2", "liquid", "Gd"): 0,
            ("Ka2", "liquid", "Dy"): 0,
            ("Ka2", "liquid", "Al"): 0,
            ("Ka2", "liquid", "Ca"): 0,
            ("Ka2", "liquid", "Fe"): 0,
        }
        self.Ka2 = Param(
            initialize=10**-1.99,
            mutable=True,
            units=units.mol / units.L,
        )

        # Assume dilute acid, density of pure water
        self.dens_mol = Param(
            initialize=1,
            units=units.kg / units.litre,
            mutable=True,
        )

        self._state_block_class = LeachSolutionStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "flow_vol": {"method": None},
                "conc_mass_comp": {"method": None},
                "conc_mol_comp": {"method": None},
                "dens_mol": {"method": "_dens_mol"},
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


class _LeachSolutionStateBlock(StateBlock):
    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        # Fix state variables
        fix_state_vars(self)

        # Deactivate inherent reactions
        for sbd in self.values():
            if not sbd.config.defined_state:
                sbd.h2o_concentration.deactivate()
                sbd.hso4_dissociation.deactivate()


@declare_process_block_class(
    "LeachSolutionStateBlock", block_class=_LeachSolutionStateBlock
)
class LeachSolutionStateBlockData(StateBlockData):
    def build(self):
        super().build()

        self.flow_vol = Var(
            units=units.L / units.hour,
            bounds=(1e-8, None),
        )
        self.conc_mass_comp = Var(
            self.params.component_list,
            units=units.mg / units.L,
            bounds=(1e-20, None),
        )
        self.conc_mol_comp = Var(
            self.params.component_list,
            units=units.mol / units.L,
            initialize=1e-5,
            bounds=(1e-20, None),
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
            # Concentration of H2O based on assumed density
            self.h2o_concentration = Constraint(
                expr=self.conc_mass_comp["H2O"] == 1e6 * units.mg / units.L
            )
            # Equilibrium for partial dissociation of HSO4
            self.hso4_dissociation = Constraint(
                expr=self.conc_mol_comp["HSO4"] * self.params.Ka2
                == self.conc_mol_comp["H"] * self.conc_mol_comp["SO4"]
            )

    def _dens_mol(self):
        add_object_reference(self, "dens_mol", self.params.dens_mol)

    def get_material_flow_terms(self, p, j):
        # Note conversion to mol/hour
        if j == "H2O":
            # Assume constant density of 1 kg/L
            return self.flow_vol * self.params.dens_mol / self.params.mw[j]
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
