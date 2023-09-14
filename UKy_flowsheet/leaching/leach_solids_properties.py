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
Initial property package for West Kentucky No. 13 coal refuse.

Authors: Andrew Lee
"""

from pyomo.environ import (
    Constraint,
    Param,
    Set,
    units,
    Var,
)

from idaes.core import (
    declare_process_block_class,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    Component,
    Phase,
    MaterialFlowBasis,
    MaterialBalanceType,
)
from idaes.core.util.exceptions import BurntToast
from idaes.core.util.initialization import fix_state_vars


# -----------------------------------------------------------------------------
# Leach solution property package
@declare_process_block_class("CoalRefuseParameters")
class CoalRefuseParameterData(PhysicalParameterBlock):
    def build(self):
        super().build()

        self.solid = Phase()

        # Solids
        self.inerts = Component()

        # REEs
        self.Sc2O3 = Component()
        self.Y2O3 = Component()
        self.La2O3 = Component()
        self.Ce2O3 = Component()
        self.Pr2O3 = Component()
        self.Nd2O3 = Component()
        self.Sm2O3 = Component()
        self.Gd2O3 = Component()
        self.Dy2O3 = Component()

        # Contaminants
        self.Al2O3 = Component()
        self.CaO = Component()
        self.Fe2O3 = Component()

        self.mw = Param(
            self.component_list,
            units=units.kg / units.mol,
            initialize={
                "inerts": 60.08e-3,
                "Sc2O3": (44.946*2 + 3*15.999)*1e-3,
                "Y2O3": (88.905*2 + 3*15.999)*1e-3,
                "La2O3": (138.905*2 + 3*15.999)*1e-3,
                "Ce2O3": (140.116*2 + 3*15.999)*1e-3,
                "Pr2O3": (140.907*2 + 3*15.999)*1e-3,
                "Nd2O3": (144.242*2 + 3*15.999)*1e-3,
                "Sm2O3": (150.36*2 + 3*15.999)*1e-3,
                "Gd2O3": (157.25*2 + 3*15.999)*1e-3,
                "Dy2O3": (162.50*2 + 3*15.999)*1e-3,
                "Al2O3": (26.982*2 + 3*15.999)*1e-3,
                "CaO": (40.078 + 15.999)*1e-3,
                "Fe2O3": (55.845*2 + 3*15.999)*1e-3,
            },
        )

        self._state_block_class = CoalRefuseStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": units.hour,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


class _CoalRefuseStateBlock(StateBlock):
    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        # Fix state variables
        fix_state_vars(self)

        if not self.config.defined_state:
            self.sum_mass_frac.deactivate()


@declare_process_block_class(
    "CoalRefuseStateBlock", block_class=_CoalRefuseStateBlock
)
class CoalRefuseStateBlockData(StateBlockData):
    def build(self):
        super().build()

        self.flow_mass = Var(
            units=units.kg / units.hour,
            bounds=(1e-8, None),
        )
        self.mass_frac_comp = Var(
            self.params.component_list,
            units=units.kg / units.kg,
            bounds=(1e-8, None),
        )

        if not self.config.defined_state:
            self.sum_mass_frac = Constraint(
                expr=1 == sum(self.mass_frac_comp[j] for j in self.component_list)
            )

    def get_material_flow_terms(self, p, j):
        return self.flow_mass * self.mass_frac_comp[j] / self.params.mw[j]

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        return {
            "flow_mass": self.flow_mass,
            "mass_frac_comp": self.mass_frac_comp,
        }
