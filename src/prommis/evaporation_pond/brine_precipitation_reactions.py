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
Example property package for precipitation reactions associated with lithium brines.
"""

# Import Pyomo libraries
from pyomo.environ import Param, Set, units

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    ReactionParameterBlock,
)


# Some more information about this module
__author__ = "Andrew Lee"


@declare_process_block_class("BrineReactionParameters")
class BrineReactionParametersData(ReactionParameterBlock):
    """
    Property parameters associated with lithium brine precipitation reactions.

    """

    def build(self):
        """
        Callable method for Block construction.
        """
        super().build()

        # For the case of evaporation ponds, we don't actually need a ReactionBlock
        # so no need to create a link here.

        # Reaction Index
        self.equilibrium_reaction_idx = Set(initialize=["P1", "P2"])

        # Reaction Stoichiometry
        self.equilibrium_reaction_stoichiometry = {
            ("P1", "liquid", "Li"): 0,
            ("P1", "liquid", "Na"): -1,
            ("P1", "liquid", "Cl"): -1,
            ("P1", "liquid", "H2O"): 0,
            ("P2", "liquid", "Li"): -1,
            ("P2", "liquid", "Na"): 0,
            ("P2", "liquid", "Cl"): -1,
            ("P2", "liquid", "H2O"): 0,
        }

        # Units of solubility constants may differ depending on stoichiometry
        # Need separate parameters for each
        self.solubility_constant_P1 = Param(
            default=1.5e6,
            units=units.mg**2 / units.L**2,
            mutable=True,
            doc="Solubility constant for reaction P1",
        )
        self.solubility_constant_P2 = Param(
            default=3e6,
            units=units.mg**2 / units.L**2,
            mutable=True,
            doc="Solubility constant for reaction P2",
        )

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


# Evaporation ponds will create their own equilibrium constraints,
# so we don't need a ReactionBlock here.
