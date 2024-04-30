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
from pyomo.environ import Constraint, log, Param, Set, units, Var

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialFlowBasis,
    ReactionParameterBlock,
    ReactionBlockDataBase,
    ReactionBlockBase,
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

        self._reaction_block_class = BrineReactionBlock

        # Reaction Index
        self.equilibrium_reaction_idx = Set(initialize=["P1"])

        # Reaction Stoichiometry
        self.equilibrium_reaction_stoichiometry = {
            ("P1", "liquid", "Li"): 0,
            ("P1", "liquid", "Na"): -1,
            ("P1", "liquid", "Cl"): -1,
            ("P1", "liquid", "H2O"): 0,
        }

        self.solubility_constant = Param(
            self.equilibrium_reaction_idx,
            default={"P1": 8e5},
            units=units.mg**2 / units.L**2,
            mutable=True,
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


class _BrineReactionBlock(ReactionBlockBase):
    """
    This Class contains methods which should be applied to Reaction Blocks as a
    whole, rather than individual elements of indexed Reaction Blocks.
    """

    pass


@declare_process_block_class("BrineReactionBlock", block_class=_BrineReactionBlock)
class BrineReactionBlockData(ReactionBlockDataBase):
    """
    An example reaction package for lithium brine precipitation
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super().build()

        @self.Constraint(self.params.equilibrium_reaction_idx)
        def equilibrium_constraint(b, r):
            return log(b.params.solubility_constant[r]) == sum(
                -b.params.equilibrium_reaction_stoichiometry[r, "liquid", j]
                * log(b.state_ref.conc_mass_comp[j])
                for j in b.state_ref.component_list
                if b.params.equilibrium_reaction_stoichiometry[r, "liquid", j] != 0.0
            )

    def get_reaction_rate_basis(self):
        return MaterialFlowBasis.molar
