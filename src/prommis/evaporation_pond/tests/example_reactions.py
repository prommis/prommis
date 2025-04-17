#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Example reaction package for testing evaporation pond unit model.

Authors: Andrew Lee
"""

# Import Pyomo libraries
from pyomo.environ import Param, Set, units

# Import IDAES cores
from idaes.core import (
    MaterialFlowBasis,
    ReactionBlockBase,
    ReactionBlockDataBase,
    ReactionParameterBlock,
    declare_process_block_class,
)
from idaes.core.util.misc import add_object_reference

# Some more information about this module
__author__ = "Andrew Lee"


@declare_process_block_class("BrineReactionParameters")
class BrineReactionParametersData(ReactionParameterBlock):
    """
    Reaction parameters associated with lithium brine precipitation.

    """

    def build(self):
        """
        Callable method for Block construction.
        """
        super().build()

        self._reaction_block_class = BrineReactionBlock

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
        self.solubility_product_P1 = Param(
            default=1.84e-3,
            units=units.mol**2 / units.L**2,
            mutable=True,
            doc="Solubility constant for reaction P1",
        )
        self.solubility_product_P2 = Param(
            default=1.22e-2,
            units=units.mol**2 / units.L**2,
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


@declare_process_block_class("BrineReactionBlock", block_class=ReactionBlockBase)
class BrineReactionData(ReactionBlockDataBase):
    def build(self):
        """
        Reaction block for brine precipitation.

        """
        super().build()

        for k in self.params.equilibrium_reaction_idx:
            Ksp = getattr(self.params, "solubility_product_" + k)
            add_object_reference(self, "solubility_product_" + k, Ksp)

    def get_reaction_rate_basis(b):
        return MaterialFlowBasis.molar
