#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

r"""
Translator block representing the precipitation/leaching interface.

========================================================================

Author: Arkoprabho Dasgupta

Model description
-----------------

This translator block helps to translate between the property packages of leaching solution phase and
precipitation liquid phase. This block is built upon the traditional Translator block with some
additional constraints.

Additional constraints
----------------------

1. eq_flow_vol_rule = Constraint which equates the inlet and the outlet volumetric flowrates.
2. eq_conc_mass_metals = Constraint which equates the inlet and the outlet concentrations of the elements.

"""

from pyomo.environ import Set

# Import IDAES cores
from idaes.core import declare_process_block_class
from idaes.models.unit_models.translator import TranslatorData

import idaes.logger as idaeslog

__author__ = "Arkoprabho Dasgupta"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("TranslatorLeachPrecip")
class TranslatorDataLeachPrecip(TranslatorData):
    """
    Translator block representing the interface of the leaching solution property
    package and the precipitation liquid property package.

    """

    def build(self):
        """
        Begin building model.
        Args:
            None
        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super().build()

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality volumetric flow equation",
        )
        def eq_flow_vol_rule(blk, t):
            return blk.properties_out[t].flow_vol == blk.properties_in[t].flow_vol

        self.components = Set(
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
                "H2O",
                "H",
                "SO4",
                "HSO4",
                "Cl",
            ]
        )

        @self.Constraint(
            self.flowsheet().time,
            self.components,
            doc="Equality equation for metal components",
        )
        def eq_conc_mass_metals(blk, t, i):
            return (
                blk.properties_out[t].conc_mass_comp[i]
                == blk.properties_in[t].conc_mass_comp[i]
            )
