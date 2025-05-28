# Import Pyomo libraries
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
    Translator block representing the SX/leaching interface
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
