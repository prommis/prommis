#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Property package for the inlet and outlet streams of the multi-component diafiltration membrane.

Author: Molly Dougher
"""

from pyomo.common.config import ConfigValue
from pyomo.environ import Var, units

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


@declare_process_block_class("MultiComponentDiafiltrationStreamParameter")
class MultiComponentDiafiltrationStreamParameterData(PhysicalParameterBlock):
    """
    Property Package for the feed and product streams in the multi-component
    diafiltration membrane.

    Currently includes the following solutes:
        Li+ (lithium ion)
        Co2+ (cobalt ion)
        Al3+ (aluminum ion)
        Cl- (chloride ion)
    """

    CONFIG = PhysicalParameterBlock.CONFIG()

    CONFIG.declare(
        "num_salts",
        ConfigValue(
            default=2,
            doc="Number of salts to be modeled",
        ),
    )

    def build(self):
        super().build()

        self.liquid = Phase()

        # add cations
        self.cation_1 = Component()
        if self.config.num_salts > 1:
            self.cation_2 = Component()
        if self.config.num_salts > 2:
            self.cation_3 = Component()

        # add anion
        self.anion = Component()

        self._state_block_class = MultiComponentDiafiltrationStreamStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "flow_vol": {"method": None},
                "conc_mol_comp": {"method": None},
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


class _MultiComponentDiafiltrationStreamStateBlock(StateBlock):
    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        fix_state_vars(self)


@declare_process_block_class(
    "MultiComponentDiafiltrationStreamStateBlock",
    block_class=_MultiComponentDiafiltrationStreamStateBlock,
)
class MultiComponentDiafiltrationStreamStateBlockData(StateBlockData):
    """
    State block for the feed and product streams in the multi-component
    diafiltration membrane.
    """

    def build(self):
        super().build()

        self.flow_vol = Var(
            units=units.m**3 / units.h,
            initialize=10,
            bounds=(1e-20, None),
        )
        self.conc_mol_comp = Var(
            self.component_list,
            units=units.mol / units.m**3,
            initialize=1e-5,
            bounds=(1e-20, None),
        )

    def get_material_flow_terms(self, p, j):
        return self.flow_vol * self.conc_mol_comp[j]

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mole

    def define_state_vars(self):
        return {"flow_vol": self.flow_vol, "conc_mol_comp": self.conc_mol_comp}
