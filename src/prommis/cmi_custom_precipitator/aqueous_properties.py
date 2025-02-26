#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Aqueous component property package for optimization-based precipitator model.

Authors: Chris Laliwala
"""

from pyomo.common.config import ConfigValue
from pyomo.environ import Set, Var
from pyomo.environ import units as pyunits

import idaes.logger as idaeslog

# Import IDAES cores
from idaes.core import (
    MaterialFlowBasis,
    Phase,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    declare_process_block_class,
)
from idaes.core.util.initialization import fix_state_vars

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("AqueousParameter")
class AqueousParameterData(PhysicalParameterBlock):
    CONFIG = PhysicalParameterBlock.CONFIG()
    CONFIG.declare(
        "aq_comp_list",
        ConfigValue(domain=list, description="List of aqueous components in process"),
    )
    CONFIG.declare(
        "eq_rxn_logkeq_dict",
        ConfigValue(
            domain=dict, description="Dictionary of equilibrium reaction constants"
        ),
    )
    CONFIG.declare(
        "eq_rxn_stoich_dict",
        ConfigValue(
            domain=dict, description="Dictionary of equilibrium reaction stoichiometry"
        ),
    )

    def build(self):
        """
        Calllable method for Block construction.
        """
        # super(AqueousParameterData, self).build()
        super().build()

        self.AqueousPhase = Phase()
        self.component_list = self.config.aq_comp_list

        ## equilibrium reaction parameters
        # equilibrium reaction index
        self.eq_rxn_set = Set(
            initialize=list(set(key for key in self.config.eq_rxn_logkeq_dict.keys()))
        )

        # stoichiometry for each equilibrium reaction
        self.eq_rxn_stoich_dict = self.config.eq_rxn_stoich_dict

        # log(keq) for each equilibrium reaction
        self.eq_rxn_logkeq_dict = self.config.eq_rxn_logkeq_dict

        self._state_block_class = AqueousStateBlock

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        obj.define_custom_properties(
            {
                "molality_aq_comp": {"method": None},
                "flow_vol": {"method": None},
            }
        )
        obj.add_default_units(
            {
                "time": pyunits.hour,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
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
        return fix_state_vars(self)


@declare_process_block_class("AqueousStateBlock", block_class=_AqueousStateBlock)
class AqueousStateBlockData(StateBlockData):
    def build(self):
        super().build()

        self.flow_vol = Var(
            units=pyunits.kg / pyunits.hour,
            initialize=1,
            bounds=(1e-5, None),
        )

        self.molality_aq_comp = Var(
            self.component_list,
            units=pyunits.mol / pyunits.kg,
            initialize=1e-20,
            bounds=(1e-20, None),
        )

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        return {
            "flow_vol": self.flow_vol,
            "molality_aq_comp": self.molality_aq_comp,
        }
