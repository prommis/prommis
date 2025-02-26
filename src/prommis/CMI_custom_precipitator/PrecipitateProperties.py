#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Precipitate property package for optimization-based precipitator model.

Authors: Chris Laliwala
"""

from pyomo.common.config import ConfigValue
from pyomo.environ import Param, Var, Set, units as pyunits

import idaes.core.util.scaling as iscale
from idaes.core import (
    MaterialFlowBasis,
    Phase,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    declare_process_block_class,
)
from idaes.core.util.initialization import fix_state_vars

# -----------------------------------------------------------------------------
# Precipitate solution property package


@declare_process_block_class("PrecipitateParameter")
class PrecipitateParameterData(PhysicalParameterBlock):
    CONFIG = PhysicalParameterBlock.CONFIG()
    CONFIG.declare(
        "precip_comp_list",
        ConfigValue(
            domain=list, description="List of precipitate components in process."
        ),
    )
    CONFIG.declare(
        "precip_eq_rxn_logkeq_dict",
        ConfigValue(
            domain=dict,
            description="Dictionary of precipitation equilibrium reaction constants",
        ),
    )
    CONFIG.declare(
        "precip_eq_rxn_stoich_dict",
        ConfigValue(
            domain=dict,
            description="Dictionary of precipitate component stoichiometry for each reaction.",
        ),
    )

    def build(self):
        """
        Callable method for block construction.
        """
        super().build()

        self.SolidPhase = Phase()
        self.component_list = self.config.precip_comp_list

        ## precipitation equilibrium reaction parameters
        # precipitation equilibrium reaction index
        self.eq_rxn_set = Set(
            initialize=list(
                set(key for key in self.config.precip_eq_rxn_logkeq_dict.keys())
            )
        )

        # stoichiometry for each precipitation equilibrium reaction
        self.eq_rxn_stoich_dict = self.config.precip_eq_rxn_stoich_dict

        # log(keq) for each equilibrium precipitation reaction
        self.eq_rxn_logkeq_dict = self.config.precip_eq_rxn_logkeq_dict

        self._state_block_class = PrecipitateStateBlock

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        obj.define_custom_properties(
            {
                "molality_precip_comp": {"method": None},
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


class _PrecipitateStateBlock(StateBlock):
    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """

        # Fix state variables
        return fix_state_vars(self)


@declare_process_block_class(
    "PrecipitateStateBlock", block_class=_PrecipitateStateBlock
)
class PrecipitateStateBlockData(StateBlockData):
    def build(self):

        self.moles_precip_comp = Var(
            self.component_list,
            units=pyunits.mol / pyunits.hour,
            initialize=1e-20,
            bounds=(1e-20, None),
        )

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        return {
            "moles_precip_comp": self.moles_precip_comp,
        }
