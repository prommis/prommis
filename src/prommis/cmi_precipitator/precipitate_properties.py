#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Precipitate property package for optimization-based precipitator model.

Authors: Chris Laliwala
"""

from pyomo.common.config import ConfigValue
from pyomo.environ import Set, Var
from pyomo.environ import units as pyunits

from idaes.core import (
    MaterialFlowBasis,
    Phase,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    declare_process_block_class,
)
from idaes.core.util.initialization import fix_state_vars


@declare_process_block_class("PrecipitateParameter")
class PrecipitateParameterData(PhysicalParameterBlock):
    """
    Property package for precipitate species.

    This property package requires that the user pass in a list of precipitate components (precipitate_comp_list),
    a dictionary of the precipitate-forming equilibrium reaction constants (logkeq_dict), and a dictionary containing
    the precipitate stoichiometry for each reaction (stoich_dict).
    """

    CONFIG = PhysicalParameterBlock.CONFIG()
    CONFIG.declare(
        "precipitate_comp_list",
        ConfigValue(
            domain=list, description="List of precipitate components in process."
        ),
    )
    CONFIG.declare(
        "logkeq_dict",
        ConfigValue(
            domain=dict,
            description="Dictionary of precipitation-forming equilibrium reaction constants",
        ),
    )
    CONFIG.declare(
        "stoich_dict",
        ConfigValue(
            domain=dict,
            description="Dictionary of precipitate stoichiometry for each reaction.",
        ),
    )

    def build(self):
        """
        Callable method for block construction.
        """
        super().build()

        self.SolidPhase = Phase()
        self.component_list = self.config.precipitate_comp_list

        ## precipitation equilibrium reaction parameters
        # precipitation equilibrium reaction index
        self.rxn_set = Set(
            initialize=list(set(key for key in self.config.logkeq_dict.keys()))
        )

        # stoichiometry for each precipitation equilibrium reaction
        self.stoich_dict = self.config.stoich_dict

        # log(keq) for each equilibrium precipitation reaction
        self.logkeq_dict = self.config.logkeq_dict

        self._state_block_class = PrecipitateStateBlock

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        obj.define_custom_properties(
            {
                "moles_precipitate_comp": {"method": None},
            }
        )
        obj.add_default_units(
            {
                "time": pyunits.s,
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
    """
    State block for the precipitate species.
    """

    def build(self):
        super().build()

        self.moles_precipitate_comp = Var(
            self.component_list,
            units=pyunits.mol / pyunits.s,
            initialize=1e-20,
            bounds=(1e-20, None),
        )

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        return {
            "moles_precipitate_comp": self.moles_precipitate_comp,
        }
