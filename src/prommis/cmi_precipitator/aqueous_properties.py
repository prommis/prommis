#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
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


@declare_process_block_class("AqueousParameter")
class AqueousParameterData(PhysicalParameterBlock):
    """
    Property package for aqueous species.

    This property package requires that the user pass in a list of aqueous components (aqueous_comp_list),
    a dictionary of the aqueous equilibrium reaction constants (logkeq_dict), and a dictionary containing
    the aqueous component stoichiometry for each reaction (stoich_dict).
    """

    CONFIG = PhysicalParameterBlock.CONFIG()
    CONFIG.declare(
        "aqueous_comp_list",
        ConfigValue(domain=list, description="List of aqueous components in process"),
    )
    CONFIG.declare(
        "logkeq_dict",
        ConfigValue(
            domain=dict, description="Dictionary of equilibrium reaction constants"
        ),
    )
    CONFIG.declare(
        "stoich_dict",
        ConfigValue(
            domain=dict, description="Dictionary of equilibrium reaction stoichiometry"
        ),
    )

    def build(self):
        """
        Calllable method for Block construction.
        """
        super().build()

        self.AqueousPhase = Phase()
        self.component_list = self.config.aqueous_comp_list

        ## equilibrium reaction parameters
        # equilibrium reaction index
        self.rxn_set = Set(
            initialize=list(set(key for key in self.config.logkeq_dict.keys()))
        )

        # stoichiometry for each equilibrium reaction
        self.stoich_dict = self.config.stoich_dict

        # log(keq) for each equilibrium reaction
        self.logkeq_dict = self.config.logkeq_dict

        self._state_block_class = AqueousStateBlock

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        obj.define_custom_properties(
            {
                "molality_aqueous_comp": {"method": None},
                "flow_vol": {"method": None},
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
    """
    State block for the aqueous species.
    """

    def build(self):
        super().build()

        self.flow_vol = Var(
            units=pyunits.kg / pyunits.s,
            initialize=1,
            bounds=(1e-5, None),
        )

        self.molality_aqueous_comp = Var(
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
            "molality_aqueous_comp": self.molality_aqueous_comp,
        }
