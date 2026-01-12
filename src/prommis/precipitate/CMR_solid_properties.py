#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Preliminary property package for CMR wastewater treatment.

Authors: Chenyu Wang
"""

from pyomo.common.config import ConfigValue
from pyomo.environ import Param, Var, units

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


def _config_blk_build(blk):
    blk.declare(
        "key_components",
        ConfigValue(
            default=None,
            domain=set,
            description="Set of key components",
            doc="(set) Set of key components",
        ),
    )


@declare_process_block_class("PrecipitateParameters")
class PrecipitateParametersData(PhysicalParameterBlock):
    """
    Solid phase property package for oxalate precipitation.

    Based on assay provided in:

    RESEARCH PERFORMANCE FINAL REPORT, Pilot-Scale Testing of an Integrated
    Circuit for the Extraction of Rare Earth Minerals and Elements from Coal
    and Coal Byproducts Using Advanced Separation Technologies,
    Honaker, R.Q., et al., DE-FE0027035

    Includes the following components:
    "CoFe2O4(s)", "Fe2O3(s)", "Zn4(OH)6SO4(s)",
    """

    CONFIG = PhysicalParameterBlock.CONFIG()
    _config_blk_build(CONFIG)

    def build(self):
        super().build()

        self.solid = Phase()

        comp_list = [
            "CoFe2O4(s)",
            "Fe2O3(s)",
            "Zn4(OH)6SO4(s)",
        ]

        self.component_list = comp_list

        self.react = {
            "CoFe2O4(s)": "Co",
            "Fe2O3(s)": "Fe",
            "Zn4(OH)6SO4(s)": "Zn",
        }

        self.stoich = Param(
            self.component_list,
            units=units.mol / units.mol,
            initialize={
                "CoFe2O4(s)": 1,
                "Fe2O3(s)": 2,
                "Zn4(OH)6SO4(s)": 4,
            },
        )

        self.mw = Param(
            self.component_list,
            units=units.kg / units.mol,
            initialize={
                "CoFe2O4(s)": 234.62e-3,
                "Fe2O3(s)": 159.69e-3,
                "Zn4(OH)6SO4(s)": 459.63e-3,
            },
        )

        self._state_block_class = PrecipitateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "flow_mol_comp": {"method": None},
                "temperature": {"method": None},
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


class _PrecipitateBlock(StateBlock):
    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        # Fix state variables
        fix_state_vars(self)


@declare_process_block_class("PrecipitateBlock", block_class=_PrecipitateBlock)
class PrecipitateStateBlockData(StateBlockData):
    """
    State block for solid REE oxalate.
    """

    def build(self):
        super().build()

        # State Variables
        self.temperature = Var(
            initialize=348.15,
            doc="Temperature",
            units=units.kelvin,
            bounds=(298.15, None),
        )

        self.flow_mol_comp = Var(
            self.params.component_list,
            units=units.mol / units.hour,
            initialize=1e-5,
            bounds=(1e-20, None),
        )

        iscale.set_scaling_factor(self.flow_mol_comp, 1e3)
        iscale.set_scaling_factor(self.temperature, 1e1)

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        return {
            "flow_mol_comp": self.flow_mol_comp,
            "temperature": self.temperature,
        }
