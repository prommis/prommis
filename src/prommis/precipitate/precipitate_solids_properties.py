#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Preliminary property package for West Kentucky No. 13 coal refuse.

Authors: Alejandro Garciadiego
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
    * Rare Earth Oxalates: "Al2(C2O4)3(s)", "Fe2(C2O4)3(s)", "Sc2(C2O4)3(s)",
    "Y2(C2O4)3(s)", "La2(C2O4)3(s)", "Ce2(C2O4)3(s)", "Pr2(C2O4)3(s)",
    "Nd2(C2O4)3(s)", "Sm2(C2O4)3(s)", "Gd2(C2O4)3(s)", "Dy2(C2O4)3(s)"

    """

    CONFIG = PhysicalParameterBlock.CONFIG()
    _config_blk_build(CONFIG)

    def build(self):
        super().build()

        self.solid = Phase()

        comp_list = [
            "Al2(C2O4)3(s)",
            "Ca(C2O4)(s)",
            "Fe2(C2O4)3(s)",
            "Sc2(C2O4)3(s)",
            "Y2(C2O4)3(s)",
            "La2(C2O4)3(s)",
            "Ce2(C2O4)3(s)",
            "Pr2(C2O4)3(s)",
            "Nd2(C2O4)3(s)",
            "Sm2(C2O4)3(s)",
            "Gd2(C2O4)3(s)",
            "Dy2(C2O4)3(s)",
        ]

        self.component_list = comp_list

        self.react = {
            "Sc2(C2O4)3(s)": "Sc",
            "Y2(C2O4)3(s)": "Y",
            "La2(C2O4)3(s)": "La",
            "Ce2(C2O4)3(s)": "Ce",
            "Pr2(C2O4)3(s)": "Pr",
            "Nd2(C2O4)3(s)": "Nd",
            "Sm2(C2O4)3(s)": "Sm",
            "Gd2(C2O4)3(s)": "Gd",
            "Dy2(C2O4)3(s)": "Dy",
            "Al2(C2O4)3(s)": "Al",
            "Ca(C2O4)(s)": "Ca",
            "Fe2(C2O4)3(s)": "Fe",
        }

        self.stoich = Param(
            self.component_list,
            units=units.mol / units.mol,
            initialize={
                "Sc2(C2O4)3(s)": 2,
                "Y2(C2O4)3(s)": 2,
                "La2(C2O4)3(s)": 2,
                "Ce2(C2O4)3(s)": 2,
                "Pr2(C2O4)3(s)": 2,
                "Nd2(C2O4)3(s)": 2,
                "Sm2(C2O4)3(s)": 2,
                "Gd2(C2O4)3(s)": 2,
                "Dy2(C2O4)3(s)": 2,
                "Al2(C2O4)3(s)": 2,
                "Ca(C2O4)(s)": 1,
                "Fe2(C2O4)3(s)": 2,
            },
        )

        self.mw = Param(
            self.component_list,
            units=units.kg / units.mol,
            initialize={
                "Sc2(C2O4)3(s)": 354 * 1e-3,
                "Y2(C2O4)3(s)": 441.87 * 1e-3,
                "La2(C2O4)3(s)": 541.87 * 1e-3,
                "Ce2(C2O4)3(s)": 544.286 * 1e-3,
                "Pr2(C2O4)3(s)": 545.87 * 1e-3,
                "Nd2(C2O4)3(s)": 552.54 * 1e-3,
                "Sm2(C2O4)3(s)": 564.77 * 1e-3,
                "Gd2(C2O4)3(s)": 578.56 * 1e-3,
                "Dy2(C2O4)3(s)": 769.21 * 1e-3,
                "Al2(C2O4)3(s)": 318.02 * 1e-3,
                "Ca(C2O4)(s)": 128.097 * 1e-3,
                "Fe2(C2O4)3(s)": 143.86 * 1e-3,
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
