#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Preliminary property package for West Kentucky No. 13 coal refuse.

Authors: Alejandro Garciadiego, Bo-Xun Wang
"""

from pyomo.common.config import ConfigValue
from pyomo.environ import Param, Var, units

from idaes.core import (
    MaterialFlowBasis,
    Phase,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    declare_process_block_class,
)
from idaes.core.util.initialization import fix_state_vars
from idaes.core.scaling import CustomScalerBase

# -----------------------------------------------------------------------------
# Precipitate solids property package


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


class PrecipitateParametersScaler(CustomScalerBase):
    """
    Scaler for the precipitate solids property package.
    """

    DEFAULT_SCALING_FACTORS = {
        "temperature": 1 / 300,
        "flow_mol_comp": 1e3,
    }

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # Scale state variables
        self.scale_variable_by_default(model.temperature, overwrite=overwrite)
        for var in model.flow_mol_comp.values():
            self.scale_variable_by_default(var, overwrite=overwrite)

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # No constraints to scale
        pass


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

        self.reaction_to_element = {
            "Sc2(C2O4)3(s)": "Sc_3+",
            "Y2(C2O4)3(s)": "Y_3+",
            "La2(C2O4)3(s)": "La_3+",
            "Ce2(C2O4)3(s)": "Ce_3+",
            "Pr2(C2O4)3(s)": "Pr_3+",
            "Nd2(C2O4)3(s)": "Nd_3+",
            "Sm2(C2O4)3(s)": "Sm_3+",
            "Gd2(C2O4)3(s)": "Gd_3+",
            "Dy2(C2O4)3(s)": "Dy_3+",
            "Al2(C2O4)3(s)": "Al_3+",
            "Ca(C2O4)(s)": "Ca_2+",
            "Fe2(C2O4)3(s)": "Fe_3+",
        }

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
    default_scaler = PrecipitateParametersScaler

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

    default_scaler = PrecipitateParametersScaler

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
            bounds=(1e-25, None),
        )

    def get_material_flow_terms(self, p, j):
        return self.flow_mol_comp[j]

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        return {
            "flow_mol_comp": self.flow_mol_comp,
            "temperature": self.temperature,
        }
