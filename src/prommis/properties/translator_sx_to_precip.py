#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

r"""
Translator block to convert from mixed acid properties without oxalates to
mixed acid properties with oxalates

========================================================================

Author: Marcus Holly

Model description
-----------------

This block takes a stream using the mixed acid properties without oxalates to mixed acid properties with oxalates.

Parameter:
----------
oxalic_acid_feed: Concentration of oxalic acid fed into the precipitator inlet.


Additional constraints
----------------------

1. eq_flow_vol_rule: Inlet and outlet volumetric flow rates are equal
2. conc_mass_comp_shared_eqn: The concentrations of the shared components are equal
3. conc_mass_comp_oxalate_eqn: Oxalate components have near-zero concentration, defined by eps_conc_mass

"""

from pyomo.common.config import ConfigValue, In
from pyomo.environ import Set, Param, value, units as pyunits

# Import IDAES cores
from idaes.core import declare_process_block_class
from idaes.core.scaling import ConstraintScalingScheme, CustomScalerBase
from idaes.models.unit_models.translator import TranslatorData

import idaes.logger as idaeslog

__author__ = "Marcus Holly"


# Set up logger
_log = idaeslog.getLogger(__name__)


class TranslatorSXPrecipScaler(CustomScalerBase):
    """
    Scaler for blocks with a single state (Feed, Product, StateJunction)
    """

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        self.call_submodel_scaler_method(
            submodel=model.properties_in,
            submodel_scalers=submodel_scalers,
            method="variable_scaling_routine",
            overwrite=overwrite,
        )
        self.call_submodel_scaler_method(
            submodel=model.properties_out,
            submodel_scalers=submodel_scalers,
            method="variable_scaling_routine",
            overwrite=overwrite,
        )

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        self.call_submodel_scaler_method(
            submodel=model.properties_in,
            submodel_scalers=submodel_scalers,
            method="constraint_scaling_routine",
            overwrite=overwrite,
        )
        self.call_submodel_scaler_method(
            submodel=model.properties_out,
            submodel_scalers=submodel_scalers,
            method="constraint_scaling_routine",
            overwrite=overwrite,
        )
        for condata in model.flow_vol_eqn.values():
            self.scale_constraint_by_nominal_value(
                condata,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )
        for condata in model.conc_mass_comp_shared_eqn.values():
            self.scale_constraint_by_nominal_value(
                condata,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )
        for condata in model.conc_mass_comp_oxalate_eqn.values():
            self.scale_constraint_by_nominal_value(
                condata,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )


@declare_process_block_class("TranslatorSXPrecip")
class TranslatorSXPrecipData(TranslatorData):
    """
    Translator block to go from the mixed acid property package
    that does not contain oxalate species to a mixed acid
    property package that contains oxalate species.

    """

    CONFIG = TranslatorData.CONFIG()

    # TODO: Decide whether this option should remain
    del CONFIG["outlet_state_defined"]
    del CONFIG["has_phase_equilibrium"]

    CONFIG.declare(
        "outlet_state_defined",
        ConfigValue(
            default=True,
            domain=In([True]),
        ),
    )
    CONFIG.declare(
        "has_phase_equilibrium",
        ConfigValue(
            default=False,
            domain=In([False]),
        ),
    )

    default_scaler = TranslatorSXPrecipScaler

    def fix_initialization_states(self):
        self.properties_in.fix_initialization_states()
        # Need to temporarily fix H2C2O4 during initialization, otherwise there will be 1 DOF
        # for t in self.flowsheet().time:
        #     self.properties_out[t].conc_mass_comp["H2C2O4"].fix(
        #         value(self.eps_conc_mass)
        #     )

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

        self.oxalic_acid_feed = Param(
            initialize=6400,
            mutable=True,
            units=pyunits.mg / pyunits.L,
            doc="Value to use for mass concentration of oxalate species "
            "in outlet stream.",
        )

        self.eps = Param(
            initialize=1e-9,
            mutable=True,
            units=pyunits.mg / pyunits.L,
            doc="Value to use for mass concentration of near-zero oxalate species "
            "in outlet stream.",
        )

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality volumetric flow equation",
        )
        def flow_vol_eqn(blk, t):
            return blk.properties_out[t].flow_vol == blk.properties_in[t].flow_vol

        self.shared_components = Set(
            initialize=[
                "Al_3+",
                "Ca_2+",
                "Fe_3+",
                "Sc_3+",
                "Y_3+",
                "La_3+",
                "Ce_3+",
                "Pr_3+",
                "Nd_3+",
                "Sm_3+",
                "Gd_3+",
                "Dy_3+",
                "H2O",
                "H_+",
                "Cl_-",
            ]
        )
        @self.Constraint(
            self.flowsheet().time,
            self.shared_components,
            doc="Defines mass concentration for the shared components",
        )
        def conc_mass_comp_shared_eqn(blk, t, i):
            return (
                blk.properties_out[t].conc_mass_comp[i]
                == blk.properties_in[t].conc_mass_comp[i]
            )

        #TODO: Consider how to make this oxalic concentration visible/editable in the flowsheet
        # and whether or not the concentrations of C2O4_2- and HC2O4_- must be specified

        @self.Constraint(
            self.flowsheet().time,
            doc="Defines mass concentration for the oxalate components",
        )
        def conc_mass_comp_oxalate_eqn(blk, t):
            return (
                blk.properties_out[t].conc_mass_comp["H2C2O4"]
                == blk.oxalic_acid_feed
            )

        self.near_zero_components = Set(
            initialize=[
                "C2O4_2-",
                "HC2O4_-",
            ]
        )
        @self.Constraint(
            self.flowsheet().time,
            self.near_zero_components,
            doc="Defines mass concentration for the oxalate components",
        )
        def conc_mass_comp_oxalate_species_eqn(blk, t, i):
            return (
                blk.properties_out[t].conc_mass_comp[i]
                == blk.eps
            )

        # @self.Constraint(
        #     self.flowsheet().time,
        #     doc="Equality temperature equation",
        # )
        # def eq_temperature_rule(blk, t):
        #     return blk.properties_out[t].temperature == blk.properties_in[t].temperature

        # @self.Constraint(
        #     self.flowsheet().time,
        #     doc="Equality pressure equation",
        # )
        # def eq_pressure_rule(blk, t):
        #     return blk.properties_out[t].pressure == blk.properties_in[t].pressure
