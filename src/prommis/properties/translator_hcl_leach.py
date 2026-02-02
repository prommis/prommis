#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

r"""
Translator block between sulfuric acid leaching and HCl stripping properties

========================================================================

Author: Douglas Allan

Model description
-----------------

This block takes a stream using the HCl stripping properties to the H2SO4 leaching properties.

Parameter:
----------
eps_conc_mass: Near-zero mass concentration (initialized at 1e-15 mg/L) to use for sulfate concentration.


Additional constraints
----------------------

1. eq_flow_vol_rule: Inlet and outlet volumetric flow rates are equal
2. conc_mass_comp_hcl_eqn: The concentrations of components in the HCl property package are equal
3. conc_mass_sulfates_eqn: Sulfate components have near-zero concentration, defined by eps_conc_mass

"""

from pyomo.environ import Set, Param, units as pyunits

# Import IDAES cores
from idaes.core import declare_process_block_class
from idaes.core.scaling import ConstraintScalingScheme, CustomScalerBase
from idaes.models.unit_models.translator import TranslatorData

import idaes.logger as idaeslog

__author__ = "Douglas Allan"


# Set up logger
_log = idaeslog.getLogger(__name__)


class TranslatorHClLeachScaler(CustomScalerBase):
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
        for condata in model.conc_mass_comp_hcl_eqn.values():
            self.scale_constraint_by_nominal_value(
                condata,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )
        for condata in model.conc_mass_sulfates_eqn.values():
            self.scale_constraint_by_nominal_value(
                condata,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )


@declare_process_block_class("TranslatorHClLeach")
class TranslatorHClLeachData(TranslatorData):
    """
    Translator block to go from the HCl Stripping property package,
    which does not contain sulfate species, to the full leaching
    property package, which does contain sulfate species.

    """

    default_scaler = TranslatorHClLeachScaler

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

        self.eps_conc_mass = Param(
            initialize=1e-15,
            mutable=True,
            units=pyunits.mg / pyunits.L,
            doc="Value to use for mass concentration of sulfate species "
            "in outlet stream.",
        )

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality volumetric flow equation",
        )
        def flow_vol_eqn(blk, t):
            return blk.properties_out[t].flow_vol == blk.properties_in[t].flow_vol

        self.HCl_components = Set(
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
                "Cl",
            ]
        )
        self.sulfate_components = Set(initialize=["HSO4", "SO4"])

        @self.Constraint(
            self.flowsheet().time,
            self.HCl_components,
            doc="Defines mass concentration of components from HCl properties",
        )
        def conc_mass_comp_hcl_eqn(blk, t, i):
            return (
                blk.properties_out[t].conc_mass_comp[i]
                == blk.properties_in[t].conc_mass_comp[i]
            )

        @self.Constraint(
            self.flowsheet().time,
            self.sulfate_components,
            doc="Defines mass concentration for sulfate species",
        )
        def conc_mass_sulfates_eqn(blk, t, i):
            return blk.properties_out[t].conc_mass_comp[i] == blk.eps_conc_mass
