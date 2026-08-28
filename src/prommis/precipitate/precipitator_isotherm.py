#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Oxalate Precipitator Unit Model
===================================

Author: Alejandro Garciadiego, Bo-Xun Wang

The Precipitator Unit Model represents an Equilibrium reactor unit model with the equilibrium equation derived from Minteq data.

Configuration Arguments
-----------------------

The precipitator unit model needs the liquid, solid, and reaction property packages which include stoichiometric values for solids being
created in the precipitator and the parameters used in the equilibrium equation.

While the liquid, solid, and reaction property packages are included in this model, users should be able to modify or replace those
property packages (e.g., by changing the values of equilibrium parameters) to customize the model for specific applications.

- ``liquid_phase``: Configuration dictionary for the liquid (aqueous) phase, described below.
- ``solid_phase``: Configuration dictionary for the solid (precipitate) phase, described below.
- ``reaction_package``: Heterogeneous reaction package to use for precipitation. This package
provides the stoichiometric and equilibrium parameters (e.g. E_D and N_D) used to
calculate conversion for each precipitating species.
- ``reaction_package_args``: Dict of arguments to be passed to the heterogeneous reaction
package when it is constructed.
- ``number_of_tanks``: Number of tanks (finite elements) to use when constructing the internal
``MSContactor`` model. Default is 1.

Each of ``liquid_phase`` and ``solid_phase`` accepts the following sub-arguments:

- ``property_package``: Property package to use for the given phase. The default package from the parent model or flowsheet is used.
- ``property_package_args``: Dict of arguments to use when constructing the property package for the given phase.
- ``has_energy_balance``: Boolean indicating whether to include an energy balance for the
given phase. Must be ``False``, as the Precipitator does not support energy balances.
- ``has_pressure_balance``: Boolean indicating whether to include a pressure balance for the
given phase. Must be ``False``, as the Precipitator does not support pressure balances.

Model Structure
---------------

The Precipitator unit model has hard coded stream names (``aqueous`` and ``precipitate`` respectively). The Precipitator
model also has one inlet and two outlets named ``aqueous_inlet``, ``aqueous_outlet`` and ``precipitate_outlet`` respectively.

Additional Constraints
----------------------

The Precipitator unit adds one additional constraint to define the conversion.

.. math:: Conversion_{c} = \exp(-\frac{\epsilon}{Oxalic Acid Dosage}^{n_{DA}})

where :math:`Conversion` is the conversion of component c, :math:`\epsilon` and :math:`n_{DA}` are the parameters
estimated based on Minteq data, :math:`Oxalic Acid Dosage` is the amount of oxalic acid added into the precipitator.


"""
from pyomo.common.collections import ComponentMap
from pyomo.common.config import ConfigDict, ConfigValue, In

from pyomo.environ import (
    Param,
    Var,
    Block,
    value,
    log,
    NonNegativeReals,
    units as pyunits,
)

from pyomo.network import Port

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
)

from idaes.core import (
    UnitModelBlockData,
    declare_process_block_class,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block

from idaes.models.unit_models.mscontactor import MSContactor
from idaes.core.initialization import ModularInitializerBase
from idaes.core.scaling import CustomScalerBase, ConstraintScalingScheme
import math

class OxalatePrecipitatorScaler(CustomScalerBase):
    """
    Scaler for the Oxalate Precipitator unit model.
    """

    DEFAULT_SCALING_FACTORS = {
        "volume": 1e-3,
        "hydraulic_retention_time": 1e0,
        "conversion": 1e1,
        "heterogeneous_reaction_extent": 1e3,
        "solid_heterogeneous_reactions_generation": 1e3,
        "liquid_heterogeneous_reactions_generation": 1e3,
    }
 
    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: ComponentMap = None
    ):
        """
        Variable scaling routine for the Oxalate Precipitator.
 
        Args:
            model: instance of OxalatePrecipitator to be scaled
            overwrite: whether to overwrite existing scaling factors
            submodel_scalers: ComponentMap of Scalers to use for sub-models,
                keyed by submodel local name
 
        Returns:
            None
        """

        self.call_submodel_scaler_method(
            submodel=model.mscontactor,
            submodel_scalers=submodel_scalers,
            method="variable_scaling_routine",
            overwrite=overwrite,
        )
 
        for t in model.flowsheet().time:
            self.scale_variable_by_default(model.volume[t], overwrite=overwrite)
            self.scale_variable_by_default(
                model.hydraulic_retention_time[t], overwrite=overwrite
            )

        for r in model.config.reaction_package.reaction_idx:
            self.scale_variable_by_default(model.conversion[r], overwrite=overwrite)

        if hasattr(model.mscontactor, "heterogeneous_reaction_extent"):
            for v in model.mscontactor.heterogeneous_reaction_extent.values():
                self.scale_variable_by_default(v, overwrite=overwrite)

        if hasattr(model.mscontactor, "solid_heterogeneous_reactions_generation"):
            for v in model.mscontactor.solid_heterogeneous_reactions_generation.values():
                self.scale_variable_by_default(v, overwrite=overwrite)

        if hasattr(model.mscontactor, "liquid_heterogeneous_reactions_generation"):
            for v in model.mscontactor.liquid_heterogeneous_reactions_generation.values():
                self.scale_variable_by_default(v, overwrite=overwrite)
 
    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: ComponentMap = None
    ):
        """
        Constraint scaling routine for the Oxalate Precipitator.
 
        Args:
            model: instance of OxalatePrecipitator to be scaled
            overwrite: whether to overwrite existing scaling factors
            submodel_scalers: ComponentMap of Scalers to use for sub-models,
                keyed by submodel local name
 
        Returns:
            None
        """
        self.call_submodel_scaler_method(
            submodel=model.mscontactor,
            submodel_scalers=submodel_scalers,
            method="constraint_scaling_routine",
            overwrite=overwrite,
        )
 
        for condata in model.eq_hydraulic_retention.values():
            self.scale_constraint_by_nominal_value(
                condata,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )

        for condata in model.heterogeneous_reaction_extent_constraint.values():
            self.scale_constraint_by_nominal_value(
                condata,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )

        for condata in model.temp_constraint.values():
            self.scale_constraint_by_nominal_value(
                condata,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )
 
        for condata in model.liq_temp_constraint.values():
            self.scale_constraint_by_nominal_value(
                condata,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )
 
        for condata in model.press_constraint.values():
            self.scale_constraint_by_nominal_value(
                condata,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )
 
        for condata in model.init_solid_constraint.values():
            self.scale_constraint_by_nominal_value(
                condata,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )

# -----------------------------------------------------------------------------
# Precipitator unit model
class OxalatePrecipitatorInitializer(ModularInitializerBase):
    """
    This is a general purpose Initializer for the Oxalate Precipitator unit model.

    This routine calls the initializer for the internal MSContactor model.

    """

    CONFIG = ModularInitializerBase.CONFIG()

    CONFIG.declare(
        "ssc_solver_options",
        ConfigDict(
            implicit=True,
            description="Dict of arguments for solver calls by ssc_solver",
        ),
    )
    CONFIG.declare(
        "calculate_variable_options",
        ConfigDict(
            implicit=True,
            description="Dict of options to pass to 1x1 block solver",
            doc="Dict of options to pass to calc_var_kwds argument in "
            "scc_solver method.",
        ),
    )

    def initialize_main_model(
        self,
        model: Block,
    ):
        """
        Initialization routine for MSContactor Blocks.

        Args:
            model: model to be initialized

        Returns:
            None
        """
        # Initialize MSContactor
        model.mscontactor.heterogeneous_reaction_extent.fix()

        msc_init = model.mscontactor.default_initializer(
            ssc_solver_options=self.config.ssc_solver_options,
            calculate_variable_options=self.config.calculate_variable_options,
        )

        msc_init.initialize(model.mscontactor)

        model.mscontactor.heterogeneous_reaction_extent.unfix()

        for t in model.flowsheet().time:
            dosage = value(model.oxalic_acid_dosage[t])
            for r in model.config.reaction_package.reaction_idx:
                if r == "Ca(C2O4)(s)":
                    model.conversion[r].set_value(value(model.min_conversion))
                else:
                    E_D = value(model.config.reaction_package.E_D[r])
                    N_D = value(model.config.reaction_package.N_D[r])
                    if dosage > 0:
                        exponent = -(E_D ** N_D) / (dosage ** N_D)
                        exponent = max(min(exponent, 0.0), -700.0)
                        model.conversion[r].set_value(
                            min(max(math.exp(exponent), 1e-20), 0.999999)
                        )

        solver = self._get_solver()
        results = solver.solve(model, tee=True)

        return results


StreamCONFIG = ConfigDict()
StreamCONFIG.declare(
    "property_package",
    ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for given stream",
        doc="""Property parameter object used to define property calculations for given stream,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
    ),
)
StreamCONFIG.declare(
    "property_package_args",
    ConfigDict(
        implicit=True,
        description="Dict of arguments to use for constructing property package",
        doc="""A ConfigDict with arguments to be passed to property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
    ),
)
StreamCONFIG.declare(
    "has_energy_balance",
    ConfigValue(
        default=False,
        domain=In([False]),
        doc="Bool indicating whether to include energy balance for stream. Must be false",
    ),
)
StreamCONFIG.declare(
    "has_pressure_balance",
    ConfigValue(
        default=False,
        domain=In([False]),
        doc="Bool indicating whether to include pressure balance for stream. Must be false",
    ),
)


@declare_process_block_class("OxalatePrecipitator")
class OxalatePrecipitatorData(UnitModelBlockData):
    """
    Oxalate Precipitator Unit Model Class
    """

    # Set default initializer
    default_initializer = OxalatePrecipitatorInitializer
    default_scaler = OxalatePrecipitatorScaler

    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "liquid_phase",
        StreamCONFIG(
            description="Liquid phase properties",
        ),
    )
    CONFIG.declare(
        "solid_phase",
        StreamCONFIG(
            description="Solid phase properties",
        ),
    )
    CONFIG.declare(
        "reaction_package",
        ConfigValue(
            # TODO: Add a domain validator for this
            description="Heterogeneous reaction package for precipitation.",
        ),
    )
    CONFIG.declare(
        "reaction_package_args",
        ConfigValue(
            default=None,
            domain=dict,
            description="Arguments for heterogeneous reaction package for precipitation.",
        ),
    )
    CONFIG.declare(
        "number_of_tanks",
        ConfigValue(
            default=1, domain=int, description="Number of tanks in precipitation"
        ),
    )

    def build(self):
        """
        Build method for OxalatePrecipitator unit model.
        """
        super().build()

        self.mscontactor = MSContactor(
            number_of_finite_elements=self.config.number_of_tanks,
            streams={
                "liquid": {
                    "property_package": self.config.liquid_phase.property_package,
                    "property_package_args": self.config.liquid_phase.property_package_args,
                    "has_energy_balance": self.config.liquid_phase.has_energy_balance,
                    "has_pressure_balance": self.config.liquid_phase.has_pressure_balance,
                },
                "solid": {
                    "property_package": self.config.solid_phase.property_package,
                    "property_package_args": self.config.solid_phase.property_package_args,
                    "has_energy_balance": self.config.solid_phase.has_energy_balance,
                    "has_pressure_balance": self.config.solid_phase.has_pressure_balance,
                },
            },
            heterogeneous_reactions=self.config.reaction_package,
            heterogeneous_reactions_args=self.config.reaction_package_args,
        )

        self.hydraulic_retention_time = Var(
            self.flowsheet().time,
            initialize=2,
            domain=NonNegativeReals,
            units=pyunits.h,
            doc="Hydraulic retention time",
        )

        self.volume = Var(
            self.flowsheet().time,
            initialize=1800,
            domain=NonNegativeReals,
            units=pyunits.l,
            doc="Volume of precipitator",
        )

        self.conversion = Var(
            self.config.reaction_package.reaction_idx,
            initialize=0.5,
            units=pyunits.dimensionless,
            bounds=(1e-20, 0.999999),
            doc="Conversion of each precipitation species",
        )

        self.min_conversion = Param(
            initialize=1e-6,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Minimum conversion for Ca(C2O4)(s)",
        )

        # Create unit level Ports
        self.aqueous_inlet = Port(extends=self.mscontactor.liquid_inlet)
        self.aqueous_outlet = Port(extends=self.mscontactor.liquid_outlet)
        self.precipitate_outlet = Port(extends=self.mscontactor.solid_outlet)

        @self.Expression(
            self.flowsheet().time,
            doc="Oxalic acid dosage",
        )
        def oxalic_acid_dosage(blk, t):
            return (blk.aqueous_inlet.conc_mass_comp[0, "H2C2O4"] +
                    blk.aqueous_inlet.conc_mass_comp[0, "HC2O4_-"] +
                    blk.aqueous_inlet.conc_mass_comp[0, "C2O4_2-"]) / (
                1000 * pyunits.mg / pyunits.l
            )

        @self.Constraint(self.flowsheet().time, doc="Hydraulic retention time equation")
        def eq_hydraulic_retention(blk, t):
            return blk.hydraulic_retention_time[t] * pyunits.convert(
                blk.aqueous_inlet.flow_vol[t], to_units=pyunits.m**3 / pyunits.hour
            ) == pyunits.convert(blk.volume[t], to_units=pyunits.m**3)

        @self.Constraint(
            self.flowsheet().time,
            self.mscontactor.elements,
            self.config.reaction_package.reaction_idx,
            doc="Reaction extent constraint",
        )
        def heterogeneous_reaction_extent_constraint(blk, t, s, r):
            return blk.mscontactor.heterogeneous_reaction_extent[t, s, r] == (
                blk.mscontactor.heterogeneous_reactions[t, s].reaction_rate[r]
                - (
                    blk.conversion[r]
                    * blk.mscontactor.liquid_inlet_state[t].flow_mol_comp[
                        blk.mscontactor.config.streams.solid.property_package.reaction_to_element[r]
                    ]
                )
            )

        @self.Constraint(
            self.flowsheet().time,
            self.mscontactor.elements,
            self.config.reaction_package.reaction_idx,
            doc="conversion constraint",
        )
        def conversion_constraint(blk, t, s, r):
            if r == "Ca(C2O4)(s)":
                return blk.conversion[r] == blk.min_conversion
            else:
                return log(blk.conversion[r]) == (
                    -(
                        (blk.config.reaction_package.E_D[r])
                        ** blk.config.reaction_package.N_D[r]
                    )
                ) / (blk.oxalic_acid_dosage[t] ** blk.config.reaction_package.N_D[r])

        @self.Constraint(self.flowsheet().time, doc="temperature equation")
        def temp_constraint(blk, t):
            return (
                blk.mscontactor.solid_inlet_state[t].temperature
                == blk.mscontactor.solid_outlet.temperature[t]
            )

        @self.Constraint(self.flowsheet().time, doc="liquid temperature equation")
        def liq_temp_constraint(blk, t):
            return blk.aqueous_inlet.temperature[t] == blk.aqueous_outlet.temperature[t]

        @self.Constraint(self.flowsheet().time, doc="pressure equation")
        def press_constraint(blk, t):
            return blk.aqueous_inlet.pressure[t] == blk.aqueous_outlet.pressure[t]

        @self.Constraint(
            self.flowsheet().time,
            self.config.solid_phase.property_package.component_list,
            doc="Initial solids",
        )
        def init_solid_constraint(blk, t, r):
            return (
                blk.mscontactor.solid_inlet_state[t].flow_mol_comp[r]
                == 1e-9 * pyunits.mole / pyunits.hour
            )

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        expr_dict = {}
        param_dict = {}
        var_dict["Unit Volume"] = self.volume[time_point]
        var_dict["Hydraulic Retention Time"] = self.hydraulic_retention_time[time_point]
        expr_dict["Oxalic Acid Dosage"] = self.oxalic_acid_dosage[time_point]
        return {"vars": var_dict, "params": param_dict, "exprs": expr_dict}
