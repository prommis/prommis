#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Preliminary Precipitator Unit Model
===================================

Author: Alejandro Garciadiego

The Precipitator Unit Model represents an Equilibrium reactor unit model with fixed partition coefficients.

Configuration Arguments
-----------------------

The precipitator unit model needs an aqueous property package which includes stoichiometric values for solids being
created in the precipitator and fixed separation coefficients of the solids.

Model Structure
---------------

The Precitator unit model has hard coded stream names (``aqueous`` and ``precipitate`` respectively). The Precipitator 
model also has one inlet and two outlets named ``aqueous_inlet``, ``aqueous_outlet`` and ``precipitate_outlet`` respectively.

Additional Constraints
----------------------

The Precipitator unit adds two additional constraint to define the stochiometry and separation.

.. math:: n_{t,prec,c} = \frac{n_{t,aq_in,c} - n_{t,aq_out,c}}{S_{comp}} 

where :math:`n_{t,prec,c}` is the outlet precipitation of component c, :math:`n_{t,aq_in,c}` is the inlet of component c in 
the aqueous phase, :math:`n_{t,aq_in,c}` is the outlet of component c in the aqueous phase at time :math:`t`, divided by the
stoichiometric parameter of component c :math:`S_{comp}`

.. math:: n_{t,aq_out,c} = n_{t,aq_in,c} * (1 - split_{c})

where :math:`split_{c}` is the fixed recovery fraction of component c; this factor can be a parameter or ideally a variable
solved by a surrogate or a model equation.

"""

# Import Pyomo libraries
from pyomo.common.config import Bool, ConfigBlock, ConfigValue

import idaes.logger as idaeslog

# Import IDAES cores
from idaes.core import (
    ControlVolume0DBlock,
    UnitModelBlockData,
    declare_process_block_class,
    useDefault,
)
from idaes.core.util.config import (
    is_physical_parameter_block,
    is_reaction_parameter_block,
)

_log = idaeslog.getLogger(__name__)


@declare_process_block_class("Precipitator")
class PrecipitatorData(UnitModelBlockData):
    """ """

    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "property_package_aqueous",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for aqueous control volume",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args_aqueous",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing aqueous property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "property_package_precipitate",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for precipitate control volume",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args_precipitate",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing precipitate property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "has_equilibrium_reactions",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Equilibrium reaction construction flag",
            doc="""Indicates whether terms for equilibrium controlled reactions
should be constructed,
**default** - True.
**Valid values:** {
**True** - include equilibrium reaction terms,
**False** - exclude equilibrium reaction terms.}""",
        ),
    )
    CONFIG.declare(
        "has_phase_equilibrium",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Phase equilibrium construction flag",
            doc="""Indicates whether terms for phase equilibrium should be
constructed,
**default** = False.
**Valid values:** {
**True** - include phase equilibrium terms
**False** - exclude phase equilibrium terms.}""",
        ),
    )
    CONFIG.declare(
        "has_heat_of_reaction",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Heat of reaction term construction flag",
            doc="""Indicates whether terms for heat of reaction terms should be
constructed,
**default** - False.
**Valid values:** {
**True** - include heat of reaction terms,
**False** - exclude heat of reaction terms.}""",
        ),
    )
    CONFIG.declare(
        "reaction_package",
        ConfigValue(
            default=None,
            domain=is_reaction_parameter_block,
            description="Reaction package to use for control volume",
            doc="""Reaction parameter object used to define reaction calculations,
**default** - None.
**Valid values:** {
**None** - no reaction package,
**ReactionParameterBlock** - a ReactionParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "reaction_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing reaction packages",
            doc="""A ConfigBlock with arguments to be passed to a reaction block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see reaction package for documentation.}""",
        ),
    )

    def build(self):
        """
        Build method for precipitator unit model.
        """
        # Call UnitModel.build to setup dynamics
        super(PrecipitatorData, self).build()

        # Add Control Volume
        self.cv_aqueous = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package_aqueous,
            property_package_args=self.config.property_package_args_aqueous,
        )
        # Add inlet and outlet state blocks to control volume
        self.cv_aqueous.add_state_blocks(has_phase_equilibrium=False)

        # ---------------------------------------------------------------------
        # Add single state block for vapor phase
        tmp_dict = dict(**self.config.property_package_args_precipitate)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["defined_state"] = False
        self.cv_precipitate = (
            self.config.property_package_precipitate.build_state_block(
                self.flowsheet().time, doc="Vapor phase properties", **tmp_dict
            )
        )

        # ---------------------------------------------------------------------
        # Check flow basis is compatible
        # TODO : Could add code to convert flow bases, but not now
        t_init = self.flowsheet().time.first()
        if (
            self.cv_precipitate[t_init].get_material_flow_basis()
            != self.cv_aqueous.properties_out[t_init].get_material_flow_basis()
        ):
            raise ConfigurationError(
                f"{self.name} Solid and aqueous property packages must use the "
                f"same material flow basis."
            )

        # add ports
        self.add_inlet_port(block=self.cv_aqueous, name="aqueous_inlet")
        self.add_outlet_port(block=self.cv_aqueous, name="aqueous_outlet")
        self.add_outlet_port(block=self.cv_precipitate, name="precipitate_outlet")

        prop_aq = self.config.property_package_aqueous
        prop_s = self.config.property_package_precipitate

        @self.Constraint(self.flowsheet().time, doc="volume balance equation.")
        def vol_balance(blk, t):
            return blk.cv_aqueous.properties_out[t].flow_vol == (
                blk.cv_aqueous.properties_in[t].flow_vol
            )

        @self.Constraint(
            self.flowsheet().time,
            prop_s.component_list,
            doc="Mass balance equations precipitate.",
        )
        def precipitate_generation(blk, t, comp):
            return blk.cv_precipitate[t].flow_mol_comp[comp] == (
                (
                    blk.cv_aqueous.properties_in[t].flow_mol_comp[prop_s.react[comp]]
                    - blk.cv_aqueous.properties_out[t].flow_mol_comp[prop_s.react[comp]]
                )
                / prop_s.stoich[comp]
            )

        @self.Constraint(
            self.flowsheet().time,
            prop_aq.dissolved_elements,
            doc="Mass balance equations aqueous.",
        )
        def aqueous_depletion(blk, t, comp):
            return blk.cv_aqueous.properties_out[t].conc_mass_comp[
                comp
            ] * blk.cv_aqueous.properties_out[t].flow_vol == (
                blk.cv_aqueous.properties_in[t].conc_mass_comp[comp]
                * blk.cv_aqueous.properties_in[t].flow_vol
                * (1 - prop_aq.split[comp] / 100)
            )
