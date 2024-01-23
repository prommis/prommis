#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Basic precipitator model
"""


# Import Pyomo libraries
from pyomo.common.config import Bool, ConfigBlock, ConfigValue
from idaes.core.util.math import smooth_max

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
        """Building model

        Args:
            None
        Returns:
            None
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
        self.cv_precipitate = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package_precipitate,
            property_package_args=self.config.property_package_args_precipitate,
        )
        # Add inlet and outlet state blocks to control volume
        self.cv_aqueous.add_state_blocks(has_phase_equilibrium=False)
        self.cv_precipitate.add_state_blocks(has_phase_equilibrium=False)
        # add ports
        self.add_inlet_port(block=self.cv_aqueous, name="aqueous_inlet")
        self.add_outlet_port(block=self.cv_aqueous, name="aqueous_outlet")
        self.add_outlet_port(block=self.cv_precipitate, name="precipitate_outlet")

        prop_aq = self.config.property_package_aqueous
        prop_s = self.config.property_package_precipitate

        @self.Constraint(self.flowsheet().time, doc="Mass balance equations.")
        def vol_balance(blk, t):
            return blk.cv_aqueous.properties_out[t].flow_vol == (
                blk.cv_aqueous.properties_in[t].flow_vol
            )

        @self.Constraint(
            self.flowsheet().time, prop_s.component_list, doc="Mass balance equations."
        )
        def in_zero(blk, t, comp):
            return self.cv_precipitate.properties_in[t].flow_mol_comp[comp] == 0

        @self.Constraint(
            self.flowsheet().time, prop_s.component_list, doc="Mass balance equations."
        )
        def generation(blk, t, comp):
            return blk.cv_precipitate.properties_out[t].flow_mol_comp[
                comp
            ] == blk.cv_precipitate.properties_in[t].flow_mol_comp[comp] + (
                (
                    blk.cv_aqueous.properties_in[t].flow_mol_comp[prop_s.react[comp]]
                    - blk.cv_aqueous.properties_out[t].flow_mol_comp[prop_s.react[comp]]
                )
                / prop_s.stoich[comp]
            )

        @self.Constraint(
            self.flowsheet().time,
            prop_aq.dissolved_elements,
            doc="Mass balance equations.",
        )
        def mass_balance(blk, t, comp):
            return blk.cv_aqueous.properties_out[t].conc_mass_comp[
                comp
            ] * blk.cv_aqueous.properties_out[t].flow_vol == (
                blk.cv_aqueous.properties_in[t].conc_mass_comp[comp]
                * blk.cv_aqueous.properties_in[t].flow_vol
                * (1 - smooth_max(0, prop_aq.split[comp], eps=1e-8) / 100)
            )

        @self.Constraint(self.flowsheet().time)
        def temperature_s_eqn(blk, t):
            return (
                blk.cv_precipitate.properties_out[t].temperature
                == blk.cv_precipitate.properties_in[t].temperature
            )
