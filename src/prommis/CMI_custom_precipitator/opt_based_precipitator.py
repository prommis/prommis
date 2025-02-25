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
Optimization-based precipitator model.
"""

# Import Pyomo libraries
import pyomo.environ as pyo
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
from idaes.core.solvers import get_solver

import idaes.core.util.scaling as iscale

from idaes.core.scaling import CustomScalerBase


@declare_process_block_class("Precipitator")
class PrecipitatorData(UnitModelBlockData):
    """"""

    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "property_package_aqueous",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for aqueous control volume",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
"Valid values:** {
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
            default=True,
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

    def build(self):
        """Building model

        Args:
            None
        Return:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(PrecipitatorData, self).build()

        # Add Control Volumes
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
        self.add_inlet_port(block=self.cv_precipitate, name="precipitate_inlet")
        self.add_outlet_port(block=self.cv_precipitate, name="precipitate_outlet")

        prop_aq = self.config.property_package_aqueous
        prop_precip = self.config.property_package_precipitate

        # Params
        # create a set containing all reaction logkeq values
        self.merged_logkeq_dict = (
            prop_aq.eq_rxn_logkeq_dict | prop_precip.eq_rxn_logkeq_dict
        )
        # create a set containing all reactions
        self.merged_rxns = self.merged_logkeq_dict.keys()

        # make a param of all the logkeq values (precipitation and aqueous equilibrium)
        self.log_k = pyo.Param(
            self.merged_logkeq_dict.keys(),
            initialize=lambda m, key: self.merged_logkeq_dict[key],
            within=pyo.Reals,
        )

        # variables
        self.rxn_extent = pyo.Var(
            self.merged_rxns,
            initialize=1,
            doc="The extent of the reactions",
            units=pyo.units.mol / pyo.units.kg,
        )
        # log(q) for precipitation reactions
        self.log_q = pyo.Var(
            prop_precip.eq_rxn_set,
            initialize=1,
            doc="log(q) var for each reaction",
            units=pyo.units.dimensionless,
        )

        # constraints
        @self.Constraint(
            self.flowsheet().time,
            prop_aq.eq_rxn_set,
            doc="equilibrium reactions (log form) constraints (non-precipitate forming)",
        )
        def log_k_equil_rxn_eqns(blk, t, r):

            return blk.log_k[r] == sum(
                prop_aq.eq_rxn_stoich_dict[r][c]
                * pyo.log10(blk.cv_aqueous.properties_out[t].molality_aq_comp[c])
                for c in prop_aq.eq_rxn_stoich_dict[r]
            )

        @self.Constraint(
            self.flowsheet().time,
            prop_precip.eq_rxn_set,
            doc="equilibrium reactions (log form) constraints (precipitate forming)",
        )
        def log_q_precip_equil_rxn_eqns(blk, t, r):

            return blk.log_q[r] == sum(
                prop_aq.eq_rxn_stoich_dict[r][c]
                * pyo.log10(blk.cv_aqueous.properties_out[t].molality_aq_comp[c])
                for c in prop_aq.eq_rxn_stoich_dict[r]
            )

        # log(q) must be <= log(k) for precipitating reactions
        @self.Constraint(
            prop_precip.eq_rxn_set,
            doc="log(q) must be less than or equal to log(k) for precipitating reactions",
        )
        def precip_rxns_log_cons(blk, r):
            return blk.log_q[r] <= self.log_k[r]

        # mole balance on all aqueous components
        @self.Constraint(
            self.flowsheet().time,
            prop_aq.component_list,
            doc="Aqueous Components Mole balance equations.",
        )
        def aq_mole_balance_eqns(blk, t, comp):
            return blk.cv_aqueous.properties_out[t].molality_aq_comp[
                comp
            ] == blk.cv_aqueous.properties_in[t].molality_aq_comp[comp] + sum(
                prop_aq.eq_rxn_stoich_dict[r][comp] * blk.rxn_extent[r]
                for r in self.merged_rxns
            )

        # mole balance on all precipitate components
        @self.Constraint(
            self.flowsheet().time,
            prop_precip.component_list,
            doc="Precipitate Components Mole balance equations.",
        )
        def precip_mole_balance_eqns(blk, t, comp):
            return blk.cv_precipitate.properties_out[t].molality_precip_comp[
                comp
            ] == blk.cv_precipitate.properties_in[t].molality_precip_comp[comp] + sum(
                prop_precip.eq_rxn_stoich_dict[r][comp] * blk.rxn_extent[r]
                for r in prop_precip.eq_rxn_set
            )

        # objective function to minimize difference between log(k) and log(q)
        @self.Objective(
            doc="objective function minimize difference between log(k) and log(q)"
        )
        def min_logs(blk, r):
            return sum(
                (self.log_k[r] - blk.log_q[r]) ** 2 for r in prop_precip.eq_rxn_set
            )

    def solve_unit(self):
        """solves unit"""

        solver_obj = get_solver(solver="ipopt_v2", writer_config={"scale_model": True})

        results = solver_obj.solve(self, tee=True)

        return results
