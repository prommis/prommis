#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Optimization-based Precipitator for the Critical Materials Innovation Hub (CMI) Process
=======================================================================================

Author: Chris Laliwala

This precipitator makes use of equilibrium constants and solubility constants to predict the final concentrations of aqueous species,
and the amount of precipitates formed at equilibrium.

Configuration Arguments
-----------------------

The model requires the user to define the aqueous species and precipitates present in the system. Additionally,
equilibrium constants for aqueous reactions, and solubility constants for precipitation/dissolution reactions are needed.
Reaction stoichiometry must also be provided by the user. The aqueous components, equilibrium constants, and a dictionary
defining the aqueous reaction stoichiometry is stored in the aqeuous property package. The precipitates, solubility constants,
and a dictionary defining the precipitation/dissolution reaction stoichiometry is stored in the precipitate property package.

Model Structure
---------------

This unit model has an inlet for aqueous ('aqueous_inlet') and precipitates ('precipitate_inlet') entering the unit, and an outlet for
aqueous ('aqueous_outlet') and precipitates ('precipitate_outlet') leaving the unit.


Additional Model Information
----------------------------

This precipitator model seeks to model aqueous systems involving precipitation and dissolution reactions as chemical equilibrium problems.
The approach taken here is to solve a system of nonlinear equations involving equilibrium constants (the law of mass action approach, LMA) [1].
Instead of utilizing saturation indices heuristics commonly used by LMA software [1], this model formulates an optimization problem where the
objective function is to minimize the square difference between the ion product, :math:`Q_{r,sp}`, defined over the actual concentration in solution,
and the solubility constant, :math:`K_{r,sp}`, defined over the equilibrium concentration in solution,

.. math:: z = \sum_{r \in N_{rxn,sp}} ( log(K_{r,sp}) - log(Q_{r,sp}) )^2 \quad (1)

where :math:`N_{rxn,sp}` is the set of precipitation/dissolution reactions. By adding in the following constraints, this objective function allows the
identification of the species that should precipitate (i.e. :math:`Q_{r,sp} = K_{r,sp}`) and those that should not (i.e. :math:`Q_{r,sp} \leq K_{r,sp}`):

.. math:: log(Q_{r,sp}) = \sum_{i \in I_{r, products}} \alpha_{i,r} log(C_i^f) \forall r \in N_{rxn,sp} \quad (2)

.. math:: log(Q_{r,sp}) \leq log(K_{r,sp}) \forall r \in N_{rxn,sp} \quad (3)

where :math:`I_{r,products}` is the set of products formed in reaction :math:`r`. Eq. (2) computes the ionic product based on the actual concentrations in
the solution; its upper bound is defined by the solubility product (Eq. 3).

Equilibrium conditions are imposed in all reactions that do not involve solids' formation as shown in Eq. (4):

.. math:: log(K_{r,aq}) = \sum_{i \in I_{r,products}} \alpha_{i,r} log(C_i^f) - \sum_{i \in I_{r, reactants}} \alpha_{i,r} log(C_i^f) \forall r \in N_{rxn,aq} \quad (4)

where :math:`\alpha_{i,r}` is the stoichiometric coefficient of species :math:`i` in reaction :math:`r`, :math:`C_i^f` is the equilibrium concentration of
species :math:`i`, and :math:`N_{rxn,aq}` is the set of all aqueous reactions. The logartihmic form(s) in Eq. (2-4) to minimize numerical issues.

The set of restrictions is completed with the inclusion of mass and concentration balances to calculate the final concentration/mass of species:

.. math:: C_i^f = C_i^0 + \sum_{r \in N_{rxn,i}} \alpha_{i,r} X_r \forall i \in I_{aq} \quad (5)

.. math:: m_i^f = m_i^0 + \sum_{r \in N_{rxn,i}} \alpha_{i,r} X_r V \forall i \in I_{sp} \quad (6)

where :math:`C_i^0` is the initial concentration of species :math:`i`, :math:`N_{rxn,i}` is the set of all reactions involving species :math:`i`, :math:`X_r`
is the extent of reaction :math:`r`, :math:`I_{aq}` is the set of all aqueous species, :math:`m_i^0` is the initial amount of solid species :math:`i`,
:math:`V` is the volumetric flowrate of solvent, and :math:`I_{sp}` is the set of all precipitates.

[1] Allan M.M. Leal, Dmitrii A. Kulik, William R. Smith, and Martin O. Saar. An overview of computational methods for chemical equilibrium and kinetic calculations for
geochemical and reactive transport modeling. *Pure Appl. Chem.*, 89:597-643, 2017.
"""
# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.common.config import Bool, ConfigBlock, ConfigValue
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import (
    ControlVolume0DBlock,
    UnitModelBlockData,
    declare_process_block_class,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.tables import create_stream_table_dataframe


@declare_process_block_class("Precipitator")
class PrecipitatorData(UnitModelBlockData):
    """
    Precipitator Unit Model Class
    """

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

        prop_aqueous = self.config.property_package_aqueous
        prop_precipitate = self.config.property_package_precipitate

        # Params
        # create a set containing all reaction logkeq values
        self.merged_logkeq_dict = (
            prop_aqueous.logkeq_dict | prop_precipitate.logkeq_dict
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
            doc="Extent of the reaction.",
            units=pyunits.mol / pyunits.kg,
        )
        # log(q) for precipitation reactions
        self.log_q = pyo.Var(
            prop_precipitate.rxn_set,
            initialize=1,
            doc="log(q) var for each reaction",
            units=pyunits.dimensionless,
        )
        self.m_ref = pyo.Param(
            initialize=1,
            doc="Reference molality of 1 to make log(molalities) dimensionless",
            units=pyunits.mol / pyunits.kg,
        )

        # constraints
        @self.Constraint(
            self.flowsheet().time,
            prop_aqueous.rxn_set,
            doc="equilibrium reactions (log form) constraints (non-precipitate forming)",
        )
        def log_k_equilibrium_rxn_eqns(blk, t, r):

            return blk.log_k[r] == sum(
                prop_aqueous.stoich_dict[r][c]
                * pyo.log10(
                    blk.cv_aqueous.properties_out[t].molality_aqueous_comp[c]
                    / self.m_ref
                )
                for c in prop_aqueous.stoich_dict[r]
            )

        @self.Constraint(
            self.flowsheet().time,
            prop_precipitate.rxn_set,
            doc="equilibrium reactions (log form) constraints (precipitate forming)",
        )
        def log_q_precipitate_equilibrium_rxn_eqns(blk, t, r):

            return blk.log_q[r] == sum(
                prop_aqueous.stoich_dict[r][c]
                * pyo.log10(
                    blk.cv_aqueous.properties_out[t].molality_aqueous_comp[c]
                    / self.m_ref
                )
                for c in prop_aqueous.stoich_dict[r]
            )

        # log(q) must be <= log(k) for precipitating reactions
        @self.Constraint(
            prop_precipitate.rxn_set,
            doc="log(q) must be less than or equal to log(k) for precipitating reactions.",
        )
        def precip_rxns_log_cons(blk, r):
            return blk.log_q[r] <= self.log_k[r]

        # mole balance on all aqueous components
        @self.Constraint(
            self.flowsheet().time,
            prop_aqueous.component_list,
            doc="Aqueous Components Mole balance equations.",
        )
        def aqueous_mole_balance_eqns(blk, t, comp):
            return blk.cv_aqueous.properties_out[t].molality_aqueous_comp[
                comp
            ] == blk.cv_aqueous.properties_in[t].molality_aqueous_comp[comp] + sum(
                prop_aqueous.stoich_dict[r][comp] * blk.rxn_extent[r]
                for r in self.merged_rxns
            )

        # balance on volume coming in and out
        @self.Constraint(self.flowsheet().time, doc="volume balance equation.")
        def vol_balance(blk, t):
            return (
                blk.cv_aqueous.properties_out[t].flow_vol
                == blk.cv_aqueous.properties_in[t].flow_vol
            )

        # mole balance on all precipitate components
        @self.Constraint(
            self.flowsheet().time,
            prop_precipitate.component_list,
            doc="Precipitate Components Mole balance equations.",
        )
        def precipitate_mole_balance_eqns(blk, t, comp):
            return blk.cv_precipitate.properties_out[t].moles_precipitate_comp[
                comp
            ] == blk.cv_precipitate.properties_in[t].moles_precipitate_comp[comp] + sum(
                prop_precipitate.stoich_dict[r][comp]
                * blk.rxn_extent[r]
                * blk.cv_aqueous.properties_out[t].flow_vol
                for r in prop_precipitate.rxn_set
            )

        # objective function to minimize difference between log(k) and log(q)
        @self.Objective(
            doc="objective function minimize difference between log(k) and log(q)"
        )
        def min_logs(blk, r):
            return sum(
                (self.log_k[r] - blk.log_q[r]) ** 2 for r in prop_precipitate.rxn_set
            )

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Aqueous Inlet": self.aqueous_inlet,
                "Aqueous Outlet": self.aqueous_outlet,
                "Precipitate Inlet": self.precipitate_inlet,
                "Precipitate Outlet": self.precipitate_outlet,
            },
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):
        var_dict = {}

        for r in self.merged_rxns:
            name = f"Reaction {r} extent"
            rxn_extent = self.rxn_extent[r]

            var_dict[name] = rxn_extent

        return {"vars": var_dict}
