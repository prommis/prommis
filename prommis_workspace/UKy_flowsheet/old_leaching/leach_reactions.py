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
Simple old_leaching model for West Kentucky No. 13 coal refuse in H2SO4.

Authors: Andrew Lee
"""
from pyomo.environ import Constraint, Expression, Param, Set, units, Var
from pyomo.common.config import ConfigValue

from idaes.core import (
    declare_process_block_class,
    ProcessBlockData,
    ProcessBlock,
)
from idaes.core.base import property_meta
from idaes.core.util.misc import add_object_reference


# -----------------------------------------------------------------------------
# Leach solution property package
@declare_process_block_class("CoalRefuseLeachingReactions")
class CoalRefuseLeachingReactionsData(ProcessBlockData, property_meta.HasPropertyClassMetadata):
    def build(self):
        super().build()

        self._reaction_block_class = CoalRefuseLeachingReactionsBlock

        self.reaction_idx = Set(
            initialize=[
                "Al2O3",
                "Fe2O3",
                "CaO",
                "Sc2O3",
                "Y2O3",
                "La2O3",
                "Ce2O3",
                "Pr2O3",
                "Nd2O3",
                "Sm2O3",
                "Gd2O3",
                "Dy2O3",
            ]
        )

        self.reaction_stoichiometry = {
            ("Al2O3", "liquid", "Al"): 2,
            ("Al2O3", "liquid", "H2O"): 3,
            ("Al2O3", "solid", "Al2O3"): -1,
            ("Al2O3", "liquid", "H"): -6,
            ("Fe2O3", "liquid", "Fe"): 2,
            ("Fe2O3", "liquid", "H2O"): 3,
            ("Fe2O3", "solid", "Fe2O3"): -1,
            ("Fe2O3", "liquid", "H"): -6,
            ("CaO", "liquid", "Ca"): 1,
            ("CaO", "liquid", "H2O"): 1,
            ("CaO", "solid", "CaO"): -1,
            ("CaO", "liquid", "H"): -2,
            ("Sc2O3", "liquid", "Sc"): 2,
            ("Sc2O3", "liquid", "H2O"): 3,
            ("Sc2O3", "solid", "Sc2O3"): -1,
            ("Sc2O3", "liquid", "H"): -6,
            ("Y2O3", "liquid", "Y"): 2,
            ("Y2O3", "liquid", "H2O"): 3,
            ("Y2O3", "solid", "Y2O3"): -1,
            ("Y2O3", "liquid", "H"): -6,
            ("La2O3", "liquid", "La"): 2,
            ("La2O3", "liquid", "H2O"): 3,
            ("La2O3", "solid", "La2O3"): -1,
            ("La2O3", "liquid", "H"): -6,
            ("Ce2O3", "liquid", "Ce"): 2,
            ("Ce2O3", "liquid", "H2O"): 3,
            ("Ce2O3", "solid", "Ce2O3"): -1,
            ("Ce2O3", "liquid", "H"): -6,
            ("Pr2O3", "liquid", "Pr"): 2,
            ("Pr2O3", "liquid", "H2O"): 3,
            ("Pr2O3", "solid", "Pr2O3"): -1,
            ("Pr2O3", "liquid", "H"): -6,
            ("Nd2O3", "liquid", "Nd"): 2,
            ("Nd2O3", "liquid", "H2O"): 3,
            ("Nd2O3", "solid", "Nd2O3"): -1,
            ("Nd2O3", "liquid", "H"): -6,
            ("Sm2O3", "liquid", "Sm"): 2,
            ("Sm2O3", "liquid", "H2O"): 3,
            ("Sm2O3", "solid", "Sm2O3"): -1,
            ("Sm2O3", "liquid", "H"): -6,
            ("Gd2O3", "liquid", "Gd"): 2,
            ("Gd2O3", "liquid", "H2O"): 3,
            ("Gd2O3", "solid", "Gd2O3"): -1,
            ("Gd2O3", "liquid", "H"): -6,
            ("Dy2O3", "liquid", "Dy"): 2,
            ("Dy2O3", "liquid", "H2O"): 3,
            ("Dy2O3", "solid", "Dy2O3"): -1,
            ("Dy2O3", "liquid", "H"): -6,
        }

        self.A = Param(
            self.reaction_idx,
            initialize={
                "Al2O3": 503.5875248,
                "Fe2O3": 322.6139099,
                "CaO": 109.3157001,
                "Sc2O3": 0.04725527,
                "Y2O3": 0.079996528,
                "La2O3": 0.202377027,
                "Ce2O3": 0.617265745,
                "Pr2O3": 0.054514632,
                "Nd2O3": 0.288997052,
                "Sm2O3": 0.031614808,
                "Gd2O3": 0.072328712,
                "Dy2O3": 0.00949131,
            },
            units=units.L/units.mol/units.hour,
            mutable=True,
        )

        self.B = Param(
            self.reaction_idx,
            initialize={
                "Al2O3": 0.0021735,
                "Fe2O3": 0.026437074,
                "CaO": 0.110283447,
                "Sc2O3": 0.009447767,
                "Y2O3": 0.028245478,
                "La2O3": 0.082176004,
                "Ce2O3": 0.084370281,
                "Pr2O3": 0.112379996,
                "Nd2O3": 0.09284766,
                "Sm2O3": 0.045429147,
                "Gd2O3": 0.10010438,
                "Dy2O3": 0.052401701,
            },
            units=units.hour**-1,
            mutable=True,
        )

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": units.hour,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )

    @property
    def reaction_block_class(self):
        if self._reaction_block_class is not None:
            return self._reaction_block_class
        else:
            raise AttributeError(
                "{} has not assigned a ReactionBlock class to be associated "
                "with this reaction package. Please contact the developer of "
                "the reaction package.".format(self.name)
            )

    def build_reaction_block(self, *args, **kwargs):
        """
        Methods to construct a ReactionBlock associated with this
        ReactionParameterBlock. This will automatically set the parameters
        construction argument for the ReactionBlock.

        Returns:
            ReactionBlock

        """
        default = kwargs.pop("default", {})
        initialize = kwargs.pop("initialize", {})

        if initialize == {}:
            default["parameters"] = self
        else:
            for i in initialize.keys():
                initialize[i]["parameters"] = self

        return self.reaction_block_class(  # pylint: disable=not-callable
            *args, **kwargs, **default, initialize=initialize
        )


class _CoalRefuseLeachingReactionsBlock(ProcessBlock):
    pass


@declare_process_block_class(
    "CoalRefuseLeachingReactionsBlock", block_class=_CoalRefuseLeachingReactionsBlock
)
class CoalRefuseLeachingReactionsData(ProcessBlockData):
    # Create Class ConfigBlock
    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare(
        "parameters",
        ConfigValue(
            # TODO
            # domain=is_reaction_parameter_block,
            description="""A reference to an instance of the Heterogeneous Reaction Parameter
    Block associated with this property package.""",
        ),
    )

    def build(self):
        super().build()

        add_object_reference(self, "_params", self.config.parameters)

        self.reaction_rate = Var(
            self.params.reaction_idx,
            initialize=0,
            units=units.mol/units.litre/units.hour
        )

        # TODO
        l_block = self.parent_block().liquid[self.index()]
        s_block = self.parent_block().solid[self.index()]

        h_conc = l_block.conc_mole_acid["H"]

        def rule_c_max(b, j):
            if j in ["H2O", "H", "HSO4", "SO4"]:
                return Expression.Skip

            if j == "Ca":
                m = "CaO"
                n = 1
            else:
                m = j+"2O3"
                n = 2
            return (
                units.convert(
                    l_block.conc_mass_metals[j] / l_block.params.mw[j],
                    to_units=units.mol/units.litre,
                ) +
                n*s_block.flow_mass*s_block.mass_frac_comp[m]/l_block.flow_vol/s_block.params.mw[m]
            )
        self.c_max = Expression(l_block.component_list, rule=rule_c_max)

        def rule_reaction_rate_eq(b, r):
            if r.startswith("Y"):
                j = "Y"
            else:
                j = r[:2]

            return b.reaction_rate[r] == units.convert(
                b.params.A[r]*h_conc**2 + b.params.B[r]*self.c_max[j],
                to_units=units.mol/units.litre/units.hour,
            )

        self.reaction_rate_eq = Constraint(self.params.reaction_idx, rule=rule_reaction_rate_eq)

    @property
    def params(self):
        return self._params