#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Simple oxalate preciptation model for West Kentucky No. 13 coal refuse.

Authors: Alejandro Garciadiego

This is an example of how to write a custom heterogeneous reaction package for use with the
precipitation unit model.

"""

from pyomo.common.config import ConfigValue
from pyomo.environ import Constraint, Param, Set, Var, units

from idaes.core import ProcessBlock, ProcessBlockData, declare_process_block_class
from idaes.core.base import property_meta
from idaes.core.util.misc import add_object_reference
import idaes.core.util.scaling as iscale


# -----------------------------------------------------------------------------
# Leach solution property package
@declare_process_block_class("OxalatePrecipitationReactions")
class OxalatePrecipitationLeachingReactionsData(
    ProcessBlockData, property_meta.HasPropertyClassMetadata
):
    """
    Reaction package for heterogeneous reactions involved in oxalate precipiration of
    REEs from solid West Kentucky No. 13 coal refuse using oxalic acid.

    This reaction package is designed to be used with the precipitator
    unit model and assumed two streams named 'liquid' and 'solid'.

    Reaction parameters fitted to a shrinking core model using data from:

    visual Minteq

    """

    def build(self):
        super().build()

        self._reaction_block_class = OxalatePrecipitationReactionsBlock

        self.reaction_idx = Set(
            initialize=[
                "Al2(C2O4)3(s)",
                "Fe2(C2O4)3(s)",
                "Ca(C2O4)(s)",
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
        )

        self.reaction_stoichiometry = {
            ("Al2(C2O4)3(s)", "liquid", "Al"): 2,
            ("Al2(C2O4)3(s)", "liquid", "H2C2O4"): 3,
            ("Al2(C2O4)3(s)", "solid", "Al2(C2O4)3(s)"): -1,
            ("Al2(C2O4)3(s)", "liquid", "H"): -6,
            ("Fe2(C2O4)3(s)", "liquid", "Fe"): 2,
            ("Fe2(C2O4)3(s)", "liquid", "H2C2O4"): 3,
            ("Fe2(C2O4)3(s)", "solid", "Fe2(C2O4)3(s)"): -1,
            ("Fe2(C2O4)3(s)", "liquid", "H"): -6,
            ("Ca(C2O4)(s)", "liquid", "Ca"): 1,
            ("Ca(C2O4)(s)", "liquid", "H2C2O4"): 1,
            ("Ca(C2O4)(s)", "solid", "Ca(C2O4)(s)"): -1,
            ("Ca(C2O4)(s)", "liquid", "H"): -2,
            ("Sc2(C2O4)3(s)", "liquid", "Sc"): 2,
            ("Sc2(C2O4)3(s)", "liquid", "H2C2O4"): 3,
            ("Sc2(C2O4)3(s)", "solid", "Sc2(C2O4)3(s)"): -1,
            ("Sc2(C2O4)3(s)", "liquid", "H"): -6,
            ("Y2(C2O4)3(s)", "liquid", "Y"): 2,
            ("Y2(C2O4)3(s)", "liquid", "H2C2O4"): 3,
            ("Y2(C2O4)3(s)", "solid", "Y2(C2O4)3(s)"): -1,
            ("Y2(C2O4)3(s)", "liquid", "H"): -6,
            ("La2(C2O4)3(s)", "liquid", "La"): 2,
            ("La2(C2O4)3(s)", "liquid", "H2C2O4"): 3,
            ("La2(C2O4)3(s)", "solid", "La2(C2O4)3(s)"): -1,
            ("La2(C2O4)3(s)", "liquid", "H"): -6,
            ("Ce2(C2O4)3(s)", "liquid", "Ce"): 2,
            ("Ce2(C2O4)3(s)", "liquid", "H2C2O4"): 3,
            ("Ce2(C2O4)3(s)", "solid", "Ce2(C2O4)3(s)"): -1,
            ("Ce2(C2O4)3(s)", "liquid", "H"): -6,
            ("Pr2(C2O4)3(s)", "liquid", "Pr"): 2,
            ("Pr2(C2O4)3(s)", "liquid", "H2C2O4"): 3,
            ("Pr2(C2O4)3(s)", "solid", "Pr2(C2O4)3(s)"): -1,
            ("Pr2(C2O4)3(s)", "liquid", "H"): -6,
            ("Nd2(C2O4)3(s)", "liquid", "Nd"): 2,
            ("Nd2(C2O4)3(s)", "liquid", "H2C2O4"): 3,
            ("Nd2(C2O4)3(s)", "solid", "Nd2(C2O4)3(s)"): -1,
            ("Nd2(C2O4)3(s)", "liquid", "H"): -6,
            ("Sm2(C2O4)3(s)", "liquid", "Sm"): 2,
            ("Sm2(C2O4)3(s)", "liquid", "H2C2O4"): 3,
            ("Sm2(C2O4)3(s)", "solid", "Sm2(C2O4)3(s)"): -1,
            ("Sm2(C2O4)3(s)", "liquid", "H"): -6,
            ("Gd2(C2O4)3(s)", "liquid", "Gd"): 2,
            ("Gd2(C2O4)3(s)", "liquid", "H2C2O4"): 3,
            ("Gd2(C2O4)3(s)", "solid", "Gd2(C2O4)3(s)"): -1,
            ("Gd2(C2O4)3(s)", "liquid", "H"): -6,
            ("Dy2(C2O4)3(s)", "liquid", "Dy"): 2,
            ("Dy2(C2O4)3(s)", "liquid", "H2C2O4"): 3,
            ("Dy2(C2O4)3(s)", "solid", "Dy2(C2O4)3(s)"): -1,
            ("Dy2(C2O4)3(s)", "liquid", "H"): -6,
        }

        # parameter based on pH 1.5
        self.E_D = Param(
            self.reaction_idx,
            initialize={
                "Al2(C2O4)3(s)": 50,
                "Fe2(C2O4)3(s)": 8.659561,
                "Ca(C2O4)(s)": 14.49274,
                "Sc2(C2O4)3(s)": 6.42030,
                "Y2(C2O4)3(s)": 4.551786,
                "La2(C2O4)3(s)": 4.3717,
                "Ce2(C2O4)3(s)": 1.18848,
                "Pr2(C2O4)3(s)": 2.09604,
                "Nd2(C2O4)3(s)": 1.01030,
                "Sm2(C2O4)3(s)": 2.29617,
                "Gd2(C2O4)3(s)": 3.07276,
                "Dy2(C2O4)3(s)": 4.8608,
            },
            units=units.dimensionless,
            mutable=True,
        )

        # parameter based on pH 1.5
        self.N_D = Param(
            self.reaction_idx,
            initialize={
                "Al2(C2O4)3(s)": 0.9,
                "Fe2(C2O4)3(s)": 4.45302,
                "Ca(C2O4)(s)": 3.6495,
                "Sc2(C2O4)3(s)": 6.42030,
                "Y2(C2O4)3(s)": 4.67403,
                "La2(C2O4)3(s)": 4.6340,
                "Ce2(C2O4)3(s)": 2.737238,
                "Pr2(C2O4)3(s)": 3.44364,
                "Nd2(C2O4)3(s)": 2.419137,
                "Sm2(C2O4)3(s)": 3.7201,
                "Gd2(C2O4)3(s)": 4.1995,
                "Dy2(C2O4)3(s)": 4.73106,
            },
            units=units.dimensionless,
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


class _OxalatePrecipitationReactionsBlock(ProcessBlock):
    pass


@declare_process_block_class(
    "OxalatePrecipitationReactionsBlock",
    block_class=_OxalatePrecipitationReactionsBlock,
)
class OxalatePrecipitationReactionsData(ProcessBlockData):
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
        """
        Reaction block for leaching of West Kentucky No. 13 coal refuse in H2SO4.

        """
        super().build()

        add_object_reference(self, "_params", self.config.parameters)

        self.reaction_rate = Var(
            self.params.reaction_idx,
            initialize=0,
            units=units.mol / units.hour,
        )

        def rule_reaction_rate_eq(b, r):
            l_block = b.parent_block().liquid[b.index()]
            s_block = b.parent_block().solid[b.index()]

            return b.reaction_rate[r] == s_block.flow_mol_comp[r]

        self.reaction_rate_eq = Constraint(
            self.params.reaction_idx, rule=rule_reaction_rate_eq
        )

        iscale.set_scaling_factor(self.reaction_rate, 1e3)

    @property
    def params(self):
        return self._params
