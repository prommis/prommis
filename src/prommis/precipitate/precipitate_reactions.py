#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Reaction package for heterogeneous reactions involved in oxalate precipitation of
REEs from solid West Kentucky No. 13 coal refuse using oxalic acid.

Authors: Alejandro Garciadiego, Bo-Xun Wang

This is an example of how to write a custom heterogeneous reaction package for use with the
precipitation unit model.

"""

from pyomo.common.config import ConfigValue
from pyomo.environ import Constraint, Param, Set, Var, units

from idaes.core.util.config import (
    is_physical_parameter_block,
    is_reaction_parameter_block,
)

from idaes.core import ProcessBlock, ProcessBlockData, declare_process_block_class
from idaes.core.base import property_meta
from idaes.core.util.misc import add_object_reference
import idaes.core.util.scaling as iscale


# -----------------------------------------------------------------------------
# Precipitation property package
@declare_process_block_class("OxalatePrecipitationReactions")
class OxalatePrecipitationLeachingReactionsData(
    ProcessBlockData, property_meta.HasPropertyClassMetadata
):
    """
    Reaction package for heterogeneous reactions involved in oxalate precipitation of
    REEs from solid West Kentucky No. 13 coal refuse using oxalic acid.

    This reaction package is designed to be used with the precipitator
    unit model and assumes two streams named 'liquid' and 'solid'.

    Reaction parameters fitted to a shrinking core model using data from:

    visual Minteq

    Assumes the equilibrium reaction can be described by the equilibrium equation.

    .. math:: Conversion_{c} = \exp(-\frac{\epsilon}{Oxalic Acid Dosage}^{n_{DA}})

    where :math:`Conversion` is the conversion of component c, :math:`\epsilon` and :math:`n_{DA}` are the parameters
    estimated based on Minteq data, :math:`Oxalic Acid Dosage` is the amount of oxalic acid added into the precipitator.

    The details for obtaining the equilibrium parameters can be found at the Jupyter notebook, "Parameter estimation tutorial.ipynb".

    """

    def build(self):
        super().build()

        self._reaction_block_class = OxalatePrecipitationReactionsBlock

        trivalent_list = ["Al", "Fe", "Sc", "Y", "La", "Ce", "Pr", "Nd", "Sm", "Gd", "Dy"]
        divalent_list = ["Ca"]
        element_list = trivalent_list + divalent_list

        index_list = [f"{e}2(C2O4)3(s)" for e in trivalent_list] + [
            f"{e}(C2O4)(s)" for e in divalent_list
        ]

        self.element_list = Set(initialize=element_list)
        self.reaction_idx = Set(initialize=index_list)

        reaction_stoichiometry = {}
        for e in trivalent_list:
            rxn = f"{e}2(C2O4)3(s)"
            reaction_stoichiometry[(rxn, "liquid", e)] = 2
            reaction_stoichiometry[(rxn, "liquid", "H2C2O4")] = 3
            reaction_stoichiometry[(rxn, "solid", rxn)] = -1
            reaction_stoichiometry[(rxn, "liquid", "H")] = -6

        for e in divalent_list:
            rxn = f"{e}(C2O4)(s)"
            reaction_stoichiometry[(rxn, "liquid", e)] = 1
            reaction_stoichiometry[(rxn, "liquid", "H2C2O4")] = 1
            reaction_stoichiometry[(rxn, "solid", rxn)] = -1
            reaction_stoichiometry[(rxn, "liquid", "H")] = -2

        self.reaction_stoichiometry = reaction_stoichiometry

        # Equilibrium parameter (\epsilon) based on pH 1.5
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

        # Equilibrium parameter (n_{DA}) based on pH 1.5
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
            description="""A reference to an instance of the Heterogeneous Reaction Parameter
    Block associated with this property package.""",
        ),
    )

    def build(self):
        """
        Reaction block for precipitation of West Kentucky No. 13 coal refuse using oxalic acid.

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