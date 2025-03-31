#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Simple leaching model for West Kentucky No. 13 coal refuse in H2SO4.

Authors: Andrew Lee

This is an example of how to write a custom heterogeneous reaction package for use with the
LeachTrain unit model.

"""
from pyomo.common.config import ConfigValue
from pyomo.environ import Constraint, Param, Set, Var, units

from idaes.core import ProcessBlock, ProcessBlockData, declare_process_block_class
from idaes.core.base import property_meta
from idaes.core.util.misc import add_object_reference


# -----------------------------------------------------------------------------
# Leach solution property package
@declare_process_block_class("CoalRefuseLeachingReactions")
class CoalRefuseLeachingReactionsData(
    ProcessBlockData, property_meta.HasPropertyClassMetadata
):
    """
    Reaction package for heterogeneous reactions involved in leaching REEs from
    solid West Kentucky No. 13 coal refuse using H2SO4.

    This reaction package is designed to be used with the LeachTrain or MSContactor
    unit models and assumed two streams named 'liquid' and 'solid'.

    Reaction parameters fitted to a shrinking core model using data from:

    RESEARCH PERFORMANCE FINAL REPORT, Pilot-Scale Testing of an Integrated
    Circuit for the Extraction of Rare Earth Minerals and Elements from Coal
    and Coal Byproducts Using Advanced Separation Technologies,
    Honaker, R.Q., et al., DE-FE0027035

    Includes reactions for the following components with H2SO4:

    * Rare Earth Oxides: Sc2O3, Y2O3, La2O3, Ce2O3, Pr2O3, Nd2O3, Sm2O3, Gd2O3, Dy2O3
    * Impurities: Al2O3, CaO, Fe2O3

    All reactions use the following form:

    rate[j] = eps*B*[H+]^A*(1-X[j])^(2/3)

    where X[j] is the solid phase conversion (i.e., recovery) of species j, eps is the
    mass-based pulp density and A and B are fitted parameters.

    """

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
                "Al2O3": 1.496716606,
                "Fe2O3": 0.902948175,
                "CaO": 0.159744406,
                "Sc2O3": 0.763763942,
                "Y2O3": 0.580700988,
                "La2O3": 0.443101432,
                "Ce2O3": 0.601391182,
                "Pr2O3": 0.501916124,
                "Nd2O3": 0.702951111,
                "Sm2O3": 0.578717372,
                "Gd2O3": 1.063666638,
                "Dy2O3": 0.428087853,
            },
            units=units.dimensionless,
            mutable=True,
        )

        self.B = Param(
            self.reaction_idx,
            initialize={
                "Al2O3": 405.5050676,
                "Fe2O3": 11.34710708,
                "CaO": 0.081844698,
                "Sc2O3": 0.000755073,
                "Y2O3": 0.000528612,
                "La2O3": 0.001028295,
                "Ce2O3": 0.007050046,
                "Pr2O3": 0.000522858,
                "Nd2O3": 0.006289942,
                "Sm2O3": 0.000252479,
                "Gd2O3": 0.013408467,
                "Dy2O3": 4.56708e-05,
            },
            units=units.mol * units.kg**-1 * units.hour**-1,
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
        """
        Reaction block for leaching of West Kentucky No. 13 coal refuse in H2SO4.

        """
        super().build()

        add_object_reference(self, "_params", self.config.parameters)

        self.reaction_rate = Var(
            self.params.reaction_idx,
            initialize=0,
            units=units.mol / units.litre / units.hour,
        )

        def rule_reaction_rate_eq(b, r):
            l_block = b.parent_block().liquid[b.index()]
            s_block = b.parent_block().solid[b.index()]

            h_conc = l_block.conc_mol_comp["H"]

            # Pulp density calculation
            eps = units.convert(
                s_block.flow_mass
                / (s_block.flow_mass / s_block.params.dens_mass + l_block.flow_vol),
                to_units=units.kg / units.litre,
            )

            # Empirical correlation with varying exponent,
            # strip units from acid concentration for simplicity
            return b.reaction_rate[r] == (
                eps
                * b.params.B[r]
                * (h_conc / (units.mol / units.L)) ** b.params.A[r]
                * (1 - s_block.conversion[r]) ** (2 / 3)
            )

        self.reaction_rate_eq = Constraint(
            self.params.reaction_idx, rule=rule_reaction_rate_eq
        )

    @property
    def params(self):
        return self._params
