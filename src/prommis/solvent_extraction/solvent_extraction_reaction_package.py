#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Reaction package for solvent extraction of rare earth elements using DEHPA as extractant
with TBP as a phase modifier.

Authors: Arkoprabho Dasgupta

This is an example of how to write a reaction package for rare earth elements involved in
solvent extraction.

"""
from pyomo.common.config import ConfigValue
from pyomo.environ import Constraint, Param, Set, Var, units, log10

from idaes.core import ProcessBlock, ProcessBlockData, declare_process_block_class
from idaes.core.base import property_meta
from idaes.core.util.misc import add_object_reference


# -----------------------------------------------------------------------------
# Solvent extraction property package
@declare_process_block_class("SolventExtractionReactions")
class SolventExtractionReactionsData(
    ProcessBlockData, property_meta.HasPropertyClassMetadata
):
    """
    Reaction package for the solvent extraction of rare earth elements from acidic leachate
    solution using organic extractant DEHPA and TBP as phase modifier.

    Rare earth elements : La, Pr, Ce, Dy, Nd, Sm, Gd, Y
    Impurities : Al, Ca, Fe, Sc

    This reaction package is used by the solvent extraction model. It assumes that there are
    two streams involved in the mass transfer, the aqueous liquid phase named as 'liquid' in
    the package, and the organic phase named as 'organic'.
    The material transfer amount is accounted through the term of the heterogeneous reaction
    extent term in the solvent extraction model. The amount of material transfer is quantified
    through the use of distribution coefficient.

    The distribution coefficient is defined as the ratio of the concentration in the organic
    phase to the concentration in the aqueous phase.

    D[i] = C_organic[i]/C_aqueous[i]   i = REEs

    This reaction package is for a system where DEHPA acts as the main extractant and TBP acts
    as the phase modifier in the system. The distribution correlation can be expressed in the
    following correlation.

    logD[i] = m[i]*pH + B[i]   i = REEs

    These correlations have been taken from 'PROCESSES RESEARCH PROJECT REPORT Production of
    salable rare earths products from coal and coal byproducts using advanced separation processes
    Phase 1'.

    The parameters for each of the components m[i],B[i] can be expressed as an empirical function
    of the extractant dosage (volume percentage) of the system.

    m[i] = m0[i] + m1[i]*dosage
    B[i] = B0[i] + B1[i]*log(dosage)

    If a system has x% DEHPA, then set `blk.extractant_dosage`=x, for example for 5% DEHPA system,
    set `blk.extractant_dosage`=5

    The data points for these empirical correlations have been taken from the phase 1 report.
    Certain rare earth elements do not have adequate data, hence certain functionalities could
    not be written here.

    Alongside rare earth elements, there has been impurities considered in the system. We do not
    have adequate data points for the impurities, so their distribution coefficient values has been
    calculated from REESim file and assumed constant.

    """

    def build(self):
        super().build()

        self._reaction_block_class = SolventExtractionReactionsBlock

        REE_list = ["La", "Y", "Pr", "Ce", "Nd", "Sm", "Gd", "Dy"]
        Impurity_list = ["Al", "Ca", "Fe", "Sc"]

        index_list = [f"{e}_mass_transfer" for e in (REE_list + Impurity_list)]
        element_list = REE_list + Impurity_list

        self.element_list = Set(initialize=element_list)
        self.reaction_idx = Set(initialize=index_list)

        reaction_stoichiometry = {}

        for e in REE_list:
            reaction_stoichiometry[(f"{e}_mass_transfer", "liquid", e)] = -1
            reaction_stoichiometry[(f"{e}_mass_transfer", "organic", f"{e}_o")] = 1
            reaction_stoichiometry[(f"{e}_mass_transfer", "liquid", "H")] = 3
            reaction_stoichiometry[(f"{e}_mass_transfer", "organic", "DEHPA")] = -3

        for e in Impurity_list:
            if e == "Ca":
                reaction_stoichiometry[(f"{e}_mass_transfer", "liquid", e)] = -1
                reaction_stoichiometry[(f"{e}_mass_transfer", "organic", f"{e}_o")] = 1
                reaction_stoichiometry[(f"{e}_mass_transfer", "liquid", "H")] = 2
                reaction_stoichiometry[(f"{e}_mass_transfer", "organic", "DEHPA")] = -2
            else:
                reaction_stoichiometry[(f"{e}_mass_transfer", "liquid", e)] = -1
                reaction_stoichiometry[(f"{e}_mass_transfer", "organic", f"{e}_o")] = 1
                reaction_stoichiometry[(f"{e}_mass_transfer", "liquid", "H")] = 3
                reaction_stoichiometry[(f"{e}_mass_transfer", "organic", "DEHPA")] = -3

        self.reaction_stoichiometry = reaction_stoichiometry

        self.extractant_dosage = Param(
            doc="Extractant dosage of the system", initialize=1, mutable=True
        )

        self.m0 = Param(
            self.element_list,
            initialize={
                "Ce": 0.30916,
                "Y": 1.63166,
                "Gd": 1.0225,
                "Dy": 1.70783,
                "Sm": 0.81233,
                "Nd": 0.31183,
                "La": 0.54,
                "Pr": 0.29,
                "Sc": 0,
                "Al": 0,
                "Ca": 0,
                "Fe": 0,
            },
        )

        self.m1 = Param(
            self.element_list,
            initialize={
                "Ce": 0.04816,
                "Y": 0.15166,
                "Gd": 0.0195,
                "Dy": 0.06443,
                "Sm": -0.02247,
                "Nd": 0.03763,
                "La": 0,
                "Pr": 0,
                "Sc": 0,
                "Al": 0,
                "Ca": 0,
                "Fe": 0,
            },
        )

        self.B0 = Param(
            self.element_list,
            initialize={
                "Ce": -1.66021,
                "Y": -2.12601,
                "Gd": -2.24143,
                "Dy": -2.42226,
                "Sm": -2.12172,
                "Nd": -1.62372,
                "La": -1.93,
                "Pr": -1.48,
                "Sc": 0,
                "Al": 0,
                "Ca": 0,
                "Fe": 0,
            },
        )

        self.B1 = Param(
            self.element_list,
            initialize={
                "Ce": -0.38599,
                "Y": 0.26612,
                "Gd": 0.03065,
                "Dy": -0.02538,
                "Sm": 0.17414,
                "Nd": -0.38096,
                "La": 0,
                "Pr": 0,
                "Sc": 0,
                "Al": 0,
                "Ca": 0,
                "Fe": 0,
            },
        )

        self.K1 = Param(
            self.element_list,
            initialize={
                "Ce": 0,
                "Y": 0,
                "Gd": 0,
                "Dy": 0,
                "Sm": 0,
                "Nd": 0,
                "La": 0,
                "Pr": 0,
                "Sc": 632.4976,
                "Al": 0.0531,
                "Ca": 0.0658,
                "Fe": 0.1496,
            },
        )

        self.K_corr = Param(
            self.element_list,
            initialize={
                "Ce": 0,
                "Y": 0,
                "Gd": 0,
                "Dy": 0,
                "Sm": 0,
                "Nd": 0,
                "La": 0,
                "Pr": 0,
                "Sc": 1,
                "Al": 1,
                "Ca": 1,
                "Fe": 1,
            },
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

        return self.reaction_block_class(
            *args, **kwargs, **default, initialize=initialize
        )


class _SolventExtractionReactionsBlock(ProcessBlock):
    pass


@declare_process_block_class(
    "SolventExtractionReactionsBlock", block_class=_SolventExtractionReactionsBlock
)
class SolventExtractionReactionsData(ProcessBlockData):
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
        Reaction block for solvent extraction of rare earth elements.

        """
        super().build()

        add_object_reference(self, "_params", self.config.parameters)

        self.distribution_coefficient = Var(
            self.params.element_list,
            initialize=1,
        )

        def distribution_expression(b, e):
            aq_block = b.parent_block().aqueous[b.index()]

            pH = aq_block.pH_phase["liquid"]
            return (b.distribution_coefficient[e]) == 10 ** (
                (b.params.m0[e] + b.params.extractant_dosage * b.params.m1[e]) * pH
                + (b.params.B0[e] + b.params.B1[e] * log10(b.params.extractant_dosage))
            ) * (1 - b.params.K_corr[e]) + b.params.K_corr[e] * b.params.K1[e]

        self.distribution_expression_constraint = Constraint(
            self.params.element_list, rule=distribution_expression
        )

    @property
    def params(self):
        return self._params
