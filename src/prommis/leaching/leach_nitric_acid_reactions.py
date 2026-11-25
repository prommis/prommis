#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Simple leaching model for coal-derived solids in nitric acid.

Authors: Bo-Xun Wang

This is an example of how to write a custom heterogeneous reaction package for use with the
LeachTrain unit model.

"""

from pyomo.common.config import ConfigDict, ConfigValue
from pyomo.environ import Expression, Param, Set, units

from idaes.core import ProcessBlock, ProcessBlockData, declare_process_block_class
from idaes.core.base import property_meta
from idaes.core.scaling import CustomScalerBase
from idaes.core.util.misc import add_object_reference

_metal_list = ["La", "Y", "Pr", "Ce", "Nd", "Sm", "Gd", "Dy", "Al", "Ca", "Fe", "Sc"]

_default_aqueous_aliases = {e: e for e in _metal_list}
_default_aqueous_aliases["H"] = "H"
_default_aqueous_aliases["H2O"] = "H2O"


class NitricAcidLeachReactionScaler(CustomScalerBase):
    """
    Scaler for the nitric acid leach reaction package.
    No variables or constraints, so no need for scaling.
    """

    DEFAULT_SCALING_FACTORS = {}

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        pass

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        pass


# -----------------------------------------------------------------------------
# Leach solution property package
@declare_process_block_class("NitricAcidLeachingReactionParameterBlock")
class NitricAcidLeachingReactionParameterData(
    ProcessBlockData, property_meta.HasPropertyClassMetadata
):
    """
    Reaction package for heterogeneous reactions involved in leaching REEs from
    solid coal-derived material using nitric acid.

    This reaction package is designed to be used with the LeachTrain or MSContactor
    unit models and assumed two streams named 'liquid' and 'solid'.

    Reaction parameters fitted to a shrinking core model using data from:

    Journal Paper: Sustainable extraction of rare earth elements from coal fly ash
    leachates using a recyclable ionic liquid. Liu, T., Huang, C.-H., npj Materials Sustainability,
    4:2, 2026

    Includes reactions for the following components with HNO3:

    * Rare Earth Oxides: Sc2O3, Y2O3, La2O3, Ce2O3, Pr2O3, Nd2O3, Sm2O3, Gd2O3, Dy2O3
    * Impurities: Al2O3, CaO, Fe2O3

    All reactions use the following form:

    rate[j] = eps*B*[H+]^A*(1-X[j])^(2/3)

    where X[j] is the solid phase conversion (i.e., recovery) of species j, eps is the
    mass-based pulp density and A and B are fitted parameters.

    """

    CONFIG = ConfigDict()
    CONFIG.declare(
        "aqueous_aliases",
        ConfigValue(
            default=_default_aqueous_aliases,
            domain=dict,
            description="Dictionary matching element names in the element "
            "list to elements in the aqueous stream. Used when, e.g., cation "
            "charges are included in the element name, as is done in the "
            "MixedAcidProperties.",
        ),
    )

    def build(self):
        super().build()

        self._reaction_block_class = NitricAcidLeachingReactionBlock

        reaction_idx = []
        reaction_stoich = {}

        H = self.config.aqueous_aliases["H"]
        H2O = self.config.aqueous_aliases["H2O"]

        for j in _metal_list:
            j_alias = self.config.aqueous_aliases[j]
            if j == "Ca":
                # CaO + 2H_+ -> H2O + Ca_2+
                oxide = "CaO"
                reaction_stoich[("CaO", "liquid", j_alias)] = 1
                reaction_stoich[("CaO", "liquid", H2O)] = 1
                reaction_stoich[("CaO", "solid", "CaO")] = -1
                reaction_stoich[("CaO", "liquid", H)] = -2
            else:
                # Reaction of the form {j}2O3 + 6H_+ -> 2{j}_3+ + 3H2O
                oxide = f"{j}2O3"
                reaction_stoich[(oxide, "liquid", j_alias)] = 2
                reaction_stoich[(oxide, "liquid", H2O)] = 3
                reaction_stoich[(oxide, "solid", oxide)] = -1
                reaction_stoich[(oxide, "liquid", H)] = -6

            reaction_idx.append(oxide)

        self.reaction_idx = Set(initialize=reaction_idx)

        self.reaction_stoichiometry = reaction_stoich

        # assumed identical for all oxides
        self.A = Param(
            self.reaction_idx,
            initialize={
                "Al2O3": 0.37,
                "Fe2O3": 0.37,
                "CaO": 0.37,
                "Sc2O3": 0.37,
                "Y2O3": 0.37,
                "La2O3": 0.37,
                "Ce2O3": 0.37,
                "Pr2O3": 0.37,
                "Nd2O3": 0.37,
                "Sm2O3": 0.37,
                "Gd2O3": 0.37,
                "Dy2O3": 0.37,
            },
            units=units.dimensionless,
            mutable=True,
        )

        self.B = Param(
            self.reaction_idx,
            initialize={
                "Al2O3": 0.183,
                "Fe2O3": 0.0440,
                "CaO": 1.160,
                "Sc2O3": 2.64e-5,
                "Y2O3": 2.23e-5,
                "La2O3": 2.86e-5,
                "Ce2O3": 3.95e-5,
                "Pr2O3": 1.19e-6,
                "Nd2O3": 2.52e-5,
                "Sm2O3": 4.68e-6,
                "Gd2O3": 3.07e-6,
                "Dy2O3": 2.29e-6,
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


class _NitricAcidLeachingReactionBlock(ProcessBlock):
    pass


@declare_process_block_class(
    "NitricAcidLeachingReactionBlock", block_class=_NitricAcidLeachingReactionBlock
)
class NitricAcidLeachingReactionData(ProcessBlockData):
    default_scaler = NitricAcidLeachReactionScaler
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
        Reaction block for leaching of coal-derived solids in HNO3.

        """
        super().build()

        add_object_reference(self, "_params", self.config.parameters)

        def rule_reaction_rate(b, r):
            l_block = b.parent_block().liquid[b.index()]
            s_block = b.parent_block().solid[b.index()]

            H_idx = b.params.config["aqueous_aliases"]["H"]
            h_conc = l_block.conc_mol_comp[H_idx]

            # Pulp density calculation
            eps = units.convert(
                s_block.flow_mass
                / (s_block.flow_mass / s_block.params.dens_mass + l_block.flow_vol),
                to_units=units.kg / units.litre,
            )

            # strip units from acid concentration for simplicity
            return (
                eps
                * b.params.B[r]
                * (h_conc / (units.mol / units.L)) ** b.params.A[r]
                * (1 - s_block.conversion_comp[r]) ** (2 / 3)
            )

        self.reaction_rate = Expression(
            self.params.reaction_idx, rule=rule_reaction_rate
        )

    @property
    def params(self):
        return self._params
