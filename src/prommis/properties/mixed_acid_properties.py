#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Property package for mixed acid streams in University of Kentucky flowsheet.

Authors: Douglas Allan, Andrew Lee

Reference:

[1] CRC Handbook of Chemistry and Physics: First Student Edition (1988)
    Third Printing (1991). CRC Pres Inc. Boca Raton, FL
"""

from pyomo.environ import (
    Constraint,
    Param,
    PositiveReals,
    Reals,
    Set,
    units,
    Var,
)
from pyomo.common.config import Bool, ConfigValue, In, ListOf

from idaes.core import (
    Component,
    MaterialFlowBasis,
    Phase,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    declare_process_block_class,
)
from idaes.core.util.initialization import fix_state_vars
from idaes.core.util.misc import add_object_reference
from idaes.core.scaling import CustomScalerBase, DefaultScalingRecommendation

_contaminant_list = ["Na_+", "Fe_3+", "Fe_2+", "Al_3+", "Ca_2+", "Mg_2+"]
_ree_list = [
    "Sc_3+",
    "Y_3+",
    "La_3+",
    "Ce_3+",
    "Pr_3+",
    "Nd_3+",
    "Sm_3+",
    "Gd_3+",
    "Dy_3+",
]
_reagent_list = ["NaOH", "H2O2"]

_all_components_mw = {
    # Base components
    "H2O": 18.015e-3,
    "H_+": 1.008e-3,
    "Na_+": 22.990e-3,
    "Cl_-": 35.453e-3,
    # Rare Earths
    "Sc_3+": 44.946e-3,
    "Y_3+": 88.905e-3,
    "La_3+": 138.905e-3,
    "Ce_3+": 140.116e-3,
    "Pr_3+": 140.907e-3,
    "Nd_3+": 144.242e-3,
    "Sm_3+": 150.36e-3,
    "Gd_3+": 157.25e-3,
    "Dy_3+": 162.50e-3,
    # Contaminants
    "Al_3+": 26.982e-3,
    "Fe_3+": 55.845e-3,
    "Fe_2+": 55.845e-3,
    "Ca_2+": 40.078e-3,
    "Mg_2+": 24.305e-3,
    # Reagents
    "NaOH": 39.997e-3,
    "H2O2": 34.014e-3,
    # Sulfates
    "HSO4_-": 97.071e-3,
    "SO4_2-": 96.064e-3,
    # Oxalates
    "H2C2O4": 90.034e-3,
    "HC2O4_-": 89.026e-3,
    "C2O4_2-": 88.018e-3,
    # Ascorbates
    "HAsc": 176.12,
    "Asc_-": 175.11,
    "HDha": 174.11,
    "Dha_-": 173.10,
}


class MixedAcidPropertiesScaler(CustomScalerBase):
    """
    Scaler for leach solution property package.
    """

    CONFIG = CustomScalerBase.CONFIG

    DEFAULT_SCALING_FACTORS = {
        "flow_vol": DefaultScalingRecommendation.userInputRequired,
        "pressure": 1e-5,
        "temperature": 1 / 300,
        "conc_mass_comp[H2O]": 1e-6,
        "conc_mass_comp[H_+]": 1e-1,
        "conc_mass_comp[Na_+]": 1e-3,
        "conc_mass_comp[Cl_-]": 1e-3,
        # Sulfates
        "conc_mass_comp[SO4_2-]": 1e-2,
        "conc_mass_comp[HSO4_-]": 1e-3,
        # TODO revisit oxalate concentration defaults
        # in precipitation PR
        "conc_mass_comp[H2C2O4]": 1e-3,
        "conc_mass_comp[HC2O4_-]": 1e-3,
        "conc_mass_comp[C2O4_2-]": 1e-3,
        # TODO revisit ascorbate concentration values later
        "conc_mass_comp[HAsc]": 1,
        "conc_mass_comp[Asc_-]": 100,
        "conc_mass_comp[HDha]": 1,
        "conc_mass_comp[Dha_-]": 100,
        # Use a large scaling factor for pH
        # to encourage the solver to take
        # smaller steps
        "pH_phase": 10,
    }
    for ree in _ree_list:
        DEFAULT_SCALING_FACTORS[f"conc_mass_comp[{ree}]"] = 10
    for contaminant in _contaminant_list:
        DEFAULT_SCALING_FACTORS[f"conc_mass_comp[{contaminant}]"] = 1e-2
    for reagent in _reagent_list:
        DEFAULT_SCALING_FACTORS[f"conc_mass_comp[{reagent}]"] = 1

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # Scale state variables
        self.scale_variable_by_default(model.flow_vol, overwrite=overwrite)
        self.scale_variable_by_default(model.pressure, overwrite=overwrite)
        self.scale_variable_by_default(model.temperature, overwrite=overwrite)
        for idx, var in model.conc_mass_comp.items():
            self.scale_variable_by_default(var, overwrite=overwrite)

        # Scale other variables
        params = model.params
        if model.is_property_constructed("pH_phase"):
            self.scale_variable_by_default(
                model.pH_phase["liquid"], overwrite=overwrite
            )
        for idx, vardata in model.conc_mol_comp.items():
            self.scale_variable_by_definition_constraint(
                vardata, model.conc_mol_comp_eqn[idx], overwrite=overwrite
            )

        for idx, vardata in model.flow_mol_comp.items():
            self.scale_variable_by_definition_constraint(
                vardata, model.flow_mol_comp_eqn[idx], overwrite=overwrite
            )

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        for idx, condata in model.conc_mol_comp_eqn.items():
            self.scale_constraint_by_component(
                condata, model.conc_mass_comp[idx], overwrite=overwrite
            )

        for idx, condata in model.flow_mol_comp_eqn.items():
            self.scale_constraint_by_component(
                condata, model.flow_mol_comp[idx], overwrite=overwrite
            )

        if model.is_property_constructed("pH_phase"):
            self.scale_constraint_by_component(
                model.pH_phase_eqn["liquid"],
                model.conc_mol_comp["H_+"],
                overwrite=overwrite,
            )

        if model.is_property_constructed("h2o_concentration_eqn"):
            self.scale_constraint_by_component(
                model.h2o_concentration_eqn,
                model.conc_mass_comp["H2O"],
                overwrite=overwrite,
            )

        if model.is_property_constructed("inherent_equilibrium_eqn"):
            for condata in model.inherent_equilibrium_eqn.values():
                self.scale_constraint_by_nominal_value(condata)


@declare_process_block_class("MixedAcidParameterBlock")
class MixedAcidParameterData(PhysicalParameterBlock):
    """
    Property package for the mixed acid streams used in the UKy flowsheet.
    This property package is most suitable for pH values from 0 to 5. For pH
    values over 5, certain acid dissociation reactions can be assumed to be
    complete, while others might need to be tracked.

    * Default components: H2O, H_+, Cl_-
    * Additional acid systems:
      * Sulfuric acid: HSO4_-, SO4_2-
      * Oxalic acid: H2C2O4, HC2O4_-, and C2O4_2-
      * Ascorbic acid: HAsc, Asc_-, HDha, Dha_-
        (ascorbic acid, ascorbate, dehydroascorbic acid, dehydroascorbate)
    * Reagents (optional): NaOH, H2O2
    * Rare Earths: Sc_3+, Y_3+, La_3+, Ce_3+, Pr_3+, Nd_3+, Sm_3+, Gd_3+, Dy_3+
    * Other metal cations: Na_+, Al_3+, Ca_3+, Fe_3+, Fe_2+

    First dissociation of H2SO4 is assumed to be complete.
    Second dissociation governed by equilibrium (Ka2) - inherent reaction.

    """

    CONFIG = PhysicalParameterBlock.CONFIG()
    CONFIG.declare(
        "metal_cation_list",
        ConfigValue(
            default=[
                "Sc_3+",
                "Y_3+",
                "La_3+",
                "Ce_3+",
                "Pr_3+",
                "Nd_3+",
                "Sm_3+",
                "Gd_3+",
                "Dy_3+",
                "Fe_3+",
                "Al_3+",
                "Ca_2+",
            ],
            domain=ListOf(str, domain=In(_ree_list + _contaminant_list)),
            doc="Which metal cations to include as components",
        ),
    )
    CONFIG.declare(
        "reagent_list",
        ConfigValue(
            default=[],
            domain=ListOf(str, domain=In(_reagent_list)),
            doc="Which reagents to include as components",
        ),
    )
    CONFIG.declare(
        "include_sulfates",
        ConfigValue(
            default=False,
            domain=Bool,
            doc="Whether to include SO4 and HSO4 as components.",
        ),
    )
    CONFIG.declare(
        "include_oxalates",
        ConfigValue(
            default=False,
            domain=Bool,
            doc="Whether to include H2C2O4, HC2O4, and C2O4 as components.",
        ),
    )
    CONFIG.declare(
        "include_ascorbates",
        ConfigValue(
            default=False,
            domain=Bool,
            doc="Whether to include HAsc, Asc_-, HDha, and Dha_- as components.",
        ),
    )

    def build(self):
        super().build()

        self.liquid = Phase()

        # Solvent
        self.H2O = Component()

        # Acid related species
        self.add_component("H_+", Component())
        self.add_component("Cl_-", Component())
        if self.config.include_sulfates:
            self.add_component("HSO4_-", Component())
            self.add_component("SO4_2-", Component())
        if self.config.include_oxalates:
            self.add_component("H2C2O4", Component())
            self.add_component("HC2O4_-", Component())
            self.add_component("C2O4_2-", Component())
        if self.config.include_ascorbates:
            self.add_component("HAsc", Component())
            self.add_component("Asc_-", Component())
            self.add_component("HDha", Component())
            self.add_component("Dha_-", Component())

        # Metal cations
        for j in self.config.metal_cation_list:
            self.add_component(j, Component())

        # Reagents
        for j in self.config.reagent_list:
            self.add_component(j, Component())

        mw_init = {}
        for j in self.component_list:
            mw_init[j] = _all_components_mw[j]

        self.mw = Param(
            self.component_list,
            units=units.kg / units.mol,
            initialize=mw_init,
        )

        if (
            self.config.include_sulfates
            or self.config.include_oxalates
            or self.config.include_ascorbates
        ):
            self._has_inherent_reactions = True
            inherent_reaction_idx = []
            self.inherent_reaction_stoichiometry = {}
            k_eq_dict = {}
        if self.config.include_sulfates:
            # Inherent reaction for partial dissociation of HSO4
            inherent_reaction_idx.append("H2SO4_Ka2")
            self.inherent_reaction_stoichiometry[("H2SO4_Ka2", "liquid", "H_+")] = 1
            self.inherent_reaction_stoichiometry[("H2SO4_Ka2", "liquid", "HSO4_-")] = -1
            self.inherent_reaction_stoichiometry[("H2SO4_Ka2", "liquid", "SO4_2-")] = 1
            k_eq_dict["H2SO4_Ka2"] = 10**-1.99 * units.mol / units.L

        if self.config.include_oxalates:
            # First hydrogen
            inherent_reaction_idx.append("H2C2O4_Ka1")
            self.inherent_reaction_stoichiometry[("H2C2O4_Ka1", "liquid", "H_+")] = 1
            self.inherent_reaction_stoichiometry[("H2C2O4_Ka1", "liquid", "H2C2O4")] = (
                -1
            )
            self.inherent_reaction_stoichiometry[
                ("H2C2O4_Ka1", "liquid", "HC2O4_-")
            ] = 1
            k_eq_dict["H2C2O4_Ka1"] = (
                5.90e-2 * units.mol / units.L
            )  # [1], page D-102 for 25 C

            # Second hydrogen
            inherent_reaction_idx.append("H2C2O4_Ka2")
            self.inherent_reaction_stoichiometry[("H2C2O4_Ka2", "liquid", "H_+")] = 1
            self.inherent_reaction_stoichiometry[
                ("H2C2O4_Ka2", "liquid", "HC2O4_-")
            ] = -1
            self.inherent_reaction_stoichiometry[
                ("H2C2O4_Ka2", "liquid", "C2O4_2-")
            ] = 1
            # Note that there's a discrepancy between the 25 C data reported on page D-102
            # and the variable temperature data reported on page D-103. The variable
            # temperature data has Ka2 at 5.91e-5 at 0 C and monotonically decreasing to
            # 3.83e-5 at 50 C. The value it reports at 25 C is 5.18e-5
            k_eq_dict["H2C2O4_Ka2"] = (
                6.40e-5 * units.mol / units.L
            )  # [1], page D-102 for 25 C

        if self.config.include_ascorbates:
            # Ascorbic acid
            inherent_reaction_idx.append("HAsc_Ka1")
            self.inherent_reaction_stoichiometry[("HAsc_Ka1", "liquid", "H_+")] = 1
            self.inherent_reaction_stoichiometry[("HAsc_Ka1", "liquid", "HAsc")] = -1

            self.inherent_reaction_stoichiometry[("HAsc_Ka1", "liquid", "Asc_-")] = 1
            k_eq_dict["HAsc_Ka1"] = (
                10**-4.04 * units.mol / units.L
            )  # from Wolfram Alpha

            # Dehydroascorbic acid
            inherent_reaction_idx.append("HDha_Ka1")
            self.inherent_reaction_stoichiometry[("HDha_Ka1", "liquid", "H_+")] = 1
            self.inherent_reaction_stoichiometry[("HDha_Ka1", "liquid", "HDha")] = -1
            self.inherent_reaction_stoichiometry[("HDha_Ka1", "liquid", "Dha_-")] = 1
            # From Drugbank, the bioinformatics and cheminformatics database.
            # https://go.drugbank.com/drugs/DB08830
            k_eq_dict["HDha_Ka1"] = 10**-3.90 * units.mol / units.L

        if self._has_inherent_reactions:
            self.inherent_reaction_idx = Set(initialize=inherent_reaction_idx)
            self.k_eq = Param(
                self.inherent_reaction_idx,
                initialize=k_eq_dict,
                units=units.mol / units.litre,
                mutable=True,
            )

            # Fill in the inherent reaction stoichiometric
            # dictionary with zeros for the species that
            # do not participate in any reactions
            for r in self.inherent_reaction_idx:
                for p in self.phase_list:
                    for j in self.component_list:
                        if (r, p, j) not in self.inherent_reaction_stoichiometry:
                            self.inherent_reaction_stoichiometry[(r, p, j)] = 0

        # Assume dilute acid, density of pure water
        self.dens_mass = Param(
            initialize=1,
            units=units.kg / units.litre,
            mutable=True,
        )

        # Heat capacity of water
        self.cp_mol = Param(
            mutable=True,
            initialize=75.327,
            doc="Molar heat capacity of water [J/mol.K]",
            units=units.J / units.mol / units.K,
        )

        self.temperature_ref = Param(
            within=PositiveReals,
            mutable=True,
            default=298.15,
            doc="Reference temperature [K]",
            units=units.K,
        )

        self._state_block_class = MixedAcidStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "flow_vol": {"method": None},
                "conc_mass_comp": {"method": None},
                "conc_mol_comp": {"method": None},
                "dens_mass": {"method": "_dens_mass"},
            }
        )
        obj.define_custom_properties(
            {
                "pH_phase": {"method": "_pH_phase"},
            }
        )
        obj.add_default_units(
            {
                "time": units.hour,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


class _MixedAcidStateBlock(StateBlock):
    default_scaler = MixedAcidPropertiesScaler

    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        # Fix state variables
        fix_state_vars(self)

        # Deactivate inherent reaction
        # and water density constraints
        for sbd in self.values():
            if not sbd.config.defined_state:
                sbd.conc_mass_comp["H2O"].unfix()

                if self.params.config.include_sulfates:
                    sbd.conc_mass_comp["HSO4_-"].unfix()

                if self.params.config.include_oxalates:
                    sbd.conc_mass_comp["HC2O4_-"].unfix()
                    sbd.conc_mass_comp["C2O4_2-"].unfix()


@declare_process_block_class("MixedAcidStateBlock", block_class=_MixedAcidStateBlock)
class MixedAcidStateBlockData(StateBlockData):
    """
    State block for mixed acid streams in UKy flowsheet.

    """

    default_scaler = MixedAcidPropertiesScaler

    def build(self):
        super().build()

        self.flow_vol = Var(
            units=units.L / units.hour,
            initialize=10,
            bounds=(1e-8, None),
        )
        self.conc_mass_comp = Var(
            self.params.component_list,
            units=units.mg / units.L,
            initialize=1e-8,
            bounds=(1e-20, None),
        )

        self.flow_mol_comp = Var(
            self.params.component_list,
            units=units.mol / units.hour,
            initialize=1e-5,
            bounds=(1e-30, None),
        )

        self.conc_mol_comp = Var(
            self.params.component_list,
            units=units.mol / units.L,
            initialize=1e-5,
            bounds=(1e-20, None),
        )

        self.temperature = Var(
            domain=Reals,
            initialize=298.15,
            bounds=(298.1, None),
            doc="State temperature [K]",
            units=units.K,
        )

        self.pressure = Var(
            domain=Reals,
            initialize=101325.0,
            bounds=(1e3, 1e6),
            doc="State pressure [Pa]",
            units=units.Pa,
        )

        # Concentration conversion constraint
        @self.Constraint(self.params.component_list)
        def conc_mol_comp_eqn(b, j):
            return (
                units.convert(
                    b.conc_mol_comp[j] * b.params.mw[j], to_units=units.mg / units.litre
                )
                == b.conc_mass_comp[j]
            )

        @self.Constraint(self.params.component_list)
        def flow_mol_comp_eqn(b, j):
            return (
                units.convert(
                    b.flow_vol * b.conc_mass_comp[j] / b.params.mw[j],
                    to_units=units.mol / units.hour,
                )
                == b.flow_mol_comp[j]
            )

        if not self.config.defined_state:
            # Concentration of H2O based on assumed density
            self.h2o_concentration_eqn = Constraint(
                expr=self.conc_mass_comp["H2O"] == 1e6 * units.mg / units.L
            )
            if self.params._has_inherent_reactions:

                @self.Constraint(self.params.inherent_reaction_idx)
                def inherent_equilibrium_eqn(b, r):
                    # For A + B <-> C + D, this rule produces
                    # conc_mol_comp["A"] * conc_mol_comp["B"] * K
                    # == conc_mol_comp["C"] * conc_mol_comp["D"]
                    #
                    # We avoid division by bringing the reactants
                    # (the species with strictly negative stoichiometric
                    # coefficients) over to the lefthand side of the
                    # equilibrium relationship
                    lhs = self.params.k_eq[r]
                    rhs = 1
                    for (
                        r2,
                        _,
                        j,
                    ), nu_j in self.params.inherent_reaction_stoichiometry.items():
                        if r == r2 and nu_j < 0:
                            if nu_j == -1:
                                lhs *= self.conc_mol_comp[j]
                            else:
                                lhs *= self.conc_mol_comp[j] ** -nu_j
                        elif r == r2 and nu_j > 0:
                            if nu_j == 1:
                                rhs *= self.conc_mol_comp[j]
                            else:
                                rhs *= self.conc_mol_comp[j] ** nu_j

                    return lhs == rhs

    def _dens_mass(self):
        add_object_reference(self, "dens_mass", self.params.dens_mass)

    def _pH_phase(self):
        self.pH_phase = Var(
            self.params.phase_list, initialize=2, doc="pH of the solution"
        )

        @self.Constraint(self.phase_list)
        def pH_phase_eqn(b, p):
            return (
                10 ** (-b.pH_phase[p]) == b.conc_mol_comp["H_+"] * units.L / units.mol
            )

    def get_material_flow_terms(self, p, j):
        return self.flow_mol_comp[j]

    def get_material_density_terms(self, p, j):
        return units.convert(
            self.conc_mol_comp[j],  #  Has units of mol/L
            to_units=units.mol / units.m**3,
        )

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        return {
            "flow_vol": self.flow_vol,
            "conc_mass_comp": self.conc_mass_comp,
            "temperature": self.temperature,
            "pressure": self.pressure,
        }
