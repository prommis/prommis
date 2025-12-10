#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Initial property package for precipitate.

Authors: Alejandro Garciadiego
"""

from pyomo.environ import Constraint, Param, Set, value, Var, units

from idaes.core import (
    Component,
    MaterialFlowBasis,
    Phase,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    declare_process_block_class,
)
from idaes.core.scaling import CustomScalerBase, DefaultScalingRecommendation
from idaes.core.util.initialization import fix_state_vars
from idaes.core.util.misc import add_object_reference

# -----------------------------------------------------------------------------
# Precipitate solution property package
_ree_list = [
    "Sc",
    "Y",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Sm",
    "Gd",
    "Dy",
]
_gangue_list = [
    "Al",
    "Ca",
    "Fe",
]
_solvent_list = ["H2O", "H", "Cl"]
_comp_list = _gangue_list + _ree_list + _solvent_list


class HClStrippingPropertiesScaler(CustomScalerBase):
    """
    Scaler for REE precipiate liquid property package.
    """

    CONFIG = CustomScalerBase.CONFIG

    DEFAULT_SCALING_FACTORS = {
        "flow_vol": DefaultScalingRecommendation.userInputRequired,
        "conc_mass_comp[H2O]": 1e-6,
        "conc_mass_comp[H]": 1e-1,
        "conc_mass_comp[SO4]": 1e-2,
        "conc_mass_comp[HSO4]": 1e-3,
        "conc_mass_comp[Cl]": 1e-3,
        # Use a large scaling factor for pH
        # to encourage the solver to take
        # smaller steps
        "pH_phase": 10,
    }
    for ree in _ree_list:
        DEFAULT_SCALING_FACTORS[f"conc_mass_comp[{ree}]"] = 0.1
    for contaminant in _gangue_list:
        DEFAULT_SCALING_FACTORS[f"conc_mass_comp[{contaminant}]"] = 0.1

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # Scale state variables
        self.scale_variable_by_default(model.flow_vol, overwrite=overwrite)
        for idx, var in model.conc_mass_comp.items():
            self.scale_variable_by_default(var, overwrite=overwrite)

        # Scale other variables
        if model.is_property_constructed("pH_phase"):
            self.scale_variable_by_default(
                model.pH_phase["liquid"], overwrite=overwrite
            )

        for idx, vardata in model.conc_mol_comp.items():
            self.scale_variable_by_definition_constraint(
                vardata, model.molar_concentration_constraint[idx], overwrite=overwrite
            )

        for idx, vardata in model.flow_mol_comp.items():
            self.scale_variable_by_definition_constraint(
                vardata, model.flow_mol_constraint[idx], overwrite=overwrite
            )

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        for idx, condata in model.molar_concentration_constraint.items():
            self.scale_constraint_by_component(
                condata, model.conc_mass_comp[idx], overwrite=overwrite
            )
        for idx, condata in model.flow_mol_constraint.items():
            self.scale_constraint_by_component(
                condata, model.flow_mol_comp[idx], overwrite=overwrite
            )

        if model.is_property_constructed("pH_phase"):
            self.scale_constraint_by_component(
                model.pH_constraint["liquid"],
                model.conc_mol_comp["H"],
                overwrite=overwrite,
            )

        # Why is this constraint not present in this property package?
        if model.is_property_constructed("h2o_concentration"):
            self.scale_constraint_by_component(
                model.h2o_concentration,
                model.conc_mass_comp["H2O"],
                overwrite=overwrite,
            )


@declare_process_block_class("HClStrippingParameterBlock")
class HClStrippingParameterData(PhysicalParameterBlock):
    """
    Property package for HCl solution generated in oxalate precipitator.

    Includes the following components:

    * Rare Earths: Sc, Y, La, Ce, Pr, Nd, Sm, Gd, Dy
    * Impurities: Al, Ca, Fe

    Assumes the equilibrium reaction has a fixed rate and fixed partition coefficients
    based on:

    Wang, Y., Ziemkiewicz, P., Noble, A., A Hybrid Experimental and Theoretical Approach
    to Optimize Recovery of Rare Earth Elements from Acid Mine Drainage,
    Minerals, 2022, 12. 236

    self.split can be substituted by surrogate model
    """

    def build(self):
        super().build()

        self.liquid = Phase()

        # Solvents
        self.H2O = Component()

        # Contaminants
        self.Al = Component()
        self.Ca = Component()
        self.Fe = Component()
        self.H = Component()
        self.Cl = Component()

        # REEs
        self.Sc = Component()
        self.Y = Component()
        self.La = Component()
        self.Ce = Component()
        self.Pr = Component()
        self.Nd = Component()
        self.Sm = Component()
        self.Gd = Component()
        self.Dy = Component()

        # parameter based on pH 1.28
        # TODO add surrogate model/equation
        self.split = Param(
            self.component_list,
            units=units.kg / units.kg,
            initialize={
                "H2O": 1e-20,
                "Sc": 31.61,
                "Y": 74.46,
                "La": 51.51,
                "Ce": 68.07,
                "Pr": 78,
                "Nd": 81.55,
                "Sm": 87.35,
                "Gd": 88.01,
                "Dy": 87.16,
                "Al": 0.9,
                "Ca": 20.50,
                "Fe": 2.44,
                "H": 1e-20,
                "Cl": 1e-20,
            },
        )

        self.mw = Param(
            self.component_list,
            units=units.kg / units.mol,
            initialize={
                "H2O": 18e-3,
                "Sc": 44.946e-3,
                "Y": 88.905e-3,
                "La": 138.905e-3,
                "Ce": 140.116e-3,
                "Pr": 140.907e-3,
                "Nd": 144.242e-3,
                "Sm": 150.36e-3,
                "Gd": 157.25e-3,
                "Dy": 162.50e-3,
                "Al": 26.982e-3,
                "Ca": 40.078e-3,
                "Fe": 55.845e-3,
                "H": 1.008e-3,
                "Cl": 35.453e-3,
            },
        )

        self.dissolved_elements = Set(initialize=_comp_list)

        # Assume dilute acid, density of pure water
        self.dens_mass = Param(
            initialize=1,
            units=units.kg / units.litre,
            mutable=True,
        )

        self._state_block_class = HClStrippingStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "flow_vol": {"method": None},
                "conc_mass_comp": {"method": None},
                "dens_mol": {"method": "_dens_mol"},
                "flow_mol_comp": {"method": None},
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


class _HClStrippingStateBlock(StateBlock):
    default_scaler = HClStrippingPropertiesScaler

    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        # Fix state variables
        fix_state_vars(self)

        for sbd in self.values():
            if not sbd.config.defined_state:
                sbd.h2o_concentration.deactivate()


@declare_process_block_class(
    "HClStrippingStateBlock", block_class=_HClStrippingStateBlock
)
class HClStrippingStateBlockkData(StateBlockData):
    """
    State block for HCl solution generated in oxalate precipitator.

    """

    default_scaler = HClStrippingPropertiesScaler

    def build(self):
        super().build()

        self.flow_vol = Var(
            units=units.L / units.hour,
            initialize=10,
            bounds=(1e-8, None),
        )
        self.conc_mass_comp = Var(
            self.params.dissolved_elements,
            units=units.mg / units.L,
            initialize=1e-5,
            bounds=(1e-30, None),
        )
        self.flow_mol_comp = Var(
            self.params.dissolved_elements,
            units=units.mol / units.hour,
            initialize=1e-5,
            bounds=(1e-30, None),
        )

        self.conc_mol_comp = Var(
            self.params.dissolved_elements,
            units=units.mol / units.L,
            initialize=1e-5,
            bounds=(1e-30, None),
        )

        # Concentration conversion constraint
        @self.Constraint(self.params.dissolved_elements)
        def molar_concentration_constraint(b, j):
            return (
                units.convert(
                    b.conc_mol_comp[j] * b.params.mw[j], to_units=units.mg / units.litre
                )
                == b.conc_mass_comp[j]
            )

        # Concentration conversion constraint
        @self.Constraint(self.params.dissolved_elements)
        def flow_mol_constraint(b, j):
            return (
                units.convert(
                    b.flow_vol * b.conc_mass_comp[j] / b.params.mw[j],
                    to_units=units.mol / units.hour,
                )
                == b.flow_mol_comp[j]
            )

        if not self.config.defined_state:
            # Concentration of H2O based on assumed density
            self.h2o_concentration = Constraint(
                expr=self.conc_mass_comp["H2O"] == 1e6 * units.mg / units.L
            )

    def _dens_mass(self):
        add_object_reference(self, "dens_mass", self.params.dens_mass)

    def _pH_phase(self):
        self.pH_phase = Var(
            self.params.phase_list, initialize=2, doc="pH of the solution"
        )

        @self.Constraint(self.phase_list)
        def pH_constraint(b, p):
            return 10 ** (-b.pH_phase[p]) == b.conc_mol_comp["H"] * units.L / units.mol

    def get_material_flow_terms(self, p, j):
        # Note conversion to mol/hour
        if j == "H2O":
            # Assume constant density of 1 kg/L
            return self.flow_vol * self.params.dens_mass / self.params.mw[j]
        else:
            # Need to convert from moles to mass
            return units.convert(
                self.flow_vol * self.conc_mass_comp[j] / self.params.mw[j],
                to_units=units.mol / units.hour,
            )

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        return {
            "flow_vol": self.flow_vol,
            "conc_mass_comp": self.conc_mass_comp,
        }
