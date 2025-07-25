#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
property package for raw solid feed material that contains REE.

Authors: Jinliang ma
"""

from pyomo.environ import Constraint, Expression, Param, Var, value, units

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
from idaes.core.scaling import CustomScalerBase, get_scaling_factor

ree_sp_list = [
    "Ree2X",
    "Sc2X",
    "Y2X",
    "La2X",
    "Ce2X",
    "Pr2X",
    "Nd2X",
    "Sm2X",
    "Gd2X",
    "Dy2X",
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
gangue_sp_list = ["Kaolinite", "Al2O3", "SiO2", "CaCO3", "CaO", "FeS2", "Fe2O3"]


class ReeFeedPropertiesScaler(CustomScalerBase):
    """
    Scaler for REE feed solid.
    """

    CONFIG = CustomScalerBase.CONFIG

    DEFAULT_SCALING_FACTORS = {
        "flow_mass": 1e-1,
        "temperature": 1e-2,
        "enth_mol": 1e-5,
        "enth_mass": 1e-6,
        "mass_frac_comp[C]": 10,
        "mass_frac_comp[H]": 20,
        "mass_frac_comp[O]": 20,
        "mass_frac_comp[N]": 30,
        "mass_frac_comp[S]": 30,
        "mass_frac_comp[H2O]": 10,
    }
    for sp in ree_sp_list:
        DEFAULT_SCALING_FACTORS[f"mass_frac_comp[{sp}]"] = 1e5
    for sp in gangue_sp_list:
        DEFAULT_SCALING_FACTORS[f"mass_frac_comp[{sp}]"] = 3

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # Scale state variables
        self.scale_variable_by_default(model.flow_mass, overwrite=overwrite)
        self.scale_variable_by_default(model.temperature, overwrite=overwrite)
        self.scale_variable_by_default(model.enth_mol, overwrite=overwrite)
        self.scale_variable_by_default(model.enth_mass, overwrite=overwrite)
        for var in model.mass_frac_comp.values():
            self.scale_variable_by_default(var, overwrite=overwrite)

        params = model.params
        for idx, var in model.flow_mol_comp.items():
            sf = (
                get_scaling_factor(model.mass_frac_comp[idx])
                * get_scaling_factor(model.flow_mass)
                * value(params.mw_comp[idx])
            )
            self.set_variable_scaling_factor(var, sf, overwrite=overwrite)
        for idx, var in model.enth_mol_comp.items():
            self.set_variable_scaling_factor(var, 1e-5, overwrite=overwrite)

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        for idx, con in model.flow_mol_comp_constraint.items():
            sf = get_scaling_factor(model.flow_mol_comp[idx])
            self.set_constraint_scaling_factor(con, sf, overwrite=overwrite)

        if model.is_property_constructed("sum_mass_frac"):
            self.set_constraint_scaling_factor(
                model.sum_mass_frac, 1, overwrite=overwrite
            )

        for idx, con in model.enth_mol_comp_constraint.items():
            sf = get_scaling_factor(model.enth_mol_comp[idx])
            self.set_constraint_scaling_factor(con, sf, overwrite=overwrite)

        sf = get_scaling_factor(model.enth_mol)
        self.set_constraint_scaling_factor(
            model.enth_mol_constraint, sf, overwrite=overwrite
        )

        sf = get_scaling_factor(model.enth_mass)
        self.set_constraint_scaling_factor(
            model.enth_mass_constraint, sf, overwrite=overwrite
        )


# -----------------------------------------------------------------------------
@declare_process_block_class("ReeFeedParameters")
class ReeFeedParameterData(PhysicalParameterBlock):
    """
    Solid phase property package for feed REE material.

    Based on assay provided in:

    RESEARCH PERFORMANCE FINAL REPORT, Pilot-Scale Testing of an Integrated
    Circuit for the Extraction of Rare Earth Minerals and Elements from Coal
    and Coal Byproducts Using Advanced Separation Technologies,
    Honaker, R.Q., et al., DE-FE0027035

    Includes the following components:

    * Organic elements: C, H, O, N, S
    * Moisture H2O
    * Dissolvable Rare Earth Oxides: Sc2O3, Y2O3, La2O3, Ce2O3, Pr2O3, Nd2O3, Sm2O3, Gd2O3, Dy2O3
    * Insoluble Rare Earth Salts: Sc2X, Y2X, La2X, Ce2X, Pr2X, Nd2X, Sm2X, Gd2X, Dy2X, Ree2X
    * Gangue Minerals: Kaolinite, Al2O3, SiO2, CaCO3, CaO, FeS2, Fe2O3

    """

    def build(self):
        super().build()

        self.solid = Phase()

        # Organic elements
        self.C = Component()
        self.H = Component()
        self.O = Component()
        self.N = Component()
        self.S = Component()

        # moisture
        self.H2O = Component()

        # Gangue or impurity species
        self.Kaolinite = Component()
        self.Al2O3 = Component()
        self.SiO2 = Component()
        self.CaCO3 = Component()
        self.CaO = Component()
        self.FeS2 = Component()
        self.Fe2O3 = Component()

        # REE Salt (Insoluble)
        # Ree2X is the product of roasting that is insoluble
        self.Ree2X = Component()
        # Other salt is assumed insoluble but can be converted to Ree2X and oxides
        self.Sc2X = Component()
        self.Y2X = Component()
        self.La2X = Component()
        self.Ce2X = Component()
        self.Pr2X = Component()
        self.Nd2X = Component()
        self.Sm2X = Component()
        self.Gd2X = Component()
        self.Dy2X = Component()

        # REE oxides (Dissovable)
        self.Sc2O3 = Component()
        self.Y2O3 = Component()
        self.La2O3 = Component()
        self.Ce2O3 = Component()
        self.Pr2O3 = Component()
        self.Nd2O3 = Component()
        self.Sm2O3 = Component()
        self.Gd2O3 = Component()
        self.Dy2O3 = Component()

        self.mw_comp = Param(
            self.component_list,
            units=units.kg / units.mol,
            initialize={
                "C": 0.012011,
                "H": 0.0010078,
                "O": 0.015999,
                "N": 0.014007,
                "S": 0.032065,
                "H2O": 0.0010078 * 2 + 0.015999,
                "Kaolinite": (26.982 * 2 + 28.086 * 2 + 1.0078 * 4 + 15.999 * 9) * 1e-3,
                "Al2O3": (26.982 * 2 + 3 * 15.999) * 1e-3,
                "SiO2": (28.086 + 2 * 15.999) * 1e-3,
                "CaCO3": (40.078 + 12.011 + 3 * 15.999) * 1e-3,
                "CaO": (40.078 + 15.999) * 1e-3,
                "FeS2": (55.845 + 2 * 32.065) * 1e-3,
                "Fe2O3": (55.845 * 2 + 3 * 15.999) * 1e-3,
                "Ree2X": (100 * 2 + 3 * 15.999) * 1e-3,
                "Sc2X": (44.946 * 2 + 3 * 15.999) * 1e-3,
                "Y2X": (88.905 * 2 + 3 * 15.999) * 1e-3,
                "La2X": (138.905 * 2 + 3 * 15.999) * 1e-3,
                "Ce2X": (140.116 * 2 + 3 * 15.999) * 1e-3,
                "Pr2X": (140.907 * 2 + 3 * 15.999) * 1e-3,
                "Nd2X": (144.242 * 2 + 3 * 15.999) * 1e-3,
                "Sm2X": (150.36 * 2 + 3 * 15.999) * 1e-3,
                "Gd2X": (157.25 * 2 + 3 * 15.999) * 1e-3,
                "Dy2X": (162.50 * 2 + 3 * 15.999) * 1e-3,
                "Sc2O3": (44.946 * 2 + 3 * 15.999) * 1e-3,
                "Y2O3": (88.905 * 2 + 3 * 15.999) * 1e-3,
                "La2O3": (138.905 * 2 + 3 * 15.999) * 1e-3,
                "Ce2O3": (140.116 * 2 + 3 * 15.999) * 1e-3,
                "Pr2O3": (140.907 * 2 + 3 * 15.999) * 1e-3,
                "Nd2O3": (144.242 * 2 + 3 * 15.999) * 1e-3,
                "Sm2O3": (150.36 * 2 + 3 * 15.999) * 1e-3,
                "Gd2O3": (157.25 * 2 + 3 * 15.999) * 1e-3,
                "Dy2O3": (162.50 * 2 + 3 * 15.999) * 1e-3,
            },
        )

        # standard enthalpy of formation, ignore REE related compounds
        self.enth0_comp = Param(
            self.component_list,
            units=units.J / units.mol,
            initialize={
                "C": 0,
                "H": 0,
                "O": 0,
                "N": 0,
                "S": 0,
                "H2O": -285830,
                "Kaolinite": -4119600,
                "Al2O3": -1675700,
                "SiO2": -910700,
                "CaCO3": -1206920,
                "CaO": -635090,
                "FeS2": -171544,
                "Fe2O3": -825500,
                "Ree2X": 0,
                "Sc2X": 0,
                "Y2X": 0,
                "La2X": 0,
                "Ce2X": 0,
                "Pr2X": 0,
                "Nd2X": 0,
                "Sm2X": 0,
                "Gd2X": 0,
                "Dy2X": 0,
                "Sc2O3": 0,
                "Y2O3": 0,
                "La2O3": 0,
                "Ce2O3": 0,
                "Pr2O3": 0,
                "Nd2O3": 0,
                "Sm2O3": 0,
                "Gd2O3": 0,
                "Dy2O3": 0,
            },
        )

        # constant part of heat capacity, ignore REE compounds, all organic elements have the same value
        self.cp0_comp = Param(
            self.component_list,
            units=units.J / units.mol / units.K,
            initialize={
                "C": 1260,
                "H": 1260,
                "O": 1260,
                "N": 1260,
                "S": 1260,
                "H2O": 4184 * 0.018,
                "Kaolinite": 246.14,
                "Al2O3": 28.039,
                "SiO2": 31.414,
                "CaCO3": 76.009,
                "CaO": 47.56,
                "FeS2": 59.09,
                "Fe2O3": 80.623,
                "Ree2X": 0,
                "Sc2X": 0,
                "Y2X": 0,
                "La2X": 0,
                "Ce2X": 0,
                "Pr2X": 0,
                "Nd2X": 0,
                "Sm2X": 0,
                "Gd2X": 0,
                "Dy2X": 0,
                "Sc2O3": 0,
                "Y2O3": 0,
                "La2O3": 0,
                "Ce2O3": 0,
                "Pr2O3": 0,
                "Nd2O3": 0,
                "Sm2O3": 0,
                "Gd2O3": 0,
                "Dy2O3": 0,
            },
        )

        # linear function coefficient of heat capacity, ignore REE compounds, all organic elements have the same value of zero
        self.cp1_comp = Param(
            self.component_list,
            units=units.J / units.mol / units.K**2,
            initialize={
                "C": 0,
                "H": 0,
                "O": 0,
                "N": 0,
                "S": 0,
                "H2O": 0,
                "Kaolinite": 0,
                "Al2O3": 0.17156,
                "SiO2": 0.054286,
                "CaCO3": 0.046296,
                "CaO": 0.0059258,
                "FeS2": 0.024232,
                "Fe2O3": 0.09936,
                "Ree2X": 0,
                "Sc2X": 0,
                "Y2X": 0,
                "La2X": 0,
                "Ce2X": 0,
                "Pr2X": 0,
                "Nd2X": 0,
                "Sm2X": 0,
                "Gd2X": 0,
                "Dy2X": 0,
                "Sc2O3": 0,
                "Y2O3": 0,
                "La2O3": 0,
                "Ce2O3": 0,
                "Pr2O3": 0,
                "Nd2O3": 0,
                "Sm2O3": 0,
                "Gd2O3": 0,
                "Dy2O3": 0,
            },
        )

        self.dens_mass = Param(
            units=units.kg / units.m**3,
            initialize=2000,
            mutable=True,
        )

        self.temperature_ref = Param(
            initialize=298.15,
            units=units.K,
        )

        self._state_block_class = ReeFeedStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "flow_mol_comp": {"method": None, "units": "mol/s"},
                "temperature": {"method": None, "units": "K"},
                "cp_mol_comp": {"method": None, "units": "J/mol/K"},
                "cp_mol": {"method": None, "units": "J/mol/K"},
                "cp_mass": {"method": None, "units": "J/kg/K"},
                "enth_mol_comp": {"method": None, "units": "J/mol"},
                "enth_mol": {"method": None, "units": "J/mol"},
                "enth_mass": {"method": None, "units": "J/kg"},
                "enth_mol_phase": {"method": None, "units": "J/mol"},
            }
        )
        obj.add_default_units(
            {
                "time": units.s,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


class _ReeFeedStateBlock(StateBlock):
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
                sbd.sum_mass_frac.deactivate()


@declare_process_block_class("ReeFeedStateBlock", block_class=_ReeFeedStateBlock)
class ReeFeedStateBlockData(StateBlockData):
    """
    State block for solid West Kentucky No. 13 coal waste.

    """

    def build(self):
        # super().build()
        super(ReeFeedStateBlockData, self).build()

        self.flow_mass = Var(
            units=units.kg / units.s,
            initialize=1,
            bounds=(1e-8, None),
            doc="Mass flow rate [kg/s]",
        )

        self.temperature = Var(
            initialize=298.15,
            units=units.K,
            bounds=(270, 3000),
            doc="Temperature [K]",
        )

        self.mass_frac_comp = Var(
            self.params.component_list,
            initialize=0.1,
            units=units.kg / units.kg,
            bounds=(1e-8, None),
            doc="Mass fraction of individual component",
        )

        self.flow_mol_comp = Var(
            self.params.component_list,
            units=units.mol / units.s,
            initialize=0.1,
            bounds=(0, None),
            doc="Component molar flowrate [mol/s]",
        )

        def rule_mole_frac_comp(b, j):
            return (
                b.mass_frac_comp[j]
                / b.params.mw_comp[j]
                / sum(
                    b.mass_frac_comp[i] / b.params.mw_comp[i]
                    for i in self.params.component_list
                )
            )

        self.mole_frac_comp = Expression(
            self.params.component_list, rule=rule_mole_frac_comp
        )

        def rule_flow_mol_comp(b, j):
            return (
                b.flow_mol_comp[j]
                == b.flow_mass * b.mass_frac_comp[j] / b.params.mw_comp[j]
            )

        self.flow_mol_comp_constraint = Constraint(
            self.params.component_list, rule=rule_flow_mol_comp
        )

        def rule_flow_mol(b):
            return sum(b.flow_mol_comp[i] for i in self.params.component_list)

        self.flow_mol = Expression(rule=rule_flow_mol)

        if not self.config.defined_state:
            self.sum_mass_frac = Constraint(
                expr=1
                == sum(self.mass_frac_comp[j] for j in self.params.component_list)
            )

        self.enth_mol_comp = Var(
            self.params.component_list,
            initialize=0,
            units=units.J / units.mol,
            doc="Component molar enthalpy [J/mol]",
        )

        self.enth_mol = Var(
            initialize=0,
            units=units.J / units.mol,
            doc="Molar enthalpy of mixture [J/mol]",
        )

        self.enth_mass = Var(
            initialize=0,
            units=units.J / units.kg,
            doc="Mass enthalpy [J/kg]",
        )

        def rule_cp_mol_comp(b, i):
            return b.params.cp0_comp[i] + b.params.cp1_comp[i] * b.temperature

        self.cp_mol_comp = Expression(self.params.component_list, rule=rule_cp_mol_comp)

        def rule_cp_mol(b):
            return sum(
                b.cp_mol_comp[i] * b.mole_frac_comp[i] for i in b.params.component_list
            )

        self.cp_mol = Expression(rule=rule_cp_mol)

        def rule_cp_mass(b):
            return sum(
                b.cp_mol_comp[i] * b.mass_frac_comp[i] / b.params.mw_comp[i]
                for i in b.params.component_list
            )

        self.cp_mass = Expression(rule=rule_cp_mass)

        def rule_enth_mol_comp(b, i):
            t = b.temperature
            tref = b.params.temperature_ref
            return b.enth_mol_comp[i] * 1e-4 == 1e-4 * (
                b.params.enth0_comp[i]
                + b.params.cp0_comp[i] * (t - tref)
                + 0.5 * b.params.cp1_comp[i] * (t * t - tref * tref)
            )

        self.enth_mol_comp_constraint = Constraint(
            self.params.component_list, rule=rule_enth_mol_comp
        )

        def rule_enth_mol(b):
            return b.enth_mol == sum(
                b.enth_mol_comp[i] * b.mole_frac_comp[i]
                for i in b.params.component_list
            )

        self.enth_mol_constraint = Constraint(rule=rule_enth_mol)

        def rule_enth_mass(b):
            return b.enth_mass == sum(
                b.enth_mol_comp[i] * b.mass_frac_comp[i] / b.params.mw_comp[i]
                for i in b.params.component_list
            )

        self.enth_mass_constraint = Constraint(rule=rule_enth_mass)

        def rule_enth_phase(b, p):
            # This property module only has one phase
            return b.enth_mol

        self.enth_mol_phase = Expression(self.params.phase_list, rule=rule_enth_phase)

    def get_material_density_terms(self, p, j):
        return self.params.dens_mass

    def get_energy_density_terms(self, p):
        if not self.is_property_constructed("energy_density_terms"):
            try:

                def rule_energy_density_terms(b, p):
                    return b.enth_mass * b.params.dens_mass

                self.energy_density_terms = Expression(
                    self.params.phase_list, rule=rule_energy_density_terms
                )
            except AttributeError:
                self.del_component(self.energy_density_terms)
        return self.energy_density_terms[p]

    def get_material_flow_terms(self, p, j):
        return self.flow_mol_comp[j]

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def get_enthalpy_flow_terms(self, p):
        if not self.is_property_constructed("enthalpy_flow_terms"):
            try:

                def rule_enthalpy_flow_terms(b, p):
                    return sum(
                        b.flow_mol_comp[i] * b.enth_mol_comp[i]
                        for i in b.params.component_list
                    )

                self.enthalpy_flow_terms = Expression(
                    self.params.phase_list, rule=rule_enthalpy_flow_terms
                )
            except AttributeError:
                self.del_component(self.enthalpy_flow_terms)
        return self.enthalpy_flow_terms[p]

    def define_state_vars(self):
        return {
            "flow_mass": self.flow_mass,
            "temperature": self.temperature,
            "mass_frac_comp": self.mass_frac_comp,
        }
