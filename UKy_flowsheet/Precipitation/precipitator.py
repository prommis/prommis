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
Basic heater/cooler models
"""

__author__ = "John Eslick"

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

# Import IDAES cores
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    EnergyBalanceType,
    MomentumBalanceType,
    MaterialBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.math import smooth_min, smooth_abs
import idaes.core.util.scaling as iscale
from pyomo.util.check_units import assert_units_consistent

import workspace.UKy_flowsheet.Precipitation.precip_prop as precip_prop
from idaes.core.util.config import is_physical_parameter_block
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


def _make_precipitator_config_block(config):
    """
    Declare configuration options for HeaterData block.
    """
    config.declare(
        "property_package_aqueous",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for aqueous control volume",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    config.declare(
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
    config.declare(
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
    config.declare(
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


@declare_process_block_class("Precipitator")
class PrecipitatorData(UnitModelBlockData):
    """ """

    CONFIG = UnitModelBlockData.CONFIG()
    _make_precipitator_config_block(CONFIG)

    def build(self):
        """Building model

        Args:
            None
        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super().build()
        # Add Control Volume

        config = self.config

        self._add_scaling_parameters()
        cv_aqueous, cv_precipitate = self._add_control_volumes()

        prop_aq = config.property_package_aqueous
        prop_s = config.property_package_precipitate

        self._add_variables_and_refs(prop_aq, prop_s)

        @self.Expression(
            self.flowsheet().time,
            prop_s.solid_set,
            doc="Calculated molality of precipitates in the inlet.",
        )
        def molality_precipitate_comp_in(blk, t, comp):
            return (
                blk.cv_precipitate.properties_in[t].flow_mol_comp[comp]
                / blk.flow_mass_h2o_key[t]
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for the aqueous phase outlet mass flow",
        )
        def flow_eqn(blk, t):
            return blk.cv_aqueous.properties_out[
                t
            ].flow_mass == blk.cv_aqueous.properties_in[t].flow_mass - sum(
                blk.molality_precipitate_comp[t, i] for i in prop_s.solid_set
            )

        # @self.Constraint(self.flowsheet().time, doc="")
        # def temperature_aq_eqn(blk, t):
        #     return (
        #         blk.cv_aqueous.properties_out[t].temperature
        #         == blk.cv_aqueous.properties_in[t].temperature
        #     )
        #
        # @self.Constraint(self.flowsheet().time)
        # def temperature_s_eqn(blk, t):
        #     return (
        #         blk.cv_precipitate.properties_out[t].temperature
        #         == blk.cv_aqueous.properties_in[t].temperature
        #     )

        # @self.Constraint(self.flowsheet().time)
        # def pressure_eqn(blk, t):
        #     return (
        #         blk.cv_aqueous.properties_out[t].pressure
        #         == blk.cv_aqueous.properties_in[t].pressure
        #     )

        @self.Constraint(
            self.flowsheet().time, prop_aq.key_set, doc="Key component concentrations"
        )
        def molality_key_comp_eqn(blk, t, comp):
            if comp == "H^+":
                # Since OH^- is not a key component but is allowed to be used in
                # reactions, we need to translate it into H^+ which is a key component
                # by basically replacing -OH with -H2O + H. Since H^+ is generated in
                # the H2O -> H + OH reaction the OH concentration counts against the
                # key concentration of H^+, so the key concentration of H^+ can even
                # be negative. A negative H^+ key concentration just means the OH^-
                # concentration is more than the H^+ concentration.
                from_oh = sum(
                    blk.cv_aqueous.properties_in[t].molality_comp[i]
                    * d["stoich"]["OH^-"]
                    for i, d in prop_aq.aq_dict.items()
                    if "OH^-" in d["stoich"]
                )
            else:
                from_oh = 0
            return blk.molality_key_comp[t, comp] == (
                blk.cv_aqueous.properties_in[t].molality_comp[comp]
                + from_oh
                - sum(  # stoichiometry coefficient is negative
                    blk.cv_aqueous.properties_in[t].molality_comp[i] * d["stoich"][comp]
                    for i, d in prop_aq.aq_dict.items()
                    if comp in d["stoich"]
                )
                - sum(  # stoichiometry coefficient is negative
                    blk.molality_precipitate_comp_in[t, i] * d["stoich"][comp]
                    for i, d in prop_s.solid_dict.items()
                    if comp in d["stoich"]
                )
            )

        @self.Constraint(
            self.flowsheet().time, prop_aq.key_set, doc="Mass balance equations."
        )
        def mass_balance(blk, t, comp):
            return blk.molality_key_comp[t, comp] == (
                blk.molality_aq_comp[t, comp]
                - sum(  # stoichiometry coefficient is negative
                    blk.molality_aq_comp[t, i] * d["stoich"][comp]
                    for i, d in prop_aq.aq_dict.items()
                    if comp in d["stoich"]
                )
                - sum(
                    blk.molality_precipitate_comp[t, i] * d["stoich"][comp]
                    for i, d in prop_s.solid_dict.items()
                    if comp in d["stoich"]
                )
            )

        @self.Constraint(
            self.flowsheet().time,
            prop_aq.aq_set,
            doc="Aqueous equilibrium equations.",
        )
        def equilibrium(blk, t, i):
            r = prop_aq.aq_dict[i]
            return r["logK"] == (
                blk.log10_molality_aq[t, i]
                + sum(
                    r["stoich"][j] * blk.log10_molality_aq[t, j]
                    for j in r["stoich"]
                    if j != "H2O"
                )
                + blk.log10_act_coeff[t, i]
                + sum(
                    r["stoich"][j] * blk.log10_act_coeff[t, j]
                    for j in r["stoich"]
                    if j != "H2O"
                )
            )

        @self.Constraint(self.flowsheet().time)
        def flow_mass_h2o_key_eqn(blk, t):
            return self.flow_mass_h2o_key[t] == self.cv_aqueous.properties_out[
                t
            ].flow_mass - blk.flow_mass_h2o_key[t] * sum(
                blk.molality_aq_comp[t, i] * prop_aq.mw_comp[i]
                for i in prop_aq.key_set | prop_aq.aq_set
            )

        @self.Expression(
            self.flowsheet().time,
            prop_s.solid_dict,
            doc="Saturation index expression",
        )
        def saturation_index_expr(blk, t, comp):
            r = prop_s.solid_dict[comp]
            return (
                r["logK"]
                - sum(
                    r["stoich"][j] * blk.log10_molality_aq[t, j]
                    for j in r["stoich"]
                    if (j not in ["H2O", comp])
                )
                - sum(
                    r["stoich"][j] * blk.log10_act_coeff[t, j]
                    for j in r["stoich"]
                    if (j not in ["H2O", comp])
                )
            )

        # @self.Constraint(self.flowsheet().time, prop_s.solid_dict)
        # def saturation_index_eqn(blk, t, comp):
        #    return blk.saturation_index[t, comp] == blk.saturation_index_expr[t, comp]

        @self.Constraint(
            self.flowsheet().time,
            prop_s.solid_dict,
            doc="To precipitate or not to precipitate that is the question.",
        )
        def precipitate_eqn(blk, t, i):
            # return blk.molality_precipitate_comp[t, i] * blk.scale_factor_precipitate_molality == 0
            return (
                smooth_min(
                    blk.molality_precipitate_comp[t, i]
                    * blk.scale_factor_precipitate_molality,
                    -blk.saturation_index_expr[t, i]
                    * blk.scale_factor_precipitate_saturation_index,
                    eps=blk.smoothing_parameter_precipitation,
                )
                == 0
            )

        @self.Constraint(self.flowsheet().time, prop_s.solid_set)
        def outlet_precipitate_molality(blk, t, comp):
            return (
                blk.cv_precipitate.properties_out[t].flow_mol_comp[comp]
                == blk.molality_precipitate_comp[t, comp] * blk.flow_mass_h2o_key[t]
            )

        @self.Expression(
            self.flowsheet().time,
            prop_aq.key_set - {"H^+", "OH^-"},
            doc="Recovery fraction of key components (fraction in solid)",
        )
        def recovery(blk, t, comp):
            sol = -sum(
                blk.molality_precipitate_comp[t, i] * d["stoich"][comp]
                for i, d in prop_s.solid_dict.items()
                if comp in d["stoich"]
            )
            return sol / self.molality_key_comp[t, comp]

    def _add_control_volumes(self):
        # we have to attach this control volume to the model for the rest of
        # the steps to work
        cv_aqueous = self.cv_aqueous = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package_aqueous,
            property_package_args=self.config.property_package_args_aqueous,
        )
        cv_precipitate = self.cv_precipitate = ControlVolume0DBlock(
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
        return self.cv_aqueous, self.cv_precipitate

    def _add_variables_and_refs(self, prop_aq, prop_s):
        self.molality_key_comp = pyo.Var(
            self.flowsheet().time,
            prop_aq.key_set,
            initialize=1e-10,
            doc="Molality of key components in all forms",
            units=pyo.units.mol / pyo.units.kg,
        )

        self.molality_precipitate_comp = pyo.Var(
            self.flowsheet().time,
            prop_s.solid_set,
            initialize=1e-10,
            doc="Molality of key components in all forms",
            units=pyo.units.mol / pyo.units.kg,
        )

        self.saturation_index = pyo.Var(
            self.flowsheet().time,
            prop_s.solid_set,
            initialize=-10.0,
        )

        self.flow_mass_h2o_key = pyo.Var(
            self.flowsheet().time, initialize=1, units=pyo.units.kg / pyo.units.s
        )

        self.log10_molality_aq = pyo.Reference(
            self.cv_aqueous.properties_out[:].log10_molality_comp[:]
        )

        self.log10_act_coeff = pyo.Reference(
            self.cv_aqueous.properties_out[:].log10_act_coeff[:]
        )

        self.molality_aq_comp = pyo.Reference(
            self.cv_aqueous.properties_out[:].molality_comp[:]
        )

        self.pH = pyo.Reference(self.cv_aqueous.properties_out[:].pH)

    def _add_scaling_parameters(self):
        """Add some mutable parameters to the model to control scaling."""
        self.smoothing_parameter_precipitation = pyo.Param(
            initialize=1e-6,
            mutable=True,
            doc="Smoothing parameter for the precipitation complementarity constraint.",
        )

        self.scale_factor_precipitate_molality = pyo.Param(
            initialize=1e4,
            mutable=True,
            doc="Precipitate constraint saturation index scaling factor",
        )

        self.scale_factor_precipitate_saturation_index = pyo.Param(
            initialize=1,
            mutable=True,
            doc="Precipitate constraint saturation index scaling factor",
        )

        # self.scale_factor_pressure = pyo.Param(
        #     initialize=1e-5, mutable=True, doc="Pressure scaling factor"
        # )
        #
        # self.scale_factor_temperature = pyo.Param(
        #     initialize=1e-2, mutable=True, doc="Temperature scaling factor"
        # )

        self.scale_factor_mass_flow = pyo.Param(
            initialize=1, mutable=True, doc="Mass flow scaling factor"
        )

        self.scale_factor_molality = pyo.Param(
            initialize=1e4, mutable=True, doc="Molality scaling factor"
        )

    def calculate_scaling_factors(self):
        # calculate variable scaling
        prop_param_aq = self.config.property_package_aqueous
        prop_param_s = self.config.property_package_precipitate

        for t in self.flowsheet().time:
            # iscale.constraint_scaling_transform(
            #     self.temperature_aq_eqn[t], pyo.value(self.scale_factor_temperature)
            # )
            # iscale.constraint_scaling_transform(
            #     self.temperature_s_eqn[t], pyo.value(self.scale_factor_temperature)
            # )
            # iscale.constraint_scaling_transform(
            #     self.pressure_eqn[t], pyo.value(self.scale_factor_pressure)
            # )
            iscale.constraint_scaling_transform(
                self.flow_mass_h2o_key_eqn[t], pyo.value(self.scale_factor_mass_flow)
            )

            for comp in prop_param_aq.key_set:
                try:
                    iscale.constraint_scaling_transform(
                        self.molality_key_comp_eqn[t, comp],
                        pyo.value(self.scale_factor_molality),
                        overwrite=True,
                    )
                except:
                    pass
                try:
                    iscale.constraint_scaling_transform(
                        self.mass_balance[t, comp],
                        pyo.value(self.scale_factor_molality),
                        overwrite=True,
                    )
                except:
                    pass

            # iscale.set_scaling_factor(
            #     self.cv_aqueous.properties_in[t].temperature,
            #     pyo.value(self.scale_factor_temperature),
            # )
            # iscale.set_scaling_factor(
            #     self.cv_aqueous.properties_out[t].temperature,
            #     pyo.value(self.scale_factor_temperature),
            # )
            # iscale.set_scaling_factor(
            #     self.cv_aqueous.properties_in[t].pressure,
            #     pyo.value(self.scale_factor_pressure),
            # )
            # iscale.set_scaling_factor(
            #     self.cv_aqueous.properties_out[t].pressure,
            #     pyo.value(self.scale_factor_pressure),
            # )
            iscale.set_scaling_factor(
                self.cv_aqueous.properties_in[t].flow_mass,
                pyo.value(self.scale_factor_mass_flow),
            )
            iscale.set_scaling_factor(
                self.cv_aqueous.properties_out[t].flow_mass,
                pyo.value(self.scale_factor_mass_flow),
            )

            # iscale.set_scaling_factor(
            #     self.cv_precipitate.properties_in[t].temperature,
            #     pyo.value(self.scale_factor_temperature),
            # )
            # iscale.set_scaling_factor(
            #     self.cv_precipitate.properties_out[t].temperature,
            #     pyo.value(self.scale_factor_temperature),
            # )

            iscale.set_scaling_factor(
                self.flow_mass_h2o_key[t], pyo.value(self.scale_factor_mass_flow)
            )

            for comp in prop_param_s.solid_set:
                iscale.constraint_scaling_transform(
                    self.outlet_precipitate_molality[t, comp],
                    pyo.value(self.scale_factor_molality * self.scale_factor_mass_flow),
                )
                iscale.set_scaling_factor(
                    self.cv_precipitate.properties_in[t].flow_mol_comp[comp],
                    pyo.value(self.scale_factor_molality * self.scale_factor_mass_flow),
                )
                iscale.set_scaling_factor(
                    self.cv_precipitate.properties_out[t].flow_mol_comp[comp],
                    pyo.value(self.scale_factor_molality * self.scale_factor_mass_flow),
                )
                iscale.set_scaling_factor(
                    self.molality_precipitate_comp[t, comp],
                    pyo.value(self.scale_factor_molality * self.scale_factor_mass_flow),
                )

        for i in self.molality_key_comp:
            iscale.set_scaling_factor(self.molality_key_comp[i], 1e3)


def make_a_test_model():
    from idaes.core import FlowsheetBlock


    key_components = {
        "H^+",
        "Ce^3+",
        "Al^3+",
        "Fe^3+",
        # "Fe^2+",
        # "Ca^2+",
        # "Mg^2+",
        "C2O4^2-",
        # "NO3^-",
        # "SO4^2-",
        # "Cl^-",
    }

    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_aq = precip_prop.AqueousStateParameterBlock(
        key_components=key_components,
    )
    m.fs.properties_solid = precip_prop.PrecipitateStateParameterBlock(
        key_components=key_components,
    )

    m.fs.unit = Precipitator(
        property_package_aqueous=m.fs.properties_aq,
        property_package_precipitate=m.fs.properties_solid,
    )
    return m


def main():
    m = make_a_test_model()
    # m.fs.unit.cv_precipitate.properties_in[0].temperature.fix(300)
    m.fs.unit.cv_precipitate.properties_in[0].flow_mol_comp.fix(0)

    # m.fs.unit.cv_aqueous.properties_in[0].temperature.fix(300)
    # m.fs.unit.cv_aqueous.properties_in[0].pressure.fix(101325)
    m.fs.unit.cv_aqueous.properties_in[0].flow_mass.fix(1)
    m.fs.unit.cv_aqueous.properties_in[0].log10_molality_comp.fix(-20)
    m.fs.unit.cv_aqueous.properties_in[0].log10_molality_comp["H^+"].fix(-1)
    m.fs.unit.cv_aqueous.properties_in[0].log10_molality_comp["OH^-"].fix(-15)
    m.fs.unit.cv_aqueous.properties_in[0].log10_molality_comp["Ce^3+"].fix(-4)
    m.fs.unit.cv_aqueous.properties_in[0].log10_molality_comp["Al^3+"].fix(-10)
    m.fs.unit.cv_aqueous.properties_in[0].log10_molality_comp["Fe^3+"].fix(-10)
    # m.fs.unit.cv_aqueous.properties_in[0].log10_molality_comp["Fe^2+"].fix(-10)
    # m.fs.unit.cv_aqueous.properties_in[0].log10_molality_comp["Ca^2+"].fix(-10)
    # m.fs.unit.cv_aqueous.properties_in[0].log10_molality_comp["Mg^2+"].fix(-10)
    m.fs.unit.cv_aqueous.properties_in[0].log10_molality_comp["C2O4^2-"].fix(-4)
    # m.fs.unit.cv_aqueous.properties_in[0].log10_molality_comp["NO3^-"].fix(-10)
    # m.fs.unit.cv_aqueous.properties_in[0].log10_molality_comp["Cl^-"].fix(-10)
    # m.fs.unit.cv_aqueous.properties_in[0].log10_molality_comp["SO4^2-"].fix(-10)

    m.fs.unit.calculate_scaling_factors()

    assert_units_consistent(m)

    solver_obj = pyo.SolverFactory(
        "ipopt",
        options={
            "nlp_scaling_method": "user-scaling",
            "linear_solver": "ma57",
            "OF_ma57_automatic_scaling": "yes",
            "ma57_pivtol": 1e-5,
            "ma57_pivtolmax": 0.1,
            "tol": 1e-8,
            # "halt_on_ampl_error": "yes",
            "max_iter": 50,
        },
    )
    solver_obj.solve(
        m,
        tee=True,
    )
    return m, solver_obj


if __name__ == "__main__":
    m, solver_obj = main()

    m.fs.unit.pH.fix(2)
    m.fs.unit.cv_aqueous.properties_in[0].log10_molality_comp["H^+"].unfix()
    solver_obj.solve(m, tee=True)

    m.fs.unit.cv_aqueous.properties_out[0].molality_comp.display()
    m.fs.unit.cv_precipitate.properties_out[0].flow_mol_comp.display()
    m.fs.unit.molality_precipitate_comp.display()
    m.fs.unit.pH.display()
