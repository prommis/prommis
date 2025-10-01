#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

#################################################################################
# This model was derived from:
# https://github.com/kurbansitterley/watertap/blob/ix_reorg/watertap/unit_models/ion_exchange/ion_exchange_multi_component.py

# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

from copy import deepcopy

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.common.config import ConfigValue, In

# Import IDAES cores
from idaes.core import declare_process_block_class
from idaes.core.util.math import smooth_min, smooth_max
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import InitializationError, ConfigurationError

from prommis.ion_exchange.ion_exchange_multicomponent_base import (
    IonExchangeBaseData,
    IonExchangeType,
)


""" This model contains the Clark model equations using Freundlich
equilibrium equations. This is a modified version of multicomponent
model from the WaterTAP developers:
https://github.com/kurbansitterley/watertap/blob/ix_reorg/watertap/unit_models/ion_exchange/ion_exchange_multi_component.py

REFERENCES

[1] Vassilis J. Inglezakis and Stavros G. Poulopoulos Adsorption, Ion
Exchange and Catalysis: Design of Operations and Environmental
Applications (2006).  doi.org/10.1016/B978-0-444-52783-7.X5000-9

[2] Croll, H. C., Adelman, M. J., Chow, S. J., Schwab, K. J., Capelle,
R., Oppenheimer, J., & Jacangelo, J. G. (2023).  Fundamental kinetic
constants for breakthrough of per- and polyfluoroalkyl substances at
varying empty bed contact times: Theoretical analysis and pilot scale
demonstration.  Chemical Engineering Journal,
464. doi:10.1016/j.cej.2023.142587

[3] Clark, R. M. (1987).  Evaluating the cost and performance of
field-scale granular activated carbon systems.  Environ Sci Technol,
21(6), 573-580. doi:10.1021/es00160a008


modified by: Soraya Rawlings

"""

__author__ = "Kurban Sitterley"


@declare_process_block_class("IonExchangeMultiComp")
class IonExchangeMultiCompData(IonExchangeBaseData):
    """
    Ion exchange multi component Clark model
    """

    CONFIG = IonExchangeBaseData.CONFIG()

    CONFIG.declare(
        "reactive_ions",
        ConfigValue(
            default=list(),
            domain=list,
            description="Designates other reactive species",
        ),
    )

    def build(self):
        super().build()

        prop_in = self.process_flow.properties_in[0]

        if self.config.regenerant != "single_use":
            # [WIP: Add regen terms when regeneration is needed]
            regen = self.regeneration_stream[0]

        comps = self.config.property_package.component_list
        target_component = self.config.target_component

        target_component = self.config.target_component
        reactive_ions = self.config.reactive_ions

        assert target_component in reactive_ions
        self.target_component_set = pyo.Set(initialize=[target_component])

        # [ESR WIP: Define if this is needed in future
        # cases. Commented for now since it is not used in our cases]
        # inerts = comps - self.target_component_set

        self.reactive_ion_set = pyo.Set(initialize=reactive_ions)

        self.inert_set = pyo.Set(initialize=comps - self.reactive_ion_set)

        if len(self.target_component_set) > 1:
            raise ConfigurationError(
                f"IonExchange0D can only accept a single target ion but {len(self.target_component_set)} were provided."
            )
        if self.config.property_package.charge_comp[target_component].value > 0:
            self.ion_exchange_type = IonExchangeType.cation
        elif self.config.property_package.charge_comp[target_component].value < 0:
            self.ion_exchange_type = IonExchangeType.anion
        else:
            raise ConfigurationError("Target ion must have non-zero charge.")

        # [ESR update: Change the inerts for the ones in the set.]
        for j in self.inert_set:
            self.process_flow.mass_transfer_term[:, "Liq", j].fix(0)

            if self.config.regenerant != "single_use":
                # [ESR WIP: Uncomment regen terms]
                regen.get_material_flow_terms("Liq", j).fix(0)

        self.num_traps = int(self.config.number_traps)

        # trap_disc is a discretization index/parameter that defines
        # how the range between c_trap_min and c_norm is broken up.
        self.trap_disc = range(self.num_traps + 1)

        self.trap_index = self.trap_disc[1:]

        self.eps = pyo.Param(initialize=1e-4, mutable=True)

        self.c_trap_min = {}
        for i in self.reactive_ion_set:
            self.c_trap_min[i] = float(self.config.c_trap_min)

        # [ESR WIP: Bring breathrough time and bv here from base model
        # since these variables depend on each ion. Also, add a
        # constraint to calculate the maximum breakthrough time. For
        # now, this maximum time is for the target component. Think if
        # there is a better way to do this.]
        self.breakthrough_time = pyo.Var(
            self.reactive_ion_set,
            initialize=1e5,  # DOW, ~7 weeks max breakthru time
            bounds=(0, None),
            # domain=NonNegativeReals,
            units=pyo.units.s,
            doc="Breakthrough time",
        )
        self.bv = pyo.Var(  # BV
            self.reactive_ion_set,
            initialize=1e5,
            bounds=(0, None),
            # domain=NonNegativeReals,
            units=pyo.units.dimensionless,
            doc="Bed volumes of feed at breakthru concentration",
        )

        @self.Constraint()
        def eq_target_breakthru_time(b):
            return (
                self.target_breakthrough_time
                == self.breakthrough_time[target_component]
            )

        # Add variables for dimensionless numbers
        self.N_Sc = pyo.Var(
            self.reactive_ion_set,
            initialize=700,
            units=pyo.units.dimensionless,
            doc="Schmidt number",
        )

        self.N_Sh = pyo.Var(
            self.reactive_ion_set,
            initialize=30,
            units=pyo.units.dimensionless,
            doc="Sherwood number",
        )

        # Add relevant variables for IX column and trapezoidal rule
        self.c_norm = pyo.Var(
            self.reactive_ion_set,
            initialize=0.5,
            bounds=(0, 1),
            units=pyo.units.dimensionless,
            doc="Dimensionless (relative) concentration [Ct/C0] of target ion",
        )

        self.c_traps = pyo.Var(
            self.reactive_ion_set,
            self.trap_disc,
            initialize=0.5,
            bounds=(0, 1),
            units=pyo.units.dimensionless,
            doc="Normalized breakthrough concentrations for estimating area under breakthrough curve",
        )

        self.tb_traps = pyo.Var(
            self.reactive_ion_set,
            self.trap_disc,
            initialize=1,
            bounds=(0, None),
            units=pyo.units.second,
            doc="Breakthrough times for estimating area under breakthrough curve",
        )

        self.traps = pyo.Var(
            self.reactive_ion_set,
            self.trap_index,
            initialize=0.01,
            bounds=(0, 1),
            units=pyo.units.dimensionless,
            doc="Trapezoid areas for estimating area under breakthrough curve",
        )

        for c in self.reactive_ion_set:
            self.c_traps[(c, 0)].fix(0)
            self.tb_traps[(c, 0)].fix(0)

        self.c_norm_avg = pyo.Var(
            self.reactive_ion_set,
            initialize=0.25,
            bounds=(0, 2),
            units=pyo.units.dimensionless,
            doc="Sum of trapezoid areas",
        )

        self.freundlich_n = pyo.Var(
            self.reactive_ion_set,
            initialize=1.5,
            bounds=(0, None),
            units=pyo.units.dimensionless,
            doc="Freundlich isotherm exponent",
        )

        self.mass_transfer_coeff = pyo.Var(  # k_T
            self.reactive_ion_set,
            initialize=0.001,
            units=pyo.units.s**-1,
            bounds=(0, None),
            doc="Mass transfer coefficient for Clark model (kT)",
        )

        self.bv_50 = pyo.Var(  # BV_50
            self.reactive_ion_set,
            initialize=2e5,
            bounds=(0, None),
            units=pyo.units.dimensionless,
            doc="Bed volumes of feed at 50 percent breakthrough",
        )

        @self.Expression(self.reactive_ion_set, doc="Breakthrough concentration")
        def c_breakthru(b, j):
            return b.c_norm[j] * prop_in.conc_mass_phase_comp["Liq", j]

        # Add constraints to calculate dimensionless numbers
        @self.Constraint(self.reactive_ion_set, doc="Schmidt number")
        def eq_Sc(b, j):  # Eq. 3.359, ref[1] Inglezakis + Poulopoulos
            return (
                b.N_Sc[j]
                == prop_in.visc_k_phase["Liq"] / prop_in.diffus_phase_comp["Liq", j]
            )

        @self.Constraint(self.reactive_ion_set, doc="Sherwood number")
        def eq_Sh(b, j):  # Eq. 3.346, ref[1] Inglezakis + Poulopoulos
            return (
                b.N_Sh[j]
                == b.Sh_A
                * b.bed_porosity**b.Sh_exp_A
                * b.N_Re**b.Sh_exp_B
                * b.N_Sc[j] ** b.Sh_exp_C
            )

        # [ESR updates: Add an if statement to consider this
        # constraint ONLY for the target component since it was adding
        # 1 degree of freedom when solving for multiple ions. [TODO:
        # Review this change once we add the regeneration stage.]

        if self.config.regenerant != "single_use":
            # [ESR WIP: Uncomment regen terms]
            @self.Constraint(
                self.reactive_ion_set, doc="Mass transfer for regeneration stream"
            )
            def eq_mass_transfer_regen(b, j):
                return (
                    regen.get_material_flow_terms("Liq", j)
                    == -b.process_flow.mass_transfer_term[0, "Liq", j]
                )

        # [ESR updates: Add multicomponent set for breakthrough_time
        # and bv]
        @self.Constraint(self.reactive_ion_set, doc="Bed volumes at breakthrough")
        def eq_bv(b, r):
            return b.breakthrough_time[r] * b.loading_rate == b.bv[r] * b.bed_depth

        @self.Constraint(
            self.reactive_ion_set, doc="Clark equation with fundamental constants"
        )  # ref[2], Croll et al (2023), Eq.9
        def eq_clark(b, j):
            left_side = pyo.exp(
                (
                    (b.mass_transfer_coeff[j] * b.bed_depth * (b.freundlich_n[j] - 1))
                    / (b.bv_50[j] * b.loading_rate)
                )
                * (b.bv_50[j] - b.bv[j])
            )
            right_side = ((1 / b.c_norm[j]) ** (b.freundlich_n[j] - 1) - 1) / (
                2 ** (b.freundlich_n[j] - 1) - 1
            )
            return left_side - right_side == 0

        # [ESR notes: This is the concentration value at the kth
        # “trapezoid” boundary for ion j.]
        @self.Constraint(
            self.reactive_ion_set,
            self.trap_index,
            doc="Evenly spaced c_norm for trapezoids",
        )
        def eq_c_traps(b, j, k):
            if k == max(b.trap_index):
                return b.c_traps[j, k] == b.c_norm[j]
            else:
                # [ESR updates: Reformulated to avoid denominators.]
                return b.c_traps[j, k] * (b.num_traps - 1) == (
                    b.c_trap_min[j] * (b.num_traps - 1)
                    + (b.trap_disc[k] - 1) * (b.c_norm[j] - b.c_trap_min[j])
                )

        # [ESR updates: Add bv_traps as an expression. This term
        # corresponds to the bed volume corresponding to the kth
        # discretization point for ion j.]
        @self.Expression(self.reactive_ion_set, self.trap_disc, doc="BV for trapezoids")
        def bv_traps(b, j, k):
            return (b.tb_traps[j, k] * b.loading_rate) / b.bed_depth

        @self.Constraint(
            self.reactive_ion_set,
            self.trap_index,
            doc="Breakthru time calc for trapezoids",
        )
        def eq_tb_traps(b, j, k):
            if k == max(self.trap_index):
                return b.tb_traps[j, k] == b.breakthrough_time[j]
            else:
                left_side = pyo.exp(
                    (
                        (
                            b.mass_transfer_coeff[j]
                            * b.bed_depth
                            * (b.freundlich_n[j] - 1)
                        )
                        / (b.bv_50[j] * b.loading_rate)
                    )
                    * (b.bv_50[j] - b.bv_traps[j, k])
                )
                right_side = ((1 / b.c_traps[j, k]) ** (b.freundlich_n[j] - 1) - 1) / (
                    2 ** (b.freundlich_n[j] - 1) - 1
                )
                return left_side - right_side == 0

        # [ESR WIP: Try to make changes in the normalization term in
        # eq_traps to see if that corrects the breakthrough time
        # differences. This normalization helps to express each
        # trapezoid’s width as a fraction of the total interval. This
        # normalization helps ensure that the integrated output (when
        # all trapezoidal areas are summed) gives a proper,
        # dimensionless measure of the overall breakthrough
        # behavior. NOTE: Re-write to avoid denominators]
        @self.Constraint(
            self.reactive_ion_set, self.trap_index, doc="Area of trapezoids"
        )
        def eq_traps(b, j, k):
            # Norm term is the normalization term

            # Original normalization value
            norm_term = b.tb_traps[j, self.num_traps]

            # # Use total interval value for normalization
            # norm_term = b.tb_traps[j, self.num_traps] - b.tb_traps[j, 1]

            # No normalization
            # norm_term = 1

            return b.traps[j, k] * norm_term == (
                # width of the trapezoid
                b.tb_traps[j, k]
                - b.tb_traps[j, k - 1]
            ) * (
                # average height between the two points
                (b.c_traps[j, k] + b.c_traps[j, k - 1])
                / 2
            )

        @self.Constraint(
            self.reactive_ion_set, doc="Average relative effluent concentration"
        )
        def eq_c_norm_avg(b, j):
            return b.c_norm_avg[j] == sum(b.traps[j, k] for k in b.trap_index)

        # [ESR WIP: Separate term]
        @self.Expression(
            self.reactive_ion_set, doc="Total mass of ion removed from the liquid phase"
        )
        def total_mass_removed(b, j):
            return (1 - b.c_norm_avg[j]) * prop_in.get_material_flow_terms("Liq", j)

        @self.Constraint(
            # self.target_component_set,
            self.reactive_ion_set,
            doc="CV mass transfer term",
        )
        def eq_mass_transfer_reactive_ions(b, j):
            # return (1 - b.c_norm_avg[j]) * prop_in.get_material_flow_terms(
            #     "Liq", j
            # ) == -b.process_flow.mass_transfer_term[0, "Liq", j]
            return (
                b.total_mass_removed[j]
                == -b.process_flow.mass_transfer_term[0, "Liq", j]
            )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # [ESR updates: Add new variables from base model.]
        if iscale.get_scaling_factor(self.breakthrough_time) is None:
            iscale.set_scaling_factor(self.breakthrough_time, 1e-4)

        if iscale.get_scaling_factor(self.bv) is None:
            iscale.set_scaling_factor(self.bv, 1e-1)

        if iscale.get_scaling_factor(self.freundlich_n) is None:
            iscale.set_scaling_factor(self.freundlich_n, 1)

        if iscale.get_scaling_factor(self.mass_transfer_coeff) is None:
            iscale.set_scaling_factor(self.mass_transfer_coeff, 1e2)

        if iscale.get_scaling_factor(self.bv_50) is None:
            iscale.set_scaling_factor(self.bv_50, 1e-3)

        if iscale.get_scaling_factor(self.tb_traps) is None:
            sf = iscale.get_scaling_factor(self.breakthrough_time)
            iscale.set_scaling_factor(self.tb_traps, sf)

        if iscale.get_scaling_factor(self.c_traps) is None:
            iscale.set_scaling_factor(self.c_traps, 1)

        if iscale.get_scaling_factor(self.traps) is None:
            iscale.set_scaling_factor(self.traps, 1e3)

        if iscale.get_scaling_factor(self.c_norm_avg) is None:
            iscale.set_scaling_factor(self.c_norm_avg, 1e2)

        for ind, c in self.eq_clark.items():
            if iscale.get_scaling_factor(c) is None:
                iscale.constraint_scaling_transform(c, 1e-2)

        for ind, c in self.eq_traps.items():
            if iscale.get_scaling_factor(c) is None:
                iscale.constraint_scaling_transform(c, 1e2)
