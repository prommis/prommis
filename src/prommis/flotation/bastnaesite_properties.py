#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Dry-solids property package for bastnaesite solids."""

from pyomo.common.config import Bool, ConfigValue
from pyomo.environ import Expression, Param, PositiveReals, Var, units

from idaes.core import (
    Component,
    MaterialFlowBasis,
    Phase,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    declare_process_block_class,
)
from idaes.core.scaling import CustomScalerBase
from idaes.core.util.initialization import fix_state_vars

__author__ = "Daison Yancy Caballero"

COMPONENTS = ("REO", "CaO", "BaO", "SrO", "inert_gangue")
# Source: Wang et al. report baseline bastnaesite ore BWI as 6.96 kWh/t
# before CO2 treatment in "A hybrid ex-situ, in-situ, and multi-step CO2
# treatment processing for enhancing the grindability of bastnaesite ore":
# https://www.sciencedirect.com/science/article/abs/pii/S0892687525007083
# kWh/t is equivalent to W*h/kg, the model's native units.
BASTNAESITE_BOND_WORK_INDEX = 6.96
# Per-component true (material) densities [kg/m^3] of each component's
# representative host mineral: quartz (inert_gangue), calcite (CaO), celestine
# (SrO), barite (BaO), and bastnaesite-(Ce) (REO). These are size-independent
# crystal densities (not bulk/apparent packing densities), used only to convert
# dry-solid mass flow to volumetric flow. Source: Handbook of Mineralogy
# (Anthony, Bideaux, Bladh & Nichols, Mineral Data Publishing, 2001-2005).
SOLID_DENSITIES = {
    "REO": 5000.0,
    "CaO": 2710.0,
    "BaO": 4480.0,
    "SrO": 3960.0,
    "inert_gangue": 2650.0,
}
REFERENCE_SCALING_FACTORS = {component: 1.0 for component in COMPONENTS}


class BastnaesitePropertiesScaler(CustomScalerBase):
    """Scaler for the bastnaesite dry-solids state block."""

    CONFIG = CustomScalerBase.CONFIG

    DEFAULT_SCALING_FACTORS = {
        f"flow_mass_comp[{component}]": REFERENCE_SCALING_FACTORS[component]
        for component in COMPONENTS
    }

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        for component in model.params.component_list:
            self.set_variable_scaling_factor(
                model.flow_mass_comp[component],
                REFERENCE_SCALING_FACTORS[component],
                overwrite=overwrite,
            )
        if model.params.config.has_particle_size_distribution:
            self.set_variable_scaling_factor(
                model.particle_size_median,
                1e-5,
                overwrite=overwrite,
            )
            self.set_variable_scaling_factor(
                model.particle_size_width,
                1.0,
                overwrite=overwrite,
            )


@declare_process_block_class("BastnaesiteParameters")
class BastnaesiteParameterData(PhysicalParameterBlock):
    """Parameter block for oxide-equivalent dry solids."""

    CONFIG = PhysicalParameterBlock.CONFIG()
    CONFIG.declare(
        "has_particle_size_distribution",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Include median-width particle size distribution state",
        ),
    )

    def build(self):
        super().build()

        self.solid = Phase()
        for component in COMPONENTS:
            setattr(self, component, Component())

        self.rho_mass_comp = Param(
            self.component_list,
            initialize=SOLID_DENSITIES,
            domain=PositiveReals,
            mutable=True,
            units=units.kg / units.m**3,
            doc="Per-component solid mass density",
        )

        if self.config.has_particle_size_distribution:
            self.bond_work_index = Param(
                initialize=BASTNAESITE_BOND_WORK_INDEX,
                mutable=True,
                units=units.W * units.hour / units.kg,
                doc="Bastnaesite Bond work index",
            )

        self._state_block_class = BastnaesiteStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "flow_mass_comp": {"method": None},
                "flow_mass": {"method": None},
                "flow_vol_phase": {"method": None},
                "mass_frac_comp": {"method": "_mass_frac_comp"},
            }
        )
        obj.define_custom_properties(
            {
                "particle_size_median": {"method": None},
                "particle_size_width": {"method": None},
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


class _BastnaesiteStateBlock(StateBlock):
    default_scaler = BastnaesitePropertiesScaler

    def fix_initialization_states(self):
        fix_state_vars(self)


@declare_process_block_class(
    "BastnaesiteStateBlock",
    block_class=_BastnaesiteStateBlock,
)
class BastnaesiteStateBlockData(StateBlockData):
    """State block with direct component dry-solid mass flows."""

    default_scaler = BastnaesitePropertiesScaler

    def build(self):
        super().build()

        self.flow_mass_comp = Var(
            self.params.component_list,
            initialize=1.0,
            units=units.kg / units.hour,
            bounds=(0, None),
            doc="Component dry-solid mass flow",
        )
        if self.params.config.has_particle_size_distribution:
            self.particle_size_median = Var(
                initialize=80000,
                units=units.um,
                bounds=(10, None),
                doc="Median particle size",
            )
            self.particle_size_width = Var(
                initialize=1.5,
                units=units.dimensionless,
                bounds=(1e-6, None),
                doc="Particle size distribution width parameter",
            )

        self.flow_mass = Expression(
            expr=sum(
                self.flow_mass_comp[component]
                for component in self.params.component_list
            ),
            doc="Total dry-solid mass flow",
        )

        def flow_vol_phase_rule(b, _phase):
            return sum(
                b.flow_mass_comp[component] / b.params.rho_mass_comp[component]
                for component in b.params.component_list
            )

        self.flow_vol_phase = Expression(
            self.params.phase_list,
            rule=flow_vol_phase_rule,
            doc="Dry-solid volumetric flow indexed by phase",
        )

    def _mass_frac_comp(self):
        # Note: this Expression is undefined for zero-flow streams (division by
        # zero). It is used only for reporting; callers that may encounter
        # zero-flow streams should guard before calling value().
        def mass_frac_comp_rule(b, j):
            return b.flow_mass_comp[j] / b.flow_mass

        self.mass_frac_comp = Expression(
            self.params.component_list,
            rule=mass_frac_comp_rule,
            doc="Component dry-solid mass fraction (undefined for zero flow)",
        )

    def get_material_flow_terms(self, p, j):
        return self.flow_mass_comp[j]

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mass

    def define_state_vars(self):
        state_vars = {"flow_mass_comp": self.flow_mass_comp}
        if self.params.config.has_particle_size_distribution:
            state_vars.update(
                {
                    "particle_size_median": self.particle_size_median,
                    "particle_size_width": self.particle_size_width,
                }
            )
        return state_vars
