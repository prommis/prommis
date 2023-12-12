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
Initial property package for solids in the filter press.

Authors: Marcus Holly
"""

from idaes.core import (
    Component,
    EnergyBalanceType,
    MaterialBalanceType,
    MaterialFlowBasis,
    Phase,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    declare_process_block_class,
)
from idaes.core.util.initialization import fix_state_vars
from pyomo.environ import Constraint, Param, Var, units


# -----------------------------------------------------------------------------
# Filter press solids property package
@declare_process_block_class("FilterPressSolidsParameters")
class FilterPressSolidsParameterData(PhysicalParameterBlock):
    def build(self):
        super().build()

        self.solid = Phase()

        # Solids - placeholder components for now
        self.solid1 = Component()
        self.solid2 = Component()

        self.mw = Param(
            self.component_list,
            units=units.kg / units.mol,
            initialize={
                "solid1": 100e-3,
                "solid2": 50e-3,
            },
        )

        self._state_block_class = FilterPressSolidsStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "temperature": {"method": None},
                "pressure": {"method": None},
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


class _FilterPressSolidsStateBlock(StateBlock):
    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        # Fix state variables
        fix_state_vars(self)

        if not self.config.defined_state:
            self.sum_mass_frac.deactivate()


@declare_process_block_class(
    "FilterPressSolidsStateBlock", block_class=_FilterPressSolidsStateBlock
)
class FilterPressSolidsStateBlockData(StateBlockData):
    def build(self):
        super().build()

        self.flow_mass = Var(
            units=units.kg / units.hour,
            bounds=(1e-8, None),
        )
        self.mass_frac_comp = Var(
            self.params.component_list,
            units=units.kg / units.kg,
            bounds=(1e-8, None),
        )
        self.temperature = Var(
            initialize=298.15,
            bounds=(273.15, 1000),
            units=units.K,
            doc="Temperature",
        )
        self.pressure = Var(
            initialize=101325,
            bounds=(1e3, 5e7),
            units=units.Pa,
            doc="Pressure",
        )

        if not self.config.defined_state:
            self.sum_mass_frac = Constraint(
                expr=1 == sum(self.mass_frac_comp[j] for j in self.component_list)
            )

    def get_material_flow_terms(self, p, j):
        return self.flow_mass * self.mass_frac_comp[j] / self.params.mw[j]

    def default_material_balance_type(self):
        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        return {
            "flow_mass": self.flow_mass,
            "mass_frac_comp": self.mass_frac_comp,
            "temperature": self.temperature,
            "pressure": self.pressure,
        }
