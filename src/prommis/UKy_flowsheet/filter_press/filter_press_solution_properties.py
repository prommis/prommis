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
Initial property package for liquids in the filter press.

Authors: Marcus Holly
"""
from pyomo.environ import (
    Param,
    units,
    Var,
)

from idaes.core import (
    declare_process_block_class,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    Component,
    Phase,
    MaterialFlowBasis,
    MaterialBalanceType,
    EnergyBalanceType,
)
from idaes.core.util.exceptions import BurntToast
from idaes.core.util.initialization import fix_state_vars


# -----------------------------------------------------------------------------
# Leach solution property package
@declare_process_block_class("FilterPressParameters")
class FilterPressParameterData(PhysicalParameterBlock):
    def build(self):
        super().build()

        self.liquid = Phase()

        # Solvent
        self.H2O = Component()

        # Liquids - placeholder components for now
        self.TDS = Component()

        self.mw = Param(
            self.component_list,
            units=units.kg / units.mol,
            initialize={
                "H2O": 18e-3,
                "TDS": 55e-3,
            },
        )

        self._state_block_class = FilterPressStateBlock

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


class _FilterPressStateBlock(StateBlock):
    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        # Fix state variables
        fix_state_vars(self)


@declare_process_block_class(
    "FilterPressStateBlock", block_class=_FilterPressStateBlock
)
class FilterPressStateBlockData(StateBlockData):
    def build(self):
        super().build()

        self.flow_vol = Var(
            units=units.L / units.hour,
            initialize=100,
            bounds=(1e-8, None),
        )
        self.conc_mole = Var(
            self.params.component_list,
            units=units.mol / units.L,
            initialize=1e-5,
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

    def get_material_flow_terms(self, p, j):
        # Note conversion to kg/hour
        if j == "H2O":
            # Assume constant density of 1 kg/L
            return self.flow_vol * (1 * units.kg / units.L) / self.params.mw[j]
        elif j in self.params.component_list:
            # Need to convert from moles to mass
            return self.flow_vol * self.conc_mole[j]
        else:
            raise BurntToast()

    def default_material_balance_type(self):
        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        return {
            "flow_vol": self.flow_vol,
            "conc_mole": self.conc_mole,
            "temperature": self.temperature,
            "pressure": self.pressure,
        }
