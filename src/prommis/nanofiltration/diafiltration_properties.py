#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Property package for sieving coefficient model for lithium and cobalt.
"""

from pyomo.environ import Param, Var, units

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


class _StateBlock(StateBlock):
    """
    Base class for state block
    """

    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.
        """
        fix_state_vars(self)


@declare_process_block_class("LiCoStateBlock", block_class=_StateBlock)
class LiCoStateBlock1Data(StateBlockData):
    """
    State block for the lithium cobalt property package
    """

    def build(self):
        super().build()

        self.flow_vol = Var(
            units=units.m**3 / units.hour,
            bounds=(1e-8, None),
        )
        self.conc_mass_solute = Var(
            ["Li", "Co"],
            units=units.kg / units.m**3,
            bounds=(1e-8, None),
        )

    def get_material_flow_terms(self, p, j):
        if j == "solvent":
            # Assume constant density of pure water
            return self.flow_vol * self.params.dens_H2O
        else:
            return self.flow_vol * self.conc_mass_solute[j]

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mass

    def define_state_vars(self):
        return {
            "flow_vol": self.flow_vol,
            "conc_mass_solute": self.conc_mass_solute,
        }


@declare_process_block_class("LiCoParameters")
class LiCoParameterData(PhysicalParameterBlock):
    """
    Parameter block for the lithium cobalt property package
    """

    def build(self):
        super().build()

        self.phase1 = Phase()

        self.solvent = Component()
        self.Li = Component()
        self.Co = Component()

        self.dens_H2O = Param(
            default=1000,
            units=units.kg / units.m**3,
        )

        self._state_block_class = LiCoStateBlock

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
