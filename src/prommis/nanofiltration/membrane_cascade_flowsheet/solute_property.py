#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Generic multi-solute property package for membrane filtration flowsheet."""

from pyomo.common.config import ConfigValue
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
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.initialization import fix_state_vars


class _StateBlock(StateBlock):
    def fix_initialization_states(self):
        """
        Fix state variables for state blocks.

        Returns:
            None
        """
        fix_state_vars(self)

    def initialization_routine(self):
        pass


@declare_process_block_class("SoluteStateBlock", block_class=_StateBlock)
class SoluteStateBlock1Data(StateBlockData):
    """Defines thermophysical property variables and constraints."""

    def build(self):
        """Construct model variables/expressions/constraints."""
        super().build()

        self.flow_vol = Var(units=units.m**3 / units.hour, bounds=(1e-8, None))

        # remove solvent from component list
        solutes = [i for i in self.component_list if i != "solvent"]

        self.flow_mass_solute = Var(
            solutes, units=units.kg / units.hour, bounds=(1e-8, None)
        )

    def get_material_flow_terms(self, p, j):
        """Get mass flow rate."""
        if j == "solvent":
            # Assume constant density of pure water
            return self.flow_vol * self.params.dens_H2O
        else:
            return self.flow_mass_solute[j]

    def get_material_flow_basis(self):
        """Get material flow basis."""
        return MaterialFlowBasis.mass

    def define_state_vars(self):
        """
        Create dict of state variables.

        Returns: dict
            a dict of state variables
        """
        return {"flow_vol": self.flow_vol, "flow_mass_solute": self.flow_mass_solute}


@declare_process_block_class("SoluteParameters")
class SoluteParameterData(PhysicalParameterBlock):
    """Defines global parameters and components of solute property package."""

    CONFIG = PhysicalParameterBlock.CONFIG()

    CONFIG.declare(
        "solutes",
        ConfigValue(domain=list, doc="List of solutes present in the membrane system"),
    )

    def build(self):
        """Construct modeling components."""
        super().build()

        # Check that list of solutes exists
        if self.config.solutes is None:
            raise ConfigurationError(
                "Solutes property package requires " + "a list of solutes."
            )

        self.phase1 = Phase()

        self.solvent = Component()

        # set solute components
        for sol in self.config.solutes:
            setattr(self, f"{sol}", Component())

        self.dens_H2O = Param(default=1000, units=units.kg / units.m**3)

        self._state_block_class = SoluteStateBlock  # noqa

    @classmethod
    def define_metadata(cls, obj):
        """Define class metadata."""
        obj.add_default_units(
            {
                "time": units.hour,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )
