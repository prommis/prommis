#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Property package for rare earth permanent magnet (REPM) materials.

Authors: Brandon Paul
"""

from pyomo.environ import Constraint, Param, Var, units

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


# -----------------------------------------------------------------------------
# Leach solution property package
@declare_process_block_class("REPMParameters")
class REPMParameterData(PhysicalParameterBlock):
    """
    Solid phase property package for rare earth permanent magnets, specifically
    recycling neodymium materials from hard drive disks (HDDs).


    Includes the following components: Nd, Nd2Fe14B

    """

    def build(self):
        super().build()

        self.solid = Phase()

        self.Nd = Component()
        # self.NdH2 = Component()
        self.Nd2Fe14B = Component()
        # self.Nd2Fe14BH2 = Component()

        self.mw = Param(
            self.component_list,
            units=units.kg / units.mol,
            initialize={
                "Nd": (144.242 * 1) * 1e-3,
                # "NdH2": (144.242 * 1 + 1.00794 * 2) * 1e-3,
                "Nd2Fe14B": (144.242 * 2 + 55.845 * 14 + 10.811 * 1) * 1e-3,
                # "Nd2Fe14BH2": (144.242 * 2 + 55.845 * 14 + 10.811 * 1 + 1.00794 * 2) * 1e-3,
            },
        )

        # TODO do i need this?
        self.dens_mass = Param(
            units=units.kg / units.litre,
            initialize=2.4,
            mutable=True,
        )

        self._state_block_class = REPMStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": units.s,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


class _REPMStateBlock(StateBlock):
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


@declare_process_block_class("REPMStateBlock", block_class=_REPMStateBlock)
class REPMStateBlockData(StateBlockData):
    """
    State block for rare earth permanent magnent (REPM) materials.

    """

    def build(self):
        super().build()

        self.flow_mass = Var(
            units=units.kg / units.s,
            bounds=(1e-8, None),
        )
        self.mass_frac_comp = Var(
            self.params.component_list,
            units=units.kg / units.kg,
            bounds=(1e-8, None),
        )

        if not self.config.defined_state:
            self.sum_mass_frac = Constraint(
                expr=1 == sum(self.mass_frac_comp[j] for j in self.component_list)
            )

    def get_material_flow_terms(self, p, j):
        return self.flow_mass * self.mass_frac_comp[j] / self.params.mw[j]

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        return {
            "flow_mass": self.flow_mass,
            "mass_frac_comp": self.mass_frac_comp,
        }
