#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Property package for the multi-salt diafiltration membrane.

Author: Molly Dougher
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


@declare_process_block_class("SoluteParameter")
class SoluteParameterData(PhysicalParameterBlock):
    """
    Property Package for the multi-salt diafiltration membrane.

    Currently includes the following solutes:
        Li+ (lithium ion)
        Co2+ (cobalt ion)
        Cl- (chlorine ion)
    """

    def build(self):
        super().build()

        self.liquid = Phase()

        # add cations
        self.Li = Component()
        self.Co = Component()

        # add anions
        self.Cl = Component()

        # add valence
        self.charge = Param(
            self.component_list,
            units=units.dimensionless,
            initialize={
                "Li": 1,
                "Co": 2,
                "Cl": -1,
            },
        )

        # add molecular weight
        self.molar_mass = Param(
            self.component_list,
            units=units.kg / units.mol,
            initialize={
                "Li": 0.006941,
                "Co": 0.05893,
                "Cl": 0.03545,
            },
        )

        # add thermal reflection coefficient, where 1 represents ideal behavior
        self.sigma = Param(
            self.component_list,
            units=units.dimensionless,
            initialize={
                "Li": 1,
                "Co": 1,
                "Cl": 1,
            },
        )

        self.num_solutes = Param(
            initialize=5,
            units=units.dimensionless,
            doc="Number of dissociated ions in solution",
        )

        self._state_block_class = SoluteStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "flow_vol": {"method": None},
                "conc_mass_comp": {"method": None},
                "flow_mol_comp": {"method": None},
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


class _SoluteStateBlock(StateBlock):
    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        fix_state_vars(self)


@declare_process_block_class("SoluteStateBlock", block_class=_SoluteStateBlock)
class SoluteStateBlockData(StateBlockData):
    """
    State block for multi-salt diafiltration membrane
    """

    def build(self):
        super().build()

        self.flow_vol = Var(
            units=units.m**3 / units.h,
            initialize=10,
            bounds=(1e-20, None),
        )
        self.conc_mass_comp = Var(
            self.component_list,
            units=units.kg / units.m**3,
            initialize=1e-5,
            bounds=(1e-20, None),
        )

    def get_material_flow_terms(self, p, j):
        return self.flow_vol * self.conc_mass_comp[j]

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mass

    def define_state_vars(self):
        return {"flow_vol": self.flow_vol, "conc_mass_comp": self.conc_mass_comp}
