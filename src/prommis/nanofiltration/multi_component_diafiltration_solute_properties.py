#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Property package for the multi-component diafiltration membrane.

Author: Molly Dougher
"""

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
from idaes.core.util.initialization import fix_state_vars


@declare_process_block_class("MultiComponentDiafiltrationSoluteParameter")
class MultiComponentDiafiltrationSoluteParameterData(PhysicalParameterBlock):
    """
    Property Package for the multi-component diafiltration membrane.

    Currently includes the following solutes:
        Li+ (lithium ion)
        Co2+ (cobalt ion)
        Cl- (chloride ion)
    """

    CONFIG = PhysicalParameterBlock.CONFIG()

    CONFIG.declare(
        "salt_system",
        ConfigValue(
            default="lithium_cobalt_chloride",
            doc="Name of the salt system to be modeled",
        ),
    )

    def build(self):
        super().build()

        self.liquid = Phase()

        # add cations
        self.cation_1 = Component()
        self.cation_2 = Component()

        # add anions
        self.anion = Component()

        # ion valence
        charge_dict = {
            "Li": 1,
            "Co": 2,
            "Cl": -1,
        }

        # infinite dilution solute diffusion coefficient
        # source: https://www.aqion.de/site/diffusion-coefficients
        diffusion_coefficient_dict = {
            "Li": 3.71,  # mm2 / h
            "Co": 2.64,  # mm2 / h
            "Cl": 7.31,  # mm2 / h
        }

        # thermal reflection coefficient, related to solute rejection
        sigma_dict = {
            "Li": 1,
            "Co": 1,
            "Cl": 1,
        }

        # partition coefficient at the solution-membrane interfaces
        # Reference: https://doi.org/10.1126/sciadv.adu8302
        # Assumptions:
        # membrane fixed charge is negative (Donnan effects are incorporated)
        # monovalent ions of similar size (i.e., Na and Li) behave similarly
        # H,Li is estimated from the data in Fig 1D (Na) of above reference at 200 mM
        # H,Co (divalent) is estimated as one order of magnitude smaller than H,Li (monovalent)
        # H,Cl is estimated from the data in Fig 1C of above reference at 200 mM
        # while H on the retentate and permeate sides can differ, we assume them to be equal for now
        partition_coefficient_dict = {
            "retentate": {
                "Li": 0.4,
                "Co": 0.04,
                "Cl": 0.01,
            },
            "permeate": {
                "Li": 0.4,
                "Co": 0.04,
                "Cl": 0.01,
            },
        }

        num_solutes_dict = {
            "lithium_cobalt_chloride": {
                "Li": 1,
                "Co": 1,
                "Cl": 3,
            }
        }

        if self.config.salt_system == "lithium_cobalt_chloride":
            cation_1 = "Li"
            cation_2 = "Co"
            anion = "Cl"

        self.charge = Param(
            self.component_list,
            units=units.dimensionless,
            initialize={
                "cation_1": charge_dict[cation_1],
                "cation_2": charge_dict[cation_2],
                "anion": charge_dict[anion],
            },
        )

        self.diffusion_coefficient = Param(
            self.component_list,
            units=units.mm**2 / units.h,
            initialize={
                "cation_1": diffusion_coefficient_dict[cation_1],
                "cation_2": diffusion_coefficient_dict[cation_2],
                "anion": diffusion_coefficient_dict[anion],
            },
        )

        self.sigma = Param(
            self.component_list,
            units=units.dimensionless,
            initialize={
                "cation_1": sigma_dict[cation_1],
                "cation_2": sigma_dict[cation_2],
                "anion": sigma_dict[anion],
            },
        )

        self.partition_coefficient_retentate = Param(
            self.component_list,
            units=units.dimensionless,
            initialize={
                "cation_1": partition_coefficient_dict["retentate"][cation_1],
                "cation_2": partition_coefficient_dict["retentate"][cation_2],
                "anion": partition_coefficient_dict["retentate"][anion],
            },
        )

        self.partition_coefficient_permeate = Param(
            self.component_list,
            units=units.dimensionless,
            initialize={
                "cation_1": partition_coefficient_dict["permeate"][cation_1],
                "cation_2": partition_coefficient_dict["permeate"][cation_2],
                "anion": partition_coefficient_dict["permeate"][anion],
            },
        )

        self.num_solutes = Param(
            self.component_list,
            units=units.dimensionless,
            initialize={
                "cation_1": num_solutes_dict[self.config.salt_system][cation_1],
                "cation_2": num_solutes_dict[self.config.salt_system][cation_2],
                "anion": num_solutes_dict[self.config.salt_system][anion],
            },
            doc="Moles of ions dissociated in solution per mole of salt(s)",
        )

        self._state_block_class = MultiComponentDiafiltrationSoluteStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "flow_vol": {"method": None},
                "conc_mol_comp": {"method": None},
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


class _MultiComponentDiafiltrationSoluteStateBlock(StateBlock):
    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        fix_state_vars(self)


@declare_process_block_class(
    "MultiComponentDiafiltrationSoluteStateBlock",
    block_class=_MultiComponentDiafiltrationSoluteStateBlock,
)
class MultiComponentDiafiltrationSoluteStateBlockData(StateBlockData):
    """
    State block for multi-component diafiltration membrane
    """

    def build(self):
        super().build()

        self.flow_vol = Var(
            units=units.m**3 / units.h,
            initialize=10,
            bounds=(1e-20, None),
        )
        self.conc_mol_comp = Var(
            self.component_list,
            units=units.mol / units.m**3,
            initialize=1e-5,
            bounds=(1e-20, None),
        )

    def get_material_flow_terms(self, p, j):
        return self.flow_vol * self.conc_mol_comp[j]

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mole

    def define_state_vars(self):
        return {"flow_vol": self.flow_vol, "conc_mol_comp": self.conc_mol_comp}
