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

from pyomo.common.config import ConfigValue, ListOf
from pyomo.environ import Param, Set, Var, units

from idaes.core import (
    MaterialFlowBasis,
    Phase,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    declare_process_block_class,
)
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.initialization import fix_state_vars


@declare_process_block_class("MultiComponentDiafiltrationSoluteParameter")
class MultiComponentDiafiltrationSoluteParameterData(PhysicalParameterBlock):
    """
    Property Package for the multi-component diafiltration membrane.

    Currently includes the following solutes:
        Li (lithium ion, +)
        Co (cobalt ion, 2+)
        Al (aluminum ion, 3+)
        Cl (chloride ion, -)
    """

    CONFIG = PhysicalParameterBlock.CONFIG()

    CONFIG.declare(
        "cation_list",
        ConfigValue(
            domain=ListOf(str),
            default=["Li", "Co"],
            doc="List of cations present in the system",
        ),
    )
    CONFIG.declare(
        "anion_list",
        ConfigValue(
            domain=ListOf(str),
            default=["Cl"],
            doc="List of anions present in the system",
        ),
    )

    def build(self):
        super().build()

        if len(self.config.anion_list) > 1:
            raise ConfigurationError(
                "The multi-component diafiltration unit model only supports systems with a common anion"
            )

        self.liquid = Phase()

        self.component_list = Set(
            initialize=self.config.cation_list + self.config.anion_list
        )

        # ion valence
        charge_dict = {
            "Li": 1,
            "Co": 2,
            "Al": 3,
            "Cl": -1,
        }

        # infinite dilution solute diffusion coefficient
        # source: https://www.aqion.de/site/diffusion-coefficients
        diffusion_coefficient_dict = {
            "Li": 3.71,  # mm2 / h
            "Co": 2.64,  # mm2 / h
            "Al": 2.01,  # mm2 / h
            "Cl": 7.31,  # mm2 / h
        }

        # thermal reflection coefficient, related to solute rejection
        sigma_dict = {
            "Li": 1,
            "Co": 1,
            "Al": 1,
            "Cl": 1,
        }

        # partition coefficient at the solution-membrane interfaces
        # Reference: https://doi.org/10.1126/sciadv.adu8302
        # Assumptions:
        # membrane fixed charge is negative (Donnan effects are incorporated)
        # monovalent ions of similar size (i.e., Na and Li) behave similarly
        # H,Li is estimated from the data in Fig 1D (Na) of above reference at 200 mM
        # H,Co (divalent) is estimated as one order of magnitude smaller than H,Li (monovalent)
        # H,Al (trivalent) is estimated as one order of magnitude smaller than H,Co (divalent)
        # H,Cl is estimated from the data in Fig 1C of above reference at 200 mM
        # while H on the retentate and permeate sides can differ, we assume them to be equal for now
        partition_coefficient_dict = {
            "retentate": {
                "Li": 0.4,
                "Co": 0.04,
                "Al": 0.004,
                "Cl": 0.01,
            },
            "permeate": {
                "Li": 0.4,
                "Co": 0.04,
                "Al": 0.004,
                "Cl": 0.01,
            },
        }

        if self.config.cation_list == ["Li"]:
            salt_system = "Li_Cl"
        elif self.config.cation_list == ["Co"]:
            salt_system = "Co_Cl2"
        elif self.config.cation_list == ["Al"]:
            salt_system = "Al_Cl3"
        elif self.config.cation_list == ["Li", "Co"]:
            salt_system = "Li_Co_Cl3"
        elif self.config.cation_list == ["Li", "Al"]:
            salt_system = "Li_Al_Cl4"
        elif self.config.cation_list == ["Co", "Al"]:
            salt_system = "Co_Al_Cl5"
        elif self.config.cation_list == ["Li", "Co", "Al"]:
            salt_system = "Li_Co_Al_Cl6"

        num_solutes_dict = {
            "Li_Cl": {
                "Li": 1,
                "Cl": 1,
            },
            "Co_Cl2": {
                "Co": 1,
                "Cl": 2,
            },
            "Al_Cl3": {
                "Al": 1,
                "Cl": 3,
            },
            "Li_Co_Cl3": {
                "Li": 1,
                "Co": 1,
                "Cl": 3,
            },
            "Li_Al_Cl4": {
                "Li": 1,
                "Al": 1,
                "Cl": 4,
            },
            "Co_Al_Cl5": {
                "Co": 1,
                "Al": 1,
                "Cl": 5,
            },
            "Li_Co_Al_Cl6": {
                "Li": 1,
                "Co": 1,
                "Al": 1,
                "Cl": 6,
            },
        }

        # initialize dictionaries for a single cation
        initialize_charge_dict = {
            self.config.cation_list[0]: charge_dict[self.config.cation_list[0]],
            self.config.anion_list[0]: charge_dict[self.config.anion_list[0]],
        }
        initialize_diffusion_coefficient_dict = {
            self.config.cation_list[0]: diffusion_coefficient_dict[
                self.config.cation_list[0]
            ],
            self.config.anion_list[0]: diffusion_coefficient_dict[
                self.config.anion_list[0]
            ],
        }
        initialize_sigma_dict = {
            self.config.cation_list[0]: sigma_dict[self.config.cation_list[0]],
            self.config.anion_list[0]: sigma_dict[self.config.anion_list[0]],
        }
        initialize_partition_coefficient_retentate_dict = {
            self.config.cation_list[0]: partition_coefficient_dict["retentate"][
                self.config.cation_list[0]
            ],
            self.config.anion_list[0]: partition_coefficient_dict["retentate"][
                self.config.anion_list[0]
            ],
        }
        initialize_partition_coefficient_permeate_dict = {
            self.config.cation_list[0]: partition_coefficient_dict["permeate"][
                self.config.cation_list[0]
            ],
            self.config.anion_list[0]: partition_coefficient_dict["permeate"][
                self.config.anion_list[0]
            ],
        }
        initialize_num_solutes_dict = {
            self.config.cation_list[0]: num_solutes_dict[salt_system][
                self.config.cation_list[0]
            ],
            self.config.anion_list[0]: num_solutes_dict[salt_system][
                self.config.anion_list[0]
            ],
        }

        # add additional cations to dictionaries
        i = 1
        while i < len(self.config.cation_list):
            initialize_charge_dict.update(
                {self.config.cation_list[i]: charge_dict[self.config.cation_list[i]]}
            )
            initialize_diffusion_coefficient_dict.update(
                {
                    self.config.cation_list[i]: diffusion_coefficient_dict[
                        self.config.cation_list[i]
                    ]
                }
            )
            initialize_sigma_dict.update(
                {self.config.cation_list[i]: sigma_dict[self.config.cation_list[i]]}
            )
            initialize_partition_coefficient_retentate_dict.update(
                {
                    self.config.cation_list[i]: partition_coefficient_dict["retentate"][
                        self.config.cation_list[i]
                    ]
                }
            )
            initialize_partition_coefficient_permeate_dict.update(
                {
                    self.config.cation_list[i]: partition_coefficient_dict["permeate"][
                        self.config.cation_list[i]
                    ]
                }
            )
            initialize_num_solutes_dict.update(
                {
                    self.config.cation_list[i]: num_solutes_dict[salt_system][
                        self.config.cation_list[i]
                    ]
                }
            )
            i += 1

        # initialize properties
        self.charge = Param(
            self.component_list,
            units=units.dimensionless,
            initialize=initialize_charge_dict,
        )

        self.diffusion_coefficient = Param(
            self.component_list,
            units=units.mm**2 / units.h,
            initialize=initialize_diffusion_coefficient_dict,
        )

        self.sigma = Param(
            self.component_list,
            units=units.dimensionless,
            initialize=initialize_sigma_dict,
        )

        self.partition_coefficient_retentate = Param(
            self.component_list,
            units=units.dimensionless,
            initialize=initialize_partition_coefficient_retentate_dict,
        )

        self.partition_coefficient_permeate = Param(
            self.component_list,
            units=units.dimensionless,
            initialize=initialize_partition_coefficient_permeate_dict,
        )

        self.num_solutes = Param(
            self.component_list,
            units=units.dimensionless,
            initialize=initialize_num_solutes_dict,
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
