#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
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
        Na (sodium ion, +)
        Ca (calcium ion, 2+)
        Co (cobalt ion, 2+)
        Al (aluminum ion, 3+)
        La (lanthanum ion, 3+)
        Cl (chloride ion, -)
    in the following salt mixtures:
        all single chloride salts (e.g., LiCl)
        LiCl + CoCl2
        LiCl + AlCl3
        CoCl2 + AlCl3
        LiCl + CoCl2 + AlCl3
        NaCl + CaCl2
        NaCl + LaCl3
        CaCl2 + LaCl3
        NaCl + CaCl2 + LaCl3
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
    CONFIG.declare(
        "non_Donnan_partition_dict",
        ConfigValue(
            default={"Li": 1, "Co": 1, "Cl": 1},
            doc="Dictionary of non-Donnan partition coefficients",
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
            "Na": 1,
            "Ca": 2,
            "Co": 2,
            "Al": 3,
            "La": 3,
            "Cl": -1,
        }

        # infinite dilution solute diffusion coefficient
        # source: https://www.aqion.de/site/diffusion-coefficients
        boundary_layer_diffusion_coefficient_dict = {
            "Li": 3.71,  # mm2 / h
            "Na": 4.79,  # mm2 / h
            "Ca": 2.85,  # mm2 / h
            "Co": 2.64,  # mm2 / h
            "Al": 2.01,  # mm2 / h
            "La": 2.23,  # mm2 / h
            "Cl": 7.31,  # mm2 / h
        }

        # membrane phase diffusion coefficient
        # assumed 3 orders of magnitude smaller than solution phase
        membrane_diffusion_coefficient_dict = {
            "Li": 3.71 * 0.001,  # mm2 / h
            "Na": 4.79 * 0.001,  # mm2 / h
            "Ca": 2.85 * 0.001,  # mm2 / h
            "Co": 2.64 * 0.001,  # mm2 / h
            "Al": 2.01 * 0.001,  # mm2 / h
            "La": 2.23 * 0.001,  # mm2 / h
            "Cl": 7.31 * 0.001,  # mm2 / h
        }

        # thermal reflection coefficient, related to solute rejection
        # assumed 1
        sigma_dict = {
            "Li": 1,  # mm2 / h
            "Na": 1,  # mm2 / h
            "Ca": 1,  # mm2 / h
            "Co": 1,  # mm2 / h
            "Al": 1,  # mm2 / h
            "La": 1,  # mm2 / h
            "Cl": 1,  # mm2 / h
        }

        if self.config.cation_list == ["Li"]:
            salt_system = "Li_Cl"
        elif self.config.cation_list == ["Na"]:
            salt_system = "Na_Cl"
        elif self.config.cation_list == ["Ca"]:
            salt_system = "Ca_Cl2"
        elif self.config.cation_list == ["Co"]:
            salt_system = "Co_Cl2"
        elif self.config.cation_list == ["Al"]:
            salt_system = "Al_Cl3"
        elif self.config.cation_list == ["La"]:
            salt_system = "La_Cl3"
        elif self.config.cation_list == ["Li", "Co"]:
            salt_system = "Li_Co_Cl3"
        elif self.config.cation_list == ["Li", "Al"]:
            salt_system = "Li_Al_Cl4"
        elif self.config.cation_list == ["Co", "Al"]:
            salt_system = "Co_Al_Cl5"
        elif self.config.cation_list == ["Li", "Co", "Al"]:
            salt_system = "Li_Co_Al_Cl6"
        elif self.config.cation_list == ["Na", "Ca"]:
            salt_system = "Na_Ca_Cl3"
        elif self.config.cation_list == ["Na", "La"]:
            salt_system = "Na_La_Cl4"
        elif self.config.cation_list == ["Ca", "La"]:
            salt_system = "Ca_La_Cl5"
        elif self.config.cation_list == ["Na", "Ca", "La"]:
            salt_system = "Na_Ca_La_Cl6"

        num_solutes_dict = {
            "Li_Cl": {"Li": 1, "Cl": 1},
            "Na_Cl": {"Na": 1, "Cl": 1},
            "Ca_Cl2": {"Ca": 1, "Cl": 2},
            "Co_Cl2": {"Co": 1, "Cl": 2},
            "Al_Cl3": {"Al": 1, "Cl": 3},
            "La_Cl3": {"La": 1, "Cl": 3},
            "Li_Co_Cl3": {"Li": 1, "Co": 1, "Cl": 3},
            "Li_Al_Cl4": {"Li": 1, "Al": 1, "Cl": 4},
            "Co_Al_Cl5": {"Co": 1, "Al": 1, "Cl": 5},
            "Li_Co_Al_Cl6": {"Li": 1, "Co": 1, "Al": 1, "Cl": 6},
            "Na_Ca_Cl3": {"Na": 1, "Ca": 1, "Cl": 3},
            "Na_La_Cl4": {"Na": 1, "La": 1, "Cl": 4},
            "Ca_La_Cl5": {"Ca": 1, "La": 1, "Cl": 5},
            "Na_Ca_La_Cl6": {"Na": 1, "Ca": 1, "La": 1, "Cl": 6},
        }

        # create subset of property dictionaries to initialize parameters
        def _subset(mapping_dict):
            return {ion: mapping_dict[ion] for ion in self.component_list}

        initialize_charge_dict = _subset(charge_dict)
        initialize_boundary_layer_diffusion_coefficient_dict = _subset(
            boundary_layer_diffusion_coefficient_dict
        )
        initialize_membrane_diffusion_coefficient_dict = _subset(
            membrane_diffusion_coefficient_dict
        )
        initialize_sigma_dict = _subset(sigma_dict)
        initialize_num_solutes_dict = _subset(num_solutes_dict[salt_system])

        # initialize properties
        self.charge = Param(
            self.component_list,
            units=units.dimensionless,
            initialize=initialize_charge_dict,
        )

        self.boundary_layer_diffusion_coefficient = Param(
            self.component_list,
            units=units.mm**2 / units.h,
            initialize=initialize_boundary_layer_diffusion_coefficient_dict,
        )

        self.membrane_diffusion_coefficient = Param(
            self.component_list,
            units=units.mm**2 / units.h,
            initialize=initialize_membrane_diffusion_coefficient_dict,
        )

        self.sigma = Param(
            self.component_list,
            units=units.dimensionless,
            initialize=initialize_sigma_dict,
        )

        self.non_Donnan_partition_coefficient = Param(
            self.component_list,
            units=units.dimensionless,
            initialize=self.config.non_Donnan_partition_dict,
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
