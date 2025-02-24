from pyomo.environ import (
    Constraint,
    Var,
    Param,
    Set,
    Expression,
    Reals,
    NonNegativeReals,
    Suffix,
)
from pyomo.environ import units as pyunits
from pyomo.common.config import ConfigValue

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    Phase,
    MaterialFlowBasis,
    PhysicalParameterBlock,
    StateBlockData,
    StateBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from idaes.core.base.components import Component, Solute, Solvent
from idaes.core.base.phases import LiquidPhase, SolidPhase, PhaseType as PT
from idaes.core.util.initialization import (
    fix_state_vars,
    revert_state_vars,
    solve_indexed_blocks,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_unfixed_variables,
)
from idaes.core.util.exceptions import PropertyPackageError
from idaes.core.util.misc import extract_data
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.solvers import get_solver

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("AqueousParameter")
class AqueousParameterData(PhysicalParameterBlock):
    CONFIG = PhysicalParameterBlock.CONFIG()
    CONFIG.declare(
        "aq_comp_list", ConfigValue(domain=list, description="List of aqueous components in process")
    )
    CONFIG.declare(
        "eq_rxn_logkeq_dict", ConfigValue(domain=dict, description="Dictionary of equilibrium reaction constants")
    )
    CONFIG.declare(
        "eq_rxn_stoich_dict", ConfigValue(domain=dict, description="Dictionary of equilibrium reaction stoichiometry")
    )

    def build(self):
        """
        Calllable method for Block construction.
        """
        # super(AqueousParameterData, self).build()
        super().build()

        self.AqueousPhase = Phase()
        self.component_list = self.config.aq_comp_list

        ## equilibrium reaction parameters
        # equilibrium reaction index
        self.eq_rxn_set = Set(initialize=list(set(key for key in self.config.eq_rxn_logkeq_dict.keys())))
        
        # stoichiometry for each equilibrium reaction
        self.eq_rxn_stoich_dict = self.config.eq_rxn_stoich_dict

        # log(keq) for each equilibrium reaction
        self.eq_rxn_logkeq_dict = self.config.eq_rxn_logkeq_dict

        self._state_block_class = AqueousStateBlock
        

        @classmethod
        def define_metadata(cls, obj):
            """Define properties supported and units."""
            obj.add_properties(
                {
                    "molality_aq_comp": {"method": None},
                }
            )
            obj.add_default_units(
                {
                "time": pyunits.hour,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
                }
            )


class _AqueousStateBlock(StateBlock):
    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        # Fix state variables
        return fix_state_vars(self)
    

@declare_process_block_class("AqueousStateBlock", block_class=_AqueousStateBlock)
class AqueousStateBlockData(StateBlockData):
    def build(self):
        super().build()

        self.molality_aq_comp = Var(
            self.params.component_list,
            units=pyunits.mol / pyunits.kg,
            initialize=1e-20,
            bounds=(1e-20, None),
        )

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        return {
            "molality_aq_comp": self.molality_aq_comp,
        }

