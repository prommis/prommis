"""
Precipitator unit model.

Modification of the IDAES Separator unit for yield based precipitation.
"""

# Pyomo import
import pyomo.environ as pyo

# IDAES imports
from idaes.models.unit_models.separator import SeparatorData
from idaes.core import declare_process_block_class
from pyomo.common.config import ConfigValue


__author__ = "Jason Yao"


@declare_process_block_class("Precipitator")
class SplitterData(SeparatorData):
    """Modification of Separator unit for precipitation."""

    CONFIG = SeparatorData.CONFIG()

    CONFIG.declare(
        "yields",
        ConfigValue(
            domain=dict,
            doc="Dictionary of precipitator yields"
        )
    )

    def build(self):
        """
        Build Separator.

        Remove all split fractions and constraints, then add project
        specific constraints.
        """
        super().build()
        # TODO add input checking to prevent incorrect setup
        self.deactivate_all_cons()
        # self.sum_split_frac[0].activate()    # activate sum sf == 1 constraint
        self.add_precipitator_constraints()

    def deactivate_all_cons(self):
        """Deactivate all current constraints."""
        for con in self.component_data_objects(pyo.Constraint):
            # print(con)
            con.deactivate()

    def add_precipitator_constraints(self):
        """
        Add Precipitator constraints.

        1. recycle flow contains all water and dilute solutes
        2. split fraction
        """
        # precipitator yield variables
        # solutes = [i for i in self.mixed_state.component_list
        #            if i != 'solvent']
        # or include generalization for some water in precipitate
        solutes = self.mixed_state.component_list
        self.yields = pyo.Var(solutes, self.outlet_idx)

        #####
        # Add bypass
        #####
        bypass_flows = ['bypass', 'ro']
        self.split_inlet = pyo.Var(bypass_flows)
        self.split_inlet_flows = pyo.Var(solutes, bypass_flows)

        # set initial bypass to 0
        self.split_inlet['bypass'].fix(1e-8)
        @self.Constraint(
            doc="Sum of bypass split frac equation"
        )
        def split_inlet_eqn(b):
            return (
                sum(b.split_inlet[i] for i in bypass_flows)
                == 1
            )

        @self.Constraint(
            self.flowsheet().time,
            solutes,
            bypass_flows,
            doc="Precipitator outlet equations"
        )
        def bypass_split_eqn(b, t, sol, flows):
            if sol == 'solvent':
                return (
                    b.split_inlet[flows]*b.mixed_state[t].flow_vol
                    == b.split_inlet_flows['solvent', flows]
                )
            return (
                b.split_inlet[flows]*b.mixed_state[t].mass_solute[sol]
                == b.split_inlet_flows[sol, flows]
            )

        #####
        # Reverse osmosis
        #####
        # set 50% of solvent to go to recycle
        self.yields['solvent', 'recycle'].fix(0.5)

        # set solute recycle outlet to be 0
        self.yields['Li', 'recycle'].fix(1e-8)
        self.yields['Co', 'recycle'].fix(1e-8)

        #####
        # Product collection
        #####
        # set solvent product outlet to be 0
        self.yields['solvent', 'solid'].fix(1e-8)

        # Sum of flows yields is 1
        @self.Constraint(
            solutes,
            doc="Sum of yields equation"
        )
        def yields_eqn(b, sol):
            return (
                sum(b.yields[sol, i] for i in self.outlet_idx)
                == 1
            )

        @self.Constraint(
            self.flowsheet().time,
            self.outlet_idx,
            solutes,
            doc="Precipitator outlet equations"
        )
        def outlet_yield_eqn(b, t, o, sol):
            o_block = getattr(self, o + "_state")
            if sol == 'solvent':
                if o == 'recycle':
                    return (
                        b.yields[sol, o]*b.split_inlet_flows['solvent', 'ro']
                        + b.split_inlet_flows['solvent', 'bypass']
                        == o_block[t].flow_vol
                    )
                else:
                    return (
                        b.yields[sol, o]*b.split_inlet_flows['solvent', 'ro']
                        == o_block[t].flow_vol
                    )
            if o == 'recycle':
                return (
                    b.yields[sol, o]*b.split_inlet_flows[sol, 'ro']
                    + b.split_inlet_flows[sol, 'bypass']
                    == o_block[t].mass_solute[sol]
                )
            else:
                return (
                    b.yields[sol, o]*b.split_inlet_flows[sol, 'ro']
                    == o_block[t].mass_solute[sol]
                )

        # precipitator volume and residence time relations
        self.tau = pyo.Var(units=pyo.units.hour)
        self.V = pyo.Var(units=pyo.units.m**3)
        self.V.fix(100)  # initial point m^3
        self.V.setub(1000)  # set UB to prevent ridiculous sizes
        sol = [i for i in self.mixed_state.component_list
               if i != 'solvent']
        self.alpha = pyo.Param(sol,
                               initialize=self.config.yields[self.index()])
        self.beta = pyo.Param(sol, initialize={s: 5 for s in sol}, units=1/pyo.units.hour)

        @self.Constraint(sol)
        def prec_res_time(b, sol):
            return (
                b.yields[sol, 'solid']
                == b.alpha[sol]
                * (1 - pyo.exp(-b.beta[sol]*b.tau))
             )

        @self.Constraint()
        def prec_vol(b):
            return (
                b.V == (1 - b.yields['solvent', 'recycle'])*b.split_inlet_flows['solvent', 'ro']*b.tau
            )

        # also adding concentration upper bounds
        @ self.Constraint([sol for sol in solutes if 'solvent' not in sol])
        def conc_limits(b, sol):
            if sol == 'Li':
                conc_ub = 20 * pyo.units.kg / pyo.units.m**3    # kg/m^3
            elif sol == 'Co':
                conc_ub = 200 * pyo.units.kg / pyo.units.m**3   # kg/m^3
            else:
                conc_ub = 200 *pyo.units.kg / pyo.units.m**3   # generic UB set to 200 kg/m^3
            return (
                conc_ub*(1 - b.yields['solvent', 'recycle'])*b.split_inlet_flows['solvent', 'ro']
                >= b.split_inlet_flows[sol, 'ro']
            )
