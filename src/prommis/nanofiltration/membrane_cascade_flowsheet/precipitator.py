#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Precipitator unit model.

Modification of the IDAES Separator unit for yield based precipitation.
"""

from pyomo.common.config import ConfigValue

# Pyomo import
from pyomo.environ import Constraint, Param, Var, exp, units

from idaes.core import declare_process_block_class
from idaes.core.util.exceptions import ConfigurationError

# IDAES imports
from idaes.models.unit_models.separator import SeparatorData

__author__ = "Jason Yao"


@declare_process_block_class("Precipitator")
class SplitterData(SeparatorData):
    """Modification of Separator unit for precipitation."""

    CONFIG = SeparatorData.CONFIG()

    CONFIG.declare(
        "yields", ConfigValue(domain=dict, doc="Dictionary of precipitator yields")
    )

    def build(self):
        """
        Build Separator.

        Remove all split fractions and constraints, then add project
        specific constraints.
        """
        super().build()
        self.deactivate_all_cons()
        self._verify_precipitator_inputs()
        self.add_precipitator_constraints()

    def deactivate_all_cons(self):
        """Deactivate all current constraints."""
        for con in self.component_data_objects(Constraint):
            con.deactivate()

    def _verify_precipitator_inputs(self):
        """Check that precipitator inputs are correct."""
        # check outlets
        if not self.config.outlet_list:
            raise ConfigurationError(
                "Precipitator unit must be provided with `outlet_list` arg."
            )
        if len(self.config.outlet_list) != 2:
            raise ConfigurationError(
                "Precipitator unit must be provided with two outlets."
            )

        # check yields
        sol = [i for i in self.mixed_state.component_list if i != "solvent"]
        if self.index():
            check_yields = self.config.yields[self.index()]
        else:
            check_yields = self.config.yields

        if len(check_yields) != len(sol):
            raise ConfigurationError(
                "Precipitator yields must have same number of solutes as inlet flow."
                f" inlet solutes are: {sol}"
            )

        for i in check_yields:
            if i not in sol:
                raise ConfigurationError(
                    "Precipitator yields must have same solutes as inlet flow."
                    f" inlet solutes are: {sol}"
                )

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
        self.yields = Var(solutes, self.outlet_idx)

        # set solvent product outlet to be 0
        self.yields["solvent", "solid"].fix(0)

        # Sum of flows yields is 1
        @self.Constraint(solutes, doc="Sum of yields equation")
        def yields_eqn(b, sol):
            return sum(b.yields[sol, i] for i in self.outlet_idx) == 1

        @self.Constraint(
            self.flowsheet().time,
            self.outlet_idx,
            solutes,
            doc="Precipitator outlet equations",
        )
        def outlet_yield_eqn(b, t, o, sol):
            o_block = getattr(self, o + "_state")
            if sol == "solvent":
                return (
                    b.yields[sol, o] * b.mixed_state[t].flow_vol == o_block[t].flow_vol
                )
            return (
                b.yields[sol, o] * b.mixed_state[t].flow_mass_solute[sol]
                == o_block[t].flow_mass_solute[sol]
            )

        # precipitator volume and residence time relations
        self.tau = Var(units=units.hour)
        # initial point 100 m^3
        self.volume = Var(units=units.m**3, bounds=(0, 1000), initialize=100)
        self.volume.fixed = True
        sol = [i for i in self.mixed_state.component_list if i != "solvent"]

        # alpha controls maximum yield of precipitator
        # if precipitator is indexed, access appropriate yield values by index
        if self.index():
            self.alpha = Param(sol, initialize=self.config.yields[self.index()])
        else:
            self.alpha = Param(sol, initialize=self.config.yields)
        self.beta = Param(sol, initialize={s: 4.6 for s in sol}, units=1 / units.hour)

        @self.Constraint(sol)
        def prec_res_time(b, sol):
            return b.yields[sol, "solid"] == b.alpha[sol] * (
                1 - exp(-b.beta[sol] * b.tau)
            )

        @self.Constraint()
        def prec_vol(b):
            return b.volume == b.mixed_state[0].flow_vol * b.tau
