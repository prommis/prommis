#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
Translator block converting from precipitator properties to solvent extraction properties.
This is copied from the IDAES Generic template for a translator block.

Assumptions:
     * Steady-state only
"""

# Import Pyomo libraries
from pyomo.environ import (
    units as pyunits,
    check_optimal_termination,
    Set,
)

# Import IDAES cores
from idaes.core import declare_process_block_class
from idaes.models.unit_models.translator import TranslatorData
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import InitializationError

import idaes.logger as idaeslog

__author__ = "Marcus Holly"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("Translator_precipitator_SX")
class TranslatorDataPrecipitatorSX(TranslatorData):
    """
    Translator block representing the precipitator/SX interface
    """

    def build(self):
        """
        Begin building model.
        Args:
            None
        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super().build()

        mw_al = 26.98154e3 * pyunits.mg / pyunits.mol
        mw_ca = 39.96259e3 * pyunits.mg / pyunits.mol
        mw_fe = 55.93494e3 * pyunits.mg / pyunits.mol
        mw_ce = 140.116e3 * pyunits.mg / pyunits.mol

        DEHPA_density = 0.9758 * pyunits.kg/pyunits.L

        @self.Expression(
            self.flowsheet().time,
            doc="Volumetric flow (L/hr)",
        )
        def volumetric_flow(blk, t):
            return (
                    pyunits.convert(blk.properties_in[t].flow_mass, to_units=pyunits.kg / pyunits.hour)
                    / DEHPA_density
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality volumetric flow equation (L/hr)",
        )
        def eq_flow_vol_rule(blk, t):
            return (
                    blk.properties_out[t].flow_vol
                    == blk.volumetric_flow[t]
            )

        @self.Expression(
            self.flowsheet().time,
            doc="Aluminum molar flow (mol/s)",
        )
        def aluminum_molar_flow(blk, t):
            return (
                10 ** blk.properties_in[t].log10_molality_comp["Al^3+"] * pyunits.mol / pyunits.kg
                * blk.properties_in[t].flow_mass
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Aluminum concentration (mg/L)",
        )
        def aluminum_conc_mass_comp(blk, t):
            return (
                blk.properties_out[t].conc_mass_comp["Al"]
                == pyunits.convert(blk.aluminum_molar_flow[t], to_units=pyunits.mol / pyunits.hour,)
                * mw_al / blk.volumetric_flow[t]
            )

        @self.Expression(
            self.flowsheet().time,
            doc="Calcium molar flow (mol/s)",
        )
        def calcium_molar_flow(blk, t):
            return (
                10 ** blk.properties_in[t].log10_molality_comp["Ca^2+"] * pyunits.mol / pyunits.kg
                * blk.properties_in[t].flow_mass
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Calcium concentration (mg/L)",
        )
        def calcium_conc_mass_comp(blk, t):
            return (
                blk.properties_out[t].flow_mass["Ca"]
                == pyunits.convert(blk.calcium_molar_flow[t], to_units=pyunits.mol / pyunits.hour,)
                * mw_ca / blk.volumetric_flow[t]
            )

        @self.Expression(
            self.flowsheet().time,
            doc="Iron molar flow (mol/s)",
        )
        def iron_molar_flow(blk, t):
            return (
                10 ** blk.properties_in[t].log10_molality_comp["Fe^3+"] * pyunits.mol / pyunits.kg
                * blk.properties_in[t].flow_mass
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Iron concentration (mg/L)",
        )
        def iron_conc_mass_comp(blk, t):
            return (
                blk.properties_out[t].flow_mass["Fe"]
                == pyunits.convert(blk.iron_molar_flow[t], to_units=pyunits.mol / pyunits.hour,)
                * mw_fe / blk.volumetric_flow[t]
            )

        @self.Expression(
            self.flowsheet().time,
            doc="Cerium molar flow (mol/s)",
        )
        def cerium_molar_flow(blk, t):
            return (
                10 ** blk.properties_in[t].log10_molality_comp["Ce^3+"] * pyunits.mol / pyunits.kg
                * blk.properties_in[t].flow_mass
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Cerium concentration (mg/L)",
        )
        def cerium_conc_mass_comp(blk, t):
            return (
                blk.properties_out[t].flow_mass["Ce"]
                == pyunits.convert(blk.iron_molar_flow[t], to_units=pyunits.mol / pyunits.hour,)
                * mw_ce / blk.volumetric_flow[t]
            )

        # TODO: Should these components be removed from the SX model?
        self.components = Set(initialize=["Y", "La", "Pr", "Nd", "Sm", "Gd", "Dy", "Sc"])

        @self.Constraint(
            self.flowsheet().time,
            self.components,
            doc="Component concentration (mg/L)",
        )
        def zero_conc_components(blk, t, i):
            return (
                blk.properties_out[t].conc_mass_comp[i]
                == 1e-7 * pyunits.g / pyunits.hour
            )

    def initialize_build(
        self,
        state_args_in=None,
        state_args_out=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        This method calls the initialization method of the state blocks.

        Keyword Arguments:
            state_args_in : a dict of arguments to be passed to the inlet
                property package (to provide an initial state for
                initialization (see documentation of the specific
                property package) (default = None).
            state_args_out : a dict of arguments to be passed to the outlet
                property package (to provide an initial state for
                initialization (see documentation of the specific
                property package) (default = None).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)

        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialize state block
        flags = self.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_in,
            hold_state=True,
        )

        self.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )

        if degrees_of_freedom(self) != 0:
            raise Exception(
                f"{self.name} degrees of freedom were not 0 at the beginning "
                f"of initialization. DoF = {degrees_of_freedom(self)}"
            )

        with idaeslog.solver_log(init_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)

        self.properties_in.release_state(flags=flags, outlvl=outlvl)

        init_log.info(f"Initialization Complete: {idaeslog.condition(res)}")

        if not check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )
