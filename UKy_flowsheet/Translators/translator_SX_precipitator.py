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
Translator block converting from solvent extraction properties to precipitator properties.
This is copied from the IDAES Generic template for a translator block.

Assumptions:
     * Steady-state only
"""

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue

# Import IDAES cores
from idaes.core import declare_process_block_class
from idaes.models.unit_models.translator import TranslatorData
from idaes.core.util.config import (
    is_reaction_parameter_block,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale

from idaes.core.util.exceptions import InitializationError

from pyomo.environ import (
    units as pyunits,
    check_optimal_termination,
    Set,
    value,
    log10 as log_10,
)

from math import log10

__author__ = "Marcus Holly"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("Translator_SX_precipitator")
class TranslatorDataSXPrecipitator(TranslatorData):
    """
    Translator block representing the SX/precipitator interface
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
        super(TranslatorDataSXPrecipitator, self).build()

        mw_al = 0.02698154 * pyunits.kg / pyunits.mol
        mw_ca = 0.03996259 * pyunits.kg / pyunits.mol
        mw_fe = 0.05593494 * pyunits.kg / pyunits.mol
        mw_ce = 0.140116 * pyunits.kg / pyunits.mol

        #TODO: Replace temperature and pressure constraints with an equality when SX has temp and pressure
        @self.Constraint(
            self.flowsheet().time,
            doc="Equality temperature equation",
        )
        def eq_temperature_rule(blk, t):
            return blk.properties_out[t].temperature == 300 * pyunits.kelvin

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality pressure equation",
        )
        def eq_pressure_rule(blk, t):
            return blk.properties_out[t].pressure == 101325 * pyunits.Pa

        @self.Expression(
            self.flowsheet().time,
            doc="Mass flow of solvent (kg/s)",
        )
        def solvent_mass_flow(blk, t):
            return (
                pyunits.convert(blk.properties_in[t].flow_vol, to_units=pyunits.L / pyunits.s,)
                * (1.840 * pyunits.kg/pyunits.L)
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality mass flow equation (kg/s)",
        )
        def eq_flow_mass_rule(blk, t):
            return (
                    blk.properties_out[t].flow_mass
                    == blk.solvent_mass_flow[t]
            )

        #TODO: may need to add oxalate, which will be added from a separate stream
        #TODO: may need a separate translator block for this or a mixer prior to this block?

        @self.Expression(
            self.flowsheet().time,
            doc="Aluminum molar flow (mol/s)",
        )
        def aluminum_molar_flow(blk, t):
            return (
                pyunits.convert(blk.properties_in[t].flow_mass["Al"], to_units=pyunits.kg / pyunits.s,)
                / mw_al
            )

        @self.Expression(
            self.flowsheet().time,
            doc="Aluminum dimensionless molality",
        )
        def aluminum_molality(blk, t):
            return (
                blk.aluminum_molar_flow[t] * pyunits.s / pyunits.mol
                / blk.solvent_mass_flow[t] * pyunits.kg / pyunits.s
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Aluminum log10 molality",
        )
        def eq_aluminum_log10molality(blk, t):
            return (
                blk.properties_out[t].log10_molality_comp["Al^3+"]
                == log_10(blk.aluminum_molality[t])
            )

        @self.Expression(
            self.flowsheet().time,
            doc="Calcium molar flow (mol/s)",
        )
        def calcium_molar_flow(blk, t):
            return (
                pyunits.convert(blk.properties_in[t].flow_mass["Ca"], to_units=pyunits.kg / pyunits.s,)
                / mw_ca
            )

        @self.Expression(
            self.flowsheet().time,
            doc="Calcium dimensionless molality",
        )
        def calcium_molality(blk, t):
            return (
                blk.calcium_molar_flow[t] * pyunits.s / pyunits.mol
                / blk.solvent_mass_flow[t] * pyunits.kg / pyunits.s
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Calcium log10 molality",
        )
        def eq_calcium_log10molality(blk, t):
            return (
                blk.properties_out[t].log10_molality_comp["Ca^2+"]
                == log_10(blk.calcium_molality[t])
            )

        @self.Expression(
            self.flowsheet().time,
            doc="Iron molar flow (mol/s)",
        )
        def iron_molar_flow(blk, t):
            return (
                pyunits.convert(blk.properties_in[t].flow_mass["Fe"], to_units=pyunits.kg / pyunits.s,)
                / mw_fe
            )

        @self.Expression(
            self.flowsheet().time,
            doc="Iron dimensionless molality",
        )
        def iron_molality(blk, t):
            return (
                blk.iron_molar_flow[t] * pyunits.s / pyunits.mol
                / blk.solvent_mass_flow[t] * pyunits.kg / pyunits.s
            )

        # Assume all iron forms iron(III) ions
        @self.Constraint(
            self.flowsheet().time,
            doc="Iron log10 molality",
        )
        def eq_iron_log10molality(blk, t):
            return (
                blk.properties_out[t].log10_molality_comp["Fe^3+"]
                == log_10(blk.iron_molality[t])
            )

        @self.Expression(
            self.flowsheet().time,
            doc="Cerium molar flow (mol/s)",
        )
        def cerium_molar_flow(blk, t):
            return (
                pyunits.convert(blk.properties_in[t].flow_mass["Ce"], to_units=pyunits.kg / pyunits.s,)
                / mw_ce
            )

        @self.Expression(
            self.flowsheet().time,
            doc="Cerium dimensionless molality",
        )
        def cerium_molality(blk, t):
            return (
                blk.cerium_molar_flow[t] * pyunits.s / pyunits.mol
                / blk.solvent_mass_flow[t] * pyunits.kg / pyunits.s
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Cerium log10 molality",
        )
        def eq_cerium_log10molality(blk, t):
            return (
                blk.properties_out[t].log10_molality_comp["Ce^3+"]
                == log_10(blk.cerium_molality[t])
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Oxalate log10 molality",
        )
        def eq_oxalate_log10molality(blk, t):
            return (
                blk.properties_out[t].log10_molality_comp["C2O4^2-"]
                == -4
            )

        self.components = Set(
            initialize=[
                "HC2O4^-",
                "H2C2O4",
                "Ce(OH)^2+",
                "Ce(OH)2^+",
                "Ce(OH)3",
                "Ce(OH)4^-",
                "Ce(C2O4)^+",
                "Ce(C2O4)2^-",
                "Ce(C2O4)3^3-",
                "Al(OH)^2+",
                "Al(OH)2^+",
                "Al(OH)3",
                "Al(OH)4^-",
                "Al2(OH)2^4+",
                "Al3(OH)4^5+",
                "Al(C2O4)^+",
                "Al(C2O4)2^-",
                "Al(C2O4)3^3-",
                "Al(HC2O4)^2+",
                "Al(OH)(C2O4)",
                "Al(OH)2(C2O4)^-",
                "Al(OH)(HC2O4)2^2-",
                "Fe(OH)^2+",
                "Fe(OH)2^+",
                "Fe(OH)3",
                "Fe(OH)4^-",
                "Fe2(OH)2^4+",
                "Fe3(OH)4^5+",
                "Fe(C2O4)^+",
                "Fe(C2O4)2^-",
                "Fe(C2O4)3^3-",
                "Ca(OH)^+",
                "Ca(C2O4)",
            ]
        )

        @self.Constraint(
            self.flowsheet().time,
            self.components,
            doc="Components",
        )
        def eq_component_log10molality(blk, t, i):
            return (
                blk.properties_out[t].log10_molality_comp[i]
                == -20
            )

        # iscale.set_scaling_factor(self.properties_out[0].flow_mass, 1e5)

        #TODO: Should components other than Fe^3+, Ca^2+, and Al^3+ be ignored/removed from the precip unit model?

        #TODO: Need to create pH adjustment block and track pH in SX stream instead of just fixing pH (as done below)
        pH = 2.5

        @self.Constraint(
            self.flowsheet().time,
            doc="Hydrogen log10 molality",
        )
        def eq_hydrogen_molality(blk, t):
            return (
                blk.properties_out[t].log10_molality_comp["H^+"]
                == log_10(10**(-pH))
            )

        # water dissociation equilibrium constant at 25C
        Kw = 1e-14

        @self.Constraint(
            self.flowsheet().time,
            doc="Hydroxide log10 molality",
        )
        def eq_hydroxide_molality(blk, t):
            return (
                blk.properties_out[t].log10_molality_comp["OH^-"]
                == log_10(Kw / 10**(-pH))
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
