# Import Pyomo libraries
from pyomo.environ import (
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

__author__ = "Arkoprabho Dasgupta"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("TranslatorLeachPrecip")
class TranslatorDataLeachPrecip(TranslatorData):
    """
    Translator block representing the SX/old_leaching interface
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

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality volumetric flow equation",
        )
        def eq_flow_vol_rule(blk, t):
            return blk.properties_out[t].flow_vol == blk.properties_in[t].flow_vol

        self.components = Set(
            initialize=[
                "Al",
                "Ca",
                "Fe",
                "Sc",
                "Y",
                "La",
                "Ce",
                "Pr",
                "Nd",
                "Sm",
                "Gd",
                "Dy",
                "H2O",
                "H",
                "SO4",
                "HSO4",
                "Cl",
            ]
        )

        @self.Constraint(
            self.flowsheet().time,
            self.components,
            doc="Equality equation for metal components",
        )
        def eq_conc_mass_metals(blk, t, i):
            return (
                blk.properties_out[t].conc_mass_comp[i]
                == blk.properties_in[t].conc_mass_comp[i]
            )

        # self.acids = Set(initialize=["H", "HSO4", "SO4"])

        # @self.Constraint(
        #     self.flowsheet().time,
        #     self.acids,
        #     doc="Equality equation for acid components",
        # )
        # def eq_conc_mole_acids(blk, t, i):
        #     return (
        #         blk.properties_out[t].conc_mole_acid[i]
        #         == 1e-7 * pyunits.mol / pyunits.L
        #     )

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
