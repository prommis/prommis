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
Translator block converting from old_leaching properties to solvent extraction properties.
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
)

__author__ = "Marcus Holly"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("Translator_leaching_SX")
class TranslatorDataLeachingSX(TranslatorData):
    """
    Translator block representing the old_leaching/SX interface
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
        super(TranslatorDataLeachingSX, self).build()

        @self.Constraint(
            self.flowsheet().time,
            doc="Equality volumetric flow equation",
        )
        def eq_flow_vol_rule(blk, t):
            return blk.properties_out[t].flow_vol == blk.properties_in[t].flow_vol

        # @self.Constraint(
        #     self.flowsheet().time,
        #     doc="Equality temperature equation",
        # )
        # def eq_temperature_rule(blk, t):
        #     return blk.properties_out[t].temperature == blk.properties_in[t].temperature

        # @self.Constraint(
        #     self.flowsheet().time,
        #     doc="Equality pressure equation",
        # )
        # def eq_pressure_rule(blk, t):
        #     return blk.properties_out[t].pressure == blk.properties_in[t].pressure

        self.metals = Set(
            initialize=["Al", "Ca", "Fe", "Sc", "Y", "La", "Ce", "Pr", "Nd", "Sm", "Gd", "Dy"]
        )

        @self.Constraint(
            self.flowsheet().time,
            self.unchanged_component,
            doc="Equality equation for metal components",
        )
        def eq_metal_mass_flow(blk, t, i):
            return (
                blk.properties_out[t].flow_mass[i]
                == blk.properties_in[t].get_material_flow_terms("Liq", i)
                * blk.properties_in[t].params.mw[i]
                * 1000 * pyunits.g / pyunits.kg
            )

        # self.acid_set = Set(
        #     initialize=[
        #         "H",
        #         "HSO4",
        #         "SO4",
        #     ]
        # )
        #
        # @self.Constraint(
        #     self.flowsheet().time,
        #     self.acid_set,
        #     doc="Components with no flow equation",
        # )
        # def eq_acid_set(blk, t, i):
        #     return blk.properties_out[t].flow_mass["H2SO4"] == sum(
        #         blk.properties_in[t].get_material_flow_terms("Liq", i) for i in blk.acid_set
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
