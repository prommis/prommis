#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Membrane unit model.

Modification of the IDAES Multi-stream contactor unit.
(Addition of ports and membrane performance constraints)
"""

from pyomo.common.config import ConfigValue

# Pyomo import
from pyomo.environ import Var, exp, units

from idaes.core import declare_process_block_class
from idaes.core.util.exceptions import ConfigurationError

# IDAES imports
from idaes.models.unit_models.mscontactor import MSContactorData

__author__ = "Jason Yao"


@declare_process_block_class("Membrane")
class MembraneData(MSContactorData):
    """Modification of the MSContactor unit.

    Contains additional side stream ports.
    Membrane length, flux, and sieving coefficient constraints are added.
    """

    CONFIG = MSContactorData.CONFIG()

    CONFIG.declare(
        "flux",
        ConfigValue(domain=float, doc="Nominal value of membrane flux parameter"),
    )

    CONFIG.declare(
        "sieving_coefficient",
        ConfigValue(domain=dict, doc="Dictionary of sieving coefficient values"),
    )

    def build(self):
        """Build Membrane Model."""
        super().build()
        self._verify_membrane_inputs()
        self.add_side_stream_ports()

        # get solutes from membrane retentate (feed side) stream
        solutes = [
            i
            for i in self.config.streams.retentate.property_package.component_list
            if i != "solvent"
        ]

        self._verify_sieving_coefficients(solutes)

        self.add_membrane_performance_variables(solutes)
        self.add_membrane_constraints(solutes)

    def _verify_membrane_inputs(self):
        """Check that all required arguments are given."""
        if self.config.number_of_finite_elements is None:
            raise ConfigurationError(
                "Membrane model must be provided with number of elements"
            )
        if self.config.flux is None:
            raise ConfigurationError(
                "Membrane model must be provided with a " + "flux [m^3 / m^2 h] value"
            )
        if self.config.sieving_coefficient is None:
            raise ConfigurationError(
                "Membrane model must be provided with a dictionary of "
                + "sieving coefficient values"
            )
        required_streams = ["permeate", "retentate"]
        for stream in self.config.streams:
            if stream not in required_streams:
                raise ConfigurationError(
                    "Membrane model must be provided with only permeate "
                    + "and retentate streams"
                )

    def _verify_sieving_coefficients(self, solutes):
        """Check that sieving coefficients are given for all solutes."""
        if len(self.config.sieving_coefficient) != len(solutes):
            raise ConfigurationError(
                "Number of sieving coefficient values must match number of "
                "solutes in property package"
            )
        # also check to make sure sieving coefficients match each solute
        for sol in self.config.sieving_coefficient:
            if sol not in solutes:
                raise ConfigurationError(
                    "Sieving coefficient must match solutes " + "in property package"
                )

    def add_side_stream_ports(self):
        """Add ports for all side streams."""
        for stream, pconfig in self.config.streams.items():
            if pconfig.side_streams is not None:
                ss_block = getattr(self, stream + "_side_stream_state")
                for i in self.elements:
                    ss_element = self.elements.at(i)
                    ss_port, _ = ss_block.build_port(
                        f"{stream} Side Stream {i} Port",
                        slice_index=(slice(None), ss_element),
                    )
                    self.add_component(stream + f"_side_stream_{i}", ss_port)

    def add_membrane_performance_variables(self, solutes):
        """Add membrane length, flux, and sieving coefficients."""
        # set small initial length for initialization
        self.length = Var(units=units.m, bounds=(0.1, 10000), initialize=100)
        # set width as 1m
        self.width = Var(units=units.m, initialize=1)
        self.flux = Var(
            self.elements, units=units.m / units.hour, initialize=self.config.flux
        )
        self.sieving_coefficient = Var(
            solutes,
            self.elements,
            units=units.dimensionless,
            initialize={
                (sol, ele): self.config.sieving_coefficient[sol]
                for sol in solutes
                for ele in self.elements
            },
        )

        # fix values
        self.length.fixed = True
        self.width.fixed = True

        for idx in self.flux:
            self.flux[idx].fixed = True
        for idx in self.sieving_coefficient:
            self.sieving_coefficient[idx].fixed = True

    def add_membrane_constraints(self, solutes):
        """Add solute sieving, solvent flux, and LB/UB constraints."""

        # add flow lower bounds
        # TODO we need to be careful of membrane length and initialization
        # with this constraint, since depending on how much flow and membrane
        # length we initialize at (and how much membrane flux is set),
        # the stage cut LB/UB might be violated.
        @self.Constraint(["permeate", "retentate"])
        def flow_limits(b, loc):
            stage_cut_lb = 0.01
            stage_cut_ub = 0.99
            if loc == "permeate":
                frac = stage_cut_lb
            elif loc == "retentate":
                frac = 1 - stage_cut_ub
            last_ele = self.elements.at(-1)
            return getattr(b, f"{loc}")[0, last_ele].flow_vol >= frac * (
                b.retentate_inlet_state[0].flow_vol
                + sum(
                    b.retentate_side_stream_state[0, ele].flow_vol
                    for ele in self.elements
                )
            )

        # also adding concentration upper bounds
        @self.Constraint(solutes, ["permeate", "retentate"])
        def conc_limits(b, sol, loc):
            last_ele = self.elements.at(-1)
            if sol == "Li":
                conc_ub = 20 * units.kg / units.m**3  # kg/m^3
            elif sol == "Co":
                conc_ub = 200 * units.kg / units.m**3  # kg/m^3
            else:
                conc_ub = 200 * units.kg / units.m**3  # generic UB set to 200 kg/m^3
            return (
                conc_ub * getattr(b, f"{loc}")[0, last_ele].flow_vol
                >= getattr(b, f"{loc}")[0, last_ele].flow_mass_solute[sol]
            )

        @self.Constraint(self.elements)
        def solvent_rule(b, ele):
            return (
                b.material_transfer_term[0, ele, "permeate", "retentate", "solvent"]
                == b.flux[ele]
                * b.length
                * b.width
                * self.config.streams.retentate.property_package.dens_H2O
                / self.config.number_of_finite_elements
            )

        # Isolate nonlinear LN into LB/UB to make original sieving eqn linear
        # add new variables for each substitution
        # set initial value of variables to be reasonable small
        self.LN_M_in = Var(solutes, self.elements, initialize=5, bounds=(-20, 15))
        self.LN_M_out = Var(solutes, self.elements, initialize=5, bounds=(-20, 15))
        self.LN_F_in = Var(self.elements, initialize=5, bounds=(-20, 15))
        self.LN_F_out = Var(self.elements, initialize=5, bounds=(-20, 15))

        # set these as exponents to remove the logs
        #######################################################################
        @self.Constraint(solutes, self.elements)
        def LN_M_in_exp(b, sol, ele):
            if ele == 1:
                m_in = (
                    b.retentate_inlet_state[0].flow_mass_solute[sol]
                    + b.retentate_side_stream_state[0, ele].flow_mass_solute[sol]
                )
            else:
                ele_prev = b.elements.prev(ele)
                m_in = (
                    b.retentate[0, ele_prev].flow_mass_solute[sol]
                    + b.retentate_side_stream_state[0, ele].flow_mass_solute[sol]
                )
            return (m_in * units.hour / units.kg) == exp(b.LN_M_in[sol, ele])

        #######################################################################

        #######################################################################
        @self.Constraint(solutes, self.elements)
        def LN_M_out_exp(b, sol, ele):
            return (
                b.retentate[0, ele].flow_mass_solute[sol] * units.hour / units.kg
            ) == exp(b.LN_M_out[sol, ele])

        #######################################################################

        #######################################################################
        @self.Constraint(self.elements)
        def LN_F_in_exp(b, ele):
            if ele == 1:
                q_in = (
                    b.retentate_inlet_state[0].flow_vol
                    + b.retentate_side_stream_state[0, ele].flow_vol
                )
            else:
                ele_prev = b.elements.prev(ele)
                q_in = (
                    b.retentate[0, ele_prev].flow_vol
                    + b.retentate_side_stream_state[0, ele].flow_vol
                )
            return (q_in * units.hour / units.m**3) == exp(b.LN_F_in[ele])

        #######################################################################

        #######################################################################
        @self.Constraint(self.elements)
        def LN_F_out_exp(b, ele):
            return (b.retentate[0, ele].flow_vol * units.hour / units.m**3) == exp(
                b.LN_F_out[ele]
            )

        #######################################################################

        # linearized original solute sieving equation
        @self.Constraint(solutes, self.elements)
        def solute_rule_lin(b, sol, ele):
            return b.LN_M_out[sol, ele] - b.LN_M_in[sol, ele] == b.sieving_coefficient[
                sol, ele
            ] * (b.LN_F_out[ele] - b.LN_F_in[ele])
