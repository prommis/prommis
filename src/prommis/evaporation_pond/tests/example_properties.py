#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Example property package for testing evaporation pond unit model.

Authors: Andrew Lee
"""

from pyomo.environ import Constraint, Param, Var, units

from idaes.core import (
    Component,
    MaterialFlowBasis,
    Phase,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    declare_process_block_class,
)
from idaes.core.util.initialization import fix_state_vars
from idaes.core.util.misc import add_object_reference


# -----------------------------------------------------------------------------
# Leach solution property package
@declare_process_block_class("BrineParameters")
class BrineParameterData(PhysicalParameterBlock):
    """
    Simple property package for lithum brine solutions.

    """

    def build(self):
        super().build()

        self.liquid = Phase()

        # Solvent
        self.H2O = Component()

        # Solutes
        self.Li = Component()
        self.Na = Component()
        self.Cl = Component()

        self.mw = Param(
            self.component_list,
            units=units.kg / units.mol,
            initialize={
                "H2O": 18e-3,
                "Li": 6.94e-3,
                "Na": 22.99e-3,
                "Cl": 35.45e-3,
            },
        )

        # Assume dilute acid, density of pure water
        self.dens_mass = Param(
            initialize=1,
            units=units.kg / units.litre,
            mutable=True,
        )

        self._state_block_class = BrineStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "flow_vol": {"method": None},
                "conc_mass_comp": {"method": None},
                "dens_mass": {"method": "_dens_mass"},
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


class _BrineStateBlock(StateBlock):
    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        # Fix state variables
        fix_state_vars(self)


@declare_process_block_class("BrineStateBlock", block_class=_BrineStateBlock)
class BrineStateBlockData(StateBlockData):
    """
    State block for lithium brine solutions.

    """

    def build(self):
        super().build()

        self.flow_vol = Var(
            units=units.L / units.hour,
            bounds=(1e-8, None),
        )
        self.conc_mass_comp = Var(
            self.params.component_list,
            units=units.mg / units.L,
            bounds=(1e-20, None),
        )

        @self.Expression(self.params.component_list, doc="Molar concentration")
        def conc_mole_comp(b, j):
            return units.convert(
                b.conc_mass_comp[j] / b.mw[j], to_units=units.mol / units.L
            )

        if not self.config.defined_state:
            # Concentration of H2O based on assumed density
            self.h2o_concentration = Constraint(
                expr=self.conc_mass_comp["H2O"] == 1e6 * units.mg / units.L
            )

    @property
    def mw(self):
        """Molecular weight of species"""
        return self.params.mw

    def _dens_mass(self):
        add_object_reference(self, "dens_mass", self.params.dens_mass)

    def get_material_flow_terms(self, _, j):
        # Note conversion to mol/hour
        if j == "H2O":
            # Assume constant density of 1 kg/L
            return self.flow_vol * self.params.dens_mass / self.params.mw[j]
        # Need to convert from moles to mass
        return units.convert(
            self.flow_vol * self.conc_mass_comp[j] / self.params.mw[j],
            to_units=units.mol / units.hour,
        )

    def get_material_density_terms(self, _, j):
        if j == "H2O":
            return units.convert(
                self.params.dens_mass / self.params.mw[j],
                to_units=units.mol / units.m**3,
            )
        return units.convert(
            self.conc_mass_comp[j] / self.params.mw[j],
            to_units=units.mol / units.m**3,
        )

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        return {
            "flow_vol": self.flow_vol,
            "conc_mass_comp": self.conc_mass_comp,
        }
