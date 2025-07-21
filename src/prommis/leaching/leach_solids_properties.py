#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Initial property package for West Kentucky No. 13 coal refuse.

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


# -----------------------------------------------------------------------------
# Leach solution property package
@declare_process_block_class("CoalRefuseParameters")
class CoalRefuseParameterData(PhysicalParameterBlock):
    """
    Solid phase property package for West Kentucky No. 13 coal waste.

    Based on assay provided in:

    RESEARCH PERFORMANCE FINAL REPORT, Pilot-Scale Testing of an Integrated
    Circuit for the Extraction of Rare Earth Minerals and Elements from Coal
    and Coal Byproducts Using Advanced Separation Technologies,
    Honaker, R.Q., et al., DE-FE0027035

    Includes the following components:

    * Inerts
    * Rare Earth Oxides: Sc2O3, Y2O3, La2O3, Ce2O3, Pr2O3, Nd2O3, Sm2O3, Gd2O3, Dy2O3
    * Impurities: Al2O3, CaO, Fe2O3

    """

    def build(self):
        super().build()

        self.solid = Phase()

        # Solids
        self.inerts = Component()

        # REEs
        self.Sc2O3 = Component()
        self.Y2O3 = Component()
        self.La2O3 = Component()
        self.Ce2O3 = Component()
        self.Pr2O3 = Component()
        self.Nd2O3 = Component()
        self.Sm2O3 = Component()
        self.Gd2O3 = Component()
        self.Dy2O3 = Component()

        # Contaminants
        self.Al2O3 = Component()
        self.CaO = Component()
        self.Fe2O3 = Component()

        self.mw = Param(
            self.component_list,
            units=units.kg / units.mol,
            initialize={
                "inerts": 60.08e-3,
                "Sc2O3": (44.946 * 2 + 3 * 15.999) * 1e-3,
                "Y2O3": (88.905 * 2 + 3 * 15.999) * 1e-3,
                "La2O3": (138.905 * 2 + 3 * 15.999) * 1e-3,
                "Ce2O3": (140.116 * 2 + 3 * 15.999) * 1e-3,
                "Pr2O3": (140.907 * 2 + 3 * 15.999) * 1e-3,
                "Nd2O3": (144.242 * 2 + 3 * 15.999) * 1e-3,
                "Sm2O3": (150.36 * 2 + 3 * 15.999) * 1e-3,
                "Gd2O3": (157.25 * 2 + 3 * 15.999) * 1e-3,
                "Dy2O3": (162.50 * 2 + 3 * 15.999) * 1e-3,
                "Al2O3": (26.982 * 2 + 3 * 15.999) * 1e-3,
                "CaO": (40.078 + 15.999) * 1e-3,
                "Fe2O3": (55.845 * 2 + 3 * 15.999) * 1e-3,
            },
        )

        self.mass_frac_comp_initial = Param(
            self.component_list,
            units=units.kg / units.kg,
            initialize={
                "inerts": 0.6952,
                "Al2O3": 0.237,
                "Fe2O3": 0.0642,
                "CaO": 3.31e-3,
                "Sc2O3": 2.77966e-05,
                "Y2O3": 3.28653e-05,
                "La2O3": 6.77769e-05,
                "Ce2O3": 0.000156161,
                "Pr2O3": 1.71438e-05,
                "Nd2O3": 6.76618e-05,
                "Sm2O3": 1.47926e-05,
                "Gd2O3": 1.0405e-05,
                "Dy2O3": 7.54827e-06,
            },
        )

        self.dens_mass = Param(
            units=units.kg / units.litre,
            initialize=2.4,
            mutable=True,
        )

        self._state_block_class = CoalRefuseStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": units.hour,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


class _CoalRefuseStateBlock(StateBlock):
    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        # Fix state variables
        fix_state_vars(self)

        for sbd in self.values():
            if not sbd.config.defined_state:
                sbd.sum_mass_frac.deactivate()


@declare_process_block_class("CoalRefuseStateBlock", block_class=_CoalRefuseStateBlock)
class CoalRefuseStateBlockData(StateBlockData):
    """
    State block for solid West Kentucky No. 13 coal waste.

    """

    def build(self):
        super().build()

        self.flow_mass = Var(
            units=units.kg / units.hour,
            bounds=(1e-8, None),
        )
        self.mass_frac_comp = Var(
            self.params.component_list,
            units=units.kg / units.kg,
            bounds=(1e-8, None),
        )

        # TODO: This really should have a _comp suffix
        self.conversion = Var(
            self.params.component_list,
            initialize=0,
            units=units.dimensionless,
            bounds=(None, 0.999999),
        )

        @self.Constraint(self.params.component_list)
        def conversion_eq(b, j):
            if j == "inerts":
                return b.conversion[j] == 0
            return (1 - b.conversion[j]) * b.params.mass_frac_comp_initial[
                j
            ] * b.mass_frac_comp["inerts"] == b.mass_frac_comp[
                j
            ] * b.params.mass_frac_comp_initial[
                "inerts"
            ]

        if not self.config.defined_state:
            self.sum_mass_frac = Constraint(
                expr=1 == sum(self.mass_frac_comp[j] for j in self.component_list)
            )

    def get_material_flow_terms(self, p, j):
        return units.convert(
            self.flow_mass * self.mass_frac_comp[j] / self.params.mw[j],
            to_units=units.mol / units.hour,
        )

    def get_material_density_terms(self, p, j):
        return units.convert(
            self.params.dens_mass * self.mass_frac_comp[j] / self.params.mw[j],
            to_units=units.mol / units.m**3,
        )

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        return {
            "flow_mass": self.flow_mass,
            "mass_frac_comp": self.mass_frac_comp,
        }
