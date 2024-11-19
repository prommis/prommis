#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Class for building the full IDAES diafiltration flowsheet."""

# Pyomo imports
import pyomo.environ as pyo
from pyomo.network import Arc

# IDAES imports
from idaes.core import FlowsheetBlock, MaterialBalanceType, MomentumBalanceType
from idaes.models.unit_models import MSContactorInitializer
from idaes.models.unit_models import (
    Mixer,
    MixingType,
    MomentumMixingType,
    MixerInitializer,
)
from idaes.models.unit_models import (
    Separator as Splitter,
    EnergySplittingType,
    SeparatorInitializer,
)
from idaes.core.util.initialization import propagate_state

# Custom imports
from solute_property import SoluteParameters
from membrane import Membrane
from precipitator import Precipitator

# other imports
import idaes.logger as idaeslog
import logging
import numpy as np


class diafiltration_model:
    """
    Multi-stage diafiltration model, built using custom IDAES units.

    Model specifications:
    NS stages
    NT tubes

    Parameters:
    solutes
    flux
    sieving coefficient
    feed
    diafiltrate
    """

    def __init__(
        self,
        NS,
        NT,
        solutes=["Li", "Co", "Ni"],
        flux=0.1,
        sieving_coefficient={"Li": 1.3, "Co": 0.5, "Ni": 0.5},
        feed={
            "solvent": 100,  # m^3/hr of water
            "Li": 1.7 * 100,  # kg/hr
            "Co": 17 * 100,  # kg/hr
            "Ni": 1e-8,
        },
        diaf={
            "solvent": 30,  # m^3/hr of water
            "Li": 0.1 * 30,  # kg/hr
            "Co": 0.2 * 30,  # kg/hr
            "Ni": 1e-8,
        },
        precipitate=False,
        # perc_precipitate={'permeate': 0.81, 'retentate': 0.99}
        perc_precipitate={
            "permeate": {"Li": 0.81, "Co": 0.05, "Ni": 0.05},
            "retentate": {"Li": 0.05, "Co": 0.99, "Ni": 0.81},
        },
    ):
        """Store model parameters."""
        self.ns = NS
        self.nt = NT
        self.solutes = solutes
        self.flux = flux
        self.sieving_coefficient = sieving_coefficient
        self.feed = feed
        self.diaf = diaf
        self.precipitate = precipitate
        self.perc_precipitate = perc_precipitate

    def build_flowsheet(self, mixing="tube"):
        """Build the multi-stage diafiltration flowsheet."""
        m = pyo.ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = SoluteParameters(solutes=self.solutes)
        m.fs.stages = pyo.RangeSet(self.ns)
        m.fs.tubes = pyo.RangeSet(self.nt)
        m.fs.solutes = pyo.Set(initialize=self.solutes)
        self.add_stages(m)
        self.add_mixers(m, mixing=mixing)
        self.add_recycles(m, mixing=mixing)
        self.add_feed(m)
        self.add_diafiltrate(m)
        self.add_isotropic_stage_length(m)
        self.add_objectives(m)

        if self.precipitate:
            self.mix_products(m)
            self.add_precipitators(m)
            self.connect_diafiltrate(m)
            self.precipitate_objectives(m)

        # resolve solver tolerance issues by setting LB of all units outside of
        # membranes to 0
        for i in m.component_data_objects(pyo.Var):
            # i.setlb(0)
            # if i.lb == 1e-8:
            #     i.setlb(0)
            if i.lb == 1e-8 and i.name.split(".")[1].split("[")[0] != "stage":
                i.setlb(0)
        #     elif i.lb == 1e-8:
        #         i.setlb(1e-2)

        return m

    def add_stages(self, m):
        """Add membrane stages."""
        # add membrane stages
        m.fs.stage = Membrane(
            pyo.RangeSet(self.ns),
            number_of_finite_elements=self.nt,
            flux=self.flux,
            sieving_coefficient=self.sieving_coefficient,
            streams={
                "retentate": {
                    "property_package": m.fs.properties,
                    "side_streams": pyo.RangeSet(self.nt),
                    "has_energy_balance": False,
                    "has_pressure_balance": False,
                },
                "permeate": {
                    "property_package": m.fs.properties,
                    "has_feed": False,
                    "has_energy_balance": False,
                    "has_pressure_balance": False,
                },
            },
        )

        # fix inlet of first stage to its LB 1e-8 (unused)
        m.fs.stage[1].retentate_inlet_state[0].flow_vol.fix(1e-8)
        for sol in self.solutes:
            m.fs.stage[1].retentate_inlet_state[0].mass_solute[sol].fix(1e-8)

        # add retentate and permeate splitters for each stage
        m.fs.split_retentate = Splitter(
            pyo.RangeSet(self.ns),
            outlet_list=["product", "recycle"],
            property_package=m.fs.properties,
            material_balance_type=MaterialBalanceType.componentTotal,
            momentum_balance_type=MomentumBalanceType.none,
            energy_split_basis=EnergySplittingType.none,
        )

        m.fs.split_permeate = Splitter(
            pyo.RangeSet(self.ns),
            outlet_list=["product", "forward"],
            property_package=m.fs.properties,
            material_balance_type=MaterialBalanceType.componentTotal,
            momentum_balance_type=MomentumBalanceType.none,
            energy_split_basis=EnergySplittingType.none,
        )

        # set up stage connections
        for i in pyo.RangeSet(self.ns):
            # connect splitters to each stage
            m.fs.add_component(
                f"connect_retentate_split_{i}",
                Arc(
                    source=m.fs.stage[i].retentate_outlet,
                    destination=m.fs.split_retentate[i].inlet,
                ),
            )
            m.fs.add_component(
                f"connect_permeate_split_{i}",
                Arc(
                    source=m.fs.stage[i].permeate_outlet,
                    destination=m.fs.split_permeate[i].inlet,
                ),
            )

            # connect subsequent stages
            if i != 1:
                m.fs.add_component(
                    f"connect_stages_{i-1}_and_{i}",
                    Arc(
                        source=m.fs.split_permeate[i - 1].forward,
                        destination=m.fs.stage[i].retentate_inlet,
                    ),
                )

        pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    def add_mixers(self, m, mixing="tube"):
        """Add side stream mixers."""
        if mixing == "tube":
            m.fs.inlet_mixers = Mixer(
                pyo.RangeSet(self.ns),
                pyo.RangeSet(self.nt),
                inlet_list=["feed", "diafiltrate", "recycle"],
                property_package=m.fs.properties,
                material_balance_type=MaterialBalanceType.componentTotal,
                energy_mixing_type=MixingType.none,
                momentum_mixing_type=MomentumMixingType.none,
            )

            # connect each mixer to a stage side stream
            for i in pyo.RangeSet(self.ns):
                for j in pyo.RangeSet(self.nt):
                    m.fs.add_component(
                        f"connect_stage_{i}_inlet_side_stream_{j}",
                        Arc(
                            source=m.fs.inlet_mixers[i, j].outlet,
                            destination=getattr(
                                m.fs.stage[i], f"retentate_side_stream_{j}"
                            ),
                        ),
                    )

                    # fix recycle inlet of last stage mixers to LB (unused)
                    if i == self.ns:
                        m.fs.inlet_mixers[i, j].recycle.flow_vol[0].fix(1e-8)
                        for sol in self.solutes:
                            m.fs.inlet_mixers[i, j].recycle.mass_solute[0, sol].fix(
                                1e-8
                            )

        elif mixing == "stage":
            m.fs.inlet_mixers = Mixer(
                pyo.RangeSet(self.ns),
                inlet_list=["feed", "diafiltrate", "recycle"],
                property_package=m.fs.properties,
                material_balance_type=MaterialBalanceType.componentTotal,
                energy_mixing_type=MixingType.none,
                momentum_mixing_type=MomentumMixingType.none,
            )

            # split the stage mixer into each stage side stream with splitter
            m.fs.splitters = Splitter(
                pyo.RangeSet(self.ns),
                num_outlets=self.nt,
                property_package=m.fs.properties,
                material_balance_type=MaterialBalanceType.componentTotal,
                momentum_balance_type=MomentumBalanceType.none,
                energy_split_basis=EnergySplittingType.none,
            )

            # connect each split stream to a stage side stream
            for i in pyo.RangeSet(self.ns):
                m.fs.add_component(
                    f"connect_inlet_splitter_{i}",
                    Arc(
                        source=m.fs.inlet_mixers[i].outlet,
                        destination=m.fs.splitters[i].inlet,
                    ),
                )
                # fix recycle inlet of last stage mixers to LB (unused)
                if i == self.ns:
                    m.fs.inlet_mixers[i].recycle.flow_vol[0].fix(1e-8)
                    for sol in self.solutes:
                        m.fs.inlet_mixers[i].recycle.mass_solute[0, sol].fix(1e-8)

                for j in pyo.RangeSet(self.nt):
                    m.fs.add_component(
                        f"connect_stage_{i}_inlet_side_stream_{j}",
                        Arc(
                            source=getattr(m.fs.splitters[i], f"outlet_{j}"),
                            destination=getattr(
                                m.fs.stage[i], f"retentate_side_stream_{j}"
                            ),
                        ),
                    )
        else:
            raise NotImplementedError(
                "Diafiltration flowsheet only supports 'stage' and 'tube' "
                "mixing superstructures."
            )

        pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    def add_recycles(self, m, mixing="tube"):
        """Add recycle streams."""
        if self.ns > 1:
            if mixing == "tube":
                m.fs.recycle_splitters = Splitter(
                    pyo.RangeSet(2, self.ns),
                    num_outlets=self.nt,
                    property_package=m.fs.properties,
                    material_balance_type=MaterialBalanceType.componentTotal,
                    momentum_balance_type=MomentumBalanceType.none,
                    energy_split_basis=EnergySplittingType.none,
                )

                # connect each recycle stream to the previous stage
                for i in pyo.RangeSet(2, self.ns):
                    m.fs.add_component(
                        f"connect_recycle_splitter_{i}",
                        Arc(
                            source=m.fs.split_retentate[i].recycle,
                            destination=m.fs.recycle_splitters[i].inlet,
                        ),
                    )
                    for j in pyo.RangeSet(self.nt):
                        m.fs.add_component(
                            f"connect_recycle_{i}_split_{j}",
                            Arc(
                                source=getattr(
                                    m.fs.recycle_splitters[i], f"outlet_{j}"
                                ),
                                destination=m.fs.inlet_mixers[i - 1, j].recycle,
                            ),
                        )

            elif mixing == "stage":
                # connect recycle stream directly to previous stage mixer
                for i in pyo.RangeSet(2, self.ns):
                    m.fs.add_component(
                        f"connect_recycle_{i}",
                        Arc(
                            source=m.fs.split_retentate[i].recycle,
                            destination=m.fs.inlet_mixers[i - 1].recycle,
                        ),
                    )
            else:
                raise NotImplementedError(
                    "Diafiltration flowsheet only supports 'stage' and 'tube' "
                    "mixing superstructures."
                )

            pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    def add_feed(self, m):
        """Add feed streams."""
        m.fs.split_feed = Splitter(
            num_outlets=len(m.fs.inlet_mixers),
            property_package=m.fs.properties,
            material_balance_type=MaterialBalanceType.componentTotal,
            momentum_balance_type=MomentumBalanceType.none,
            energy_split_basis=EnergySplittingType.none,
        )

        # TODO change this to directly use state block instead of reference
        # TODO look for other instances with this issue.
        # set feed values
        m.fs.split_feed.inlet.flow_vol[0].fix(self.feed["solvent"])
        for sol in self.solutes:
            m.fs.split_feed.inlet.mass_solute[0, sol].fix(self.feed[sol])

        # connect feeds to inlet mixers
        for i in m.fs.inlet_mixers:
            if type(i) is tuple:
                num = (i[0] - 1) * self.nt + i[1]
            else:
                num = i
            m.fs.add_component(
                f"connect_feed_{num}",
                Arc(
                    source=getattr(m.fs.split_feed, f"outlet_{num}"),
                    destination=m.fs.inlet_mixers[i].feed,
                ),
            )

        pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    def add_diafiltrate(self, m):
        """Add diafiltrate streams."""
        m.fs.split_diafiltrate = Splitter(
            num_outlets=len(m.fs.inlet_mixers),
            property_package=m.fs.properties,
            material_balance_type=MaterialBalanceType.componentTotal,
            momentum_balance_type=MomentumBalanceType.none,
            energy_split_basis=EnergySplittingType.none,
        )

        # set diafiltrate values
        m.fs.split_diafiltrate.inlet.flow_vol[0].fix(self.diaf["solvent"])
        for sol in self.solutes:
            m.fs.split_diafiltrate.inlet.mass_solute[0, sol].fix(self.diaf[sol])

        # connect diafiltrates to inlet mixers
        for i in m.fs.inlet_mixers:
            if type(i) is tuple:
                num = (i[0] - 1) * self.nt + i[1]
            else:
                num = i
            m.fs.add_component(
                f"connect_diafiltrate_{num}",
                Arc(
                    source=getattr(m.fs.split_diafiltrate, f"outlet_{num}"),
                    destination=m.fs.inlet_mixers[i].diafiltrate,
                ),
            )

        pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    def add_isotropic_stage_length(self, m):
        """Add isotropic stage length rule."""

        @m.fs.Constraint(pyo.RangeSet(2, self.ns))
        def isotropic_rule(b, i):
            return b.stage[1].length == b.stage[i].length

    def add_objectives(self, m):
        """Add recovery objectives."""

        # mass recovery of solutes
        @m.Expression()
        def rec_mass_co(b):
            return sum(
                m.fs.split_retentate[i].product.mass_solute[0, "Co"]
                for i in pyo.RangeSet(self.ns)
            )

        @m.Expression()
        def rec_mass_li(b):
            return sum(
                m.fs.split_permeate[i].product.mass_solute[0, "Li"]
                for i in pyo.RangeSet(self.ns)
            )

        # percent recovery of solutes
        @m.Expression()
        def rec_perc_co(b):
            return b.rec_mass_co / (
                (self.feed["Co"] + self.diaf["Co"]) * pyo.units.kg / pyo.units.hour
            )

        @m.Expression()
        def rec_perc_li(b):
            return b.rec_mass_li / (
                (self.feed["Li"] + self.diaf["Li"]) * pyo.units.kg / pyo.units.hour
            )

        # mass recovery objectives
        m.co_obj = pyo.Objective(expr=m.rec_mass_co, sense=pyo.maximize)

        m.li_obj = pyo.Objective(expr=m.rec_mass_li, sense=pyo.maximize)
        m.li_obj.deactivate()  # normally using Co as objective

        # percent recovery lower bound constraints
        m.R = pyo.Param(initialize=0.95, mutable=True)

        @m.Constraint()
        def li_lb(b):
            return b.rec_perc_li >= m.R

        @m.Constraint()
        def co_lb(b):
            return b.rec_perc_co >= m.R

        m.co_lb.deactivate()  # normally use Co as objective

        # product stream purities
        @m.Expression()
        def purity_co(b):
            return sum(
                m.fs.split_retentate[i].product.mass_solute[0, "Co"]
                for i in pyo.RangeSet(self.ns)
            ) / sum(
                m.fs.split_retentate[i].product.mass_solute[0, j]
                for i in pyo.RangeSet(self.ns)
                for j in self.solutes
            )

        @m.Expression()
        def purity_li(b):
            return sum(
                m.fs.split_permeate[i].product.mass_solute[0, "Li"]
                for i in pyo.RangeSet(self.ns)
            ) / sum(
                m.fs.split_permeate[i].product.mass_solute[0, j]
                for i in pyo.RangeSet(self.ns)
                for j in self.solutes
            )

        m.pure = pyo.Param(initialize=0.1, mutable=True)

        @m.Constraint()
        def purity_co_lb(b):
            return sum(
                m.fs.split_retentate[i].product.mass_solute[0, "Co"]
                for i in pyo.RangeSet(self.ns)
            ) >= m.pure * sum(
                m.fs.split_retentate[i].product.mass_solute[0, j]
                for i in pyo.RangeSet(self.ns)
                for j in self.solutes
            )

        @m.Constraint()
        def purity_li_lb(b):
            return sum(
                m.fs.split_permeate[i].product.mass_solute[0, "Li"]
                for i in pyo.RangeSet(self.ns)
            ) >= m.pure * sum(
                m.fs.split_permeate[i].product.mass_solute[0, j]
                for i in pyo.RangeSet(self.ns)
                for j in self.solutes
            )

        m.purity_co_lb.deactivate()
        m.purity_li_lb.deactivate()

        # product stream mass impurities
        @m.Expression()
        def impurity_retentate(b):
            return sum(
                m.fs.split_retentate[i].product.mass_solute[0, "Li"]
                for i in pyo.RangeSet(self.ns)
            )

        @m.Expression()
        def impurity_permeate(b):
            return sum(
                m.fs.split_permeate[i].product.mass_solute[0, "Co"]
                for i in pyo.RangeSet(self.ns)
            )

    def mix_products(self, m):
        """Mix process products."""
        m.fs.mix_product_retentate = Mixer(
            num_inlets=self.ns,
            property_package=m.fs.properties,
            material_balance_type=MaterialBalanceType.componentTotal,
            energy_mixing_type=MixingType.none,
            momentum_mixing_type=MomentumMixingType.none,
        )

        m.fs.mix_product_permeate = Mixer(
            num_inlets=self.ns,
            property_package=m.fs.properties,
            material_balance_type=MaterialBalanceType.componentTotal,
            energy_mixing_type=MixingType.none,
            momentum_mixing_type=MomentumMixingType.none,
        )
        # connect product streams to product mixer
        for i in pyo.RangeSet(self.ns):
            m.fs.add_component(
                f"connect_prod_retentate_{i}",
                Arc(
                    source=m.fs.split_retentate[i].product,
                    destination=getattr(m.fs.mix_product_retentate, f"inlet_{i}"),
                ),
            )
            m.fs.add_component(
                f"connect_prod_permeate_{i}",
                Arc(
                    source=m.fs.split_permeate[i].product,
                    destination=getattr(m.fs.mix_product_permeate, f"inlet_{i}"),
                ),
            )
        pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    def add_precipitators(self, m):
        """Add precipitator units."""
        m.fs.precipitator = Precipitator(
            ["retentate", "permeate"],
            outlet_list=["solid", "recycle"],
            yields=self.perc_precipitate,
            property_package=m.fs.properties,
            material_balance_type=MaterialBalanceType.componentTotal,
            momentum_balance_type=MomentumBalanceType.none,
            energy_split_basis=EnergySplittingType.none,
        )

        # connect mixed products
        for prod in ["retentate", "permeate"]:
            m.fs.add_component(
                f"connect_precipitator_{prod}",
                Arc(
                    source=getattr(m.fs, f"mix_product_{prod}").outlet,
                    destination=m.fs.precipitator[prod].inlet,
                ),
            )

            # set split fraction
            # for sol in self.solutes:
            # m.fs.precipitator[prod].yields[sol, 'solid']\
            #                        .fix(self.perc_precipitate[prod][sol])

        pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    def connect_diafiltrate(self, m):
        """Connect precipitator recycle to diafiltrate inlet."""
        # add mixer and splitter
        m.fs.mix_precipitate_recycle = Mixer(
            inlet_list=["retentate", "permeate"],
            property_package=m.fs.properties,
            material_balance_type=MaterialBalanceType.componentTotal,
            energy_mixing_type=MixingType.none,
            momentum_mixing_type=MomentumMixingType.none,
        )
        m.fs.split_precipitate_recycle = Splitter(
            outlet_list=["waste", "recycle"],
            property_package=m.fs.properties,
            material_balance_type=MaterialBalanceType.componentTotal,
            momentum_balance_type=MomentumBalanceType.none,
            energy_split_basis=EnergySplittingType.none,
        )

        # connect precipitator outlets to mixer
        for prod in ["retentate", "permeate"]:
            m.fs.add_component(
                f"connect_precipitator_{prod}_mixer",
                Arc(
                    source=m.fs.precipitator[prod].recycle,
                    destination=getattr(m.fs.mix_precipitate_recycle, f"{prod}"),
                ),
            )
        m.fs.connect_precipitate_mix_split = Arc(
            source=m.fs.mix_precipitate_recycle.outlet,
            destination=m.fs.split_precipitate_recycle.inlet,
        )

        # connect to diafiltrate inlet
        m.fs.connect_precipitate_diafiltrate = Arc(
            source=m.fs.split_precipitate_recycle.recycle,
            destination=m.fs.split_diafiltrate.inlet,
        )

        # unfix diafiltrate inlet
        m.fs.split_diafiltrate.inlet.flow_vol[0].unfix()
        for sol in self.solutes:
            m.fs.split_diafiltrate.inlet.mass_solute[0, sol].unfix()

        # This should be the same as the set diafiltrate UB
        m.fs.split_precipitate_recycle.split_fraction[0, "recycle"].fix(
            self.diaf["solvent"] / (self.feed["solvent"] + self.diaf["solvent"])
        )

        pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    def precipitate_objectives(self, m):
        """Set up objectives for flowsheet with precipitators."""
        m.co_obj.deactivate()
        m.li_lb.deactivate()

        @m.Expression()
        def prec_mass_co(b):
            return m.fs.precipitator["retentate"].solid.mass_solute[0, "Co"]

        @m.Expression()
        def prec_mass_li(b):
            return m.fs.precipitator["permeate"].solid.mass_solute[0, "Li"]

        # percent recovery of solutes
        @m.Expression()
        def prec_perc_co(b):
            return b.prec_mass_co / (
                (m.fs.split_feed.mixed_state[0].mass_solute["Co"])
                * pyo.units.kg
                / pyo.units.hour
            )

        @m.Expression()
        def prec_perc_li(b):
            return b.prec_mass_li / (
                (m.fs.split_feed.mixed_state[0].mass_solute["Li"])
                * pyo.units.kg
                / pyo.units.hour
            )

        m.prec_co_obj = pyo.Objective(expr=m.prec_mass_co, sense=pyo.maximize)
        m.prec_li_obj = pyo.Objective(expr=m.prec_mass_li, sense=pyo.maximize)
        m.prec_li_obj.deactivate()
        m.prec_li_lb = pyo.Constraint(expr=m.prec_perc_li >= m.R)
        m.Rco = pyo.Param(initialize=0.95, mutable=True)
        m.prec_co_lb = pyo.Constraint(expr=m.prec_perc_co >= m.Rco)
        m.prec_co_lb.deactivate()

    def initialize_stage_splitters(self, m, stage, mixing="tube"):
        """Initialize splitters of a stage."""
        split_initializer = SeparatorInitializer()
        i = stage
        if i == 1:
            m.fs.split_retentate[i].split_fraction[0, "recycle"].fix(0)
            m.fs.split_permeate[i].split_fraction[0, "product"].fix(0.5)
        elif i == self.ns:
            m.fs.split_retentate[i].split_fraction[0, "product"].fix(0.5)
            m.fs.split_permeate[i].split_fraction[0, "forward"].fix(0)
        else:
            m.fs.split_retentate[i].split_fraction[0, "product"].fix(0.5)
            m.fs.split_permeate[i].split_fraction[0, "product"].fix(0.5)

        split_initializer.initialize(m.fs.split_retentate[i])
        split_initializer.initialize(m.fs.split_permeate[i])

        # propagate recycle if not first stage
        if i != 1:
            if mixing == "tube":
                propagate_state(
                    source=m.fs.split_retentate[i].recycle,
                    destination=m.fs.recycle_splitters[i].inlet,
                )
                for j in pyo.RangeSet(self.nt - 1):
                    m.fs.recycle_splitters[i].split_fraction[0, f"outlet_{j}"].fix(
                        1 / self.nt
                    )
                split_initializer.initialize(m.fs.recycle_splitters[i])

                # propagate recycle to the previous stage
                for j in pyo.RangeSet(self.nt):
                    m.fs.inlet_mixers[i - 1, j].recycle.flow_vol[0].unfix()
                    for sol in self.solutes:
                        m.fs.inlet_mixers[i - 1, j].recycle.mass_solute[0, sol].unfix()

                    propagate_state(
                        source=getattr(m.fs.recycle_splitters[i], f"outlet_{j}"),
                        destination=m.fs.inlet_mixers[i - 1, j].recycle,
                    )
            elif mixing == "stage":
                propagate_state(
                    source=m.fs.split_retentate[i].recycle,
                    destination=m.fs.inlet_mixers[i - 1].recycle,
                )

    def initialize_single_stage(self, m, stage, mixing="tube"):
        """Initialize a single stage."""
        # set up initializers
        stage_initializer = MSContactorInitializer()

        # propagate side stream inlets
        i = stage
        for j in pyo.RangeSet(self.nt):
            m.fs.stage[i].retentate_side_stream_state[0, j].flow_vol.unfix()
            for sol in self.solutes:
                m.fs.stage[i].retentate_side_stream_state[0, j].mass_solute[sol].unfix()
            if mixing == "tube":
                propagate_state(
                    source=m.fs.inlet_mixers[i, j].outlet,
                    destination=getattr(m.fs.stage[i], f"retentate_side_stream_{j}"),
                )
            elif mixing == "stage":
                propagate_state(
                    source=getattr(m.fs.splitters[i], f"outlet_{j}"),
                    destination=getattr(m.fs.stage[i], f"retentate_side_stream_{j}"),
                )
            m.fs.stage[i].retentate_side_stream_state[0, j].flow_vol.fix()
            for sol in self.solutes:
                m.fs.stage[i].retentate_side_stream_state[0, j].mass_solute[sol].fix()

            # the mscontactor initializer sets the initial value of the stream
            # to the sum of its inlets. Set material transfer to permeate to a
            # nonzero value to avoid warnings of permeate values being lower
            # than the LB of 1e-8.
            m.fs.stage[i].material_transfer_term[
                0, j, "permeate", "retentate", "solvent"
            ].set_value(1e-4)
            # m.fs.stage[i].permeate[0.0, j].flow_vol.setlb(0)
            for sol in self.solutes:
                m.fs.stage[i].material_transfer_term[
                    0, j, "permeate", "retentate", sol
                ].set_value(1e-4)
                # m.fs.stage[i].permeate[0.0, j].mass_solute[sol].setlb(0)
        # propagate retentate inlet if not first stage
        if i != 1:
            propagate_state(
                source=m.fs.split_permeate[i - 1].forward,
                destination=m.fs.stage[i].retentate_inlet,
            )

        # initialize membrane unit and associated product splitters
        stage_initializer.initialize(m.fs.stage[i])

        propagate_state(
            source=m.fs.stage[i].retentate_outlet,
            destination=m.fs.split_retentate[i].inlet,
        )
        propagate_state(
            source=m.fs.stage[i].permeate_outlet,
            destination=m.fs.split_permeate[i].inlet,
        )

        self.initialize_stage_splitters(m, stage, mixing)

    def initialize_precipitators(self, m, mixing, precipitate):
        """Initialize precipitators."""
        mixer_initializer = MixerInitializer()
        split_initializer = SeparatorInitializer()

        # initialize precipitator related units
        if precipitate:
            for i in pyo.RangeSet(self.ns):
                propagate_state(
                    source=m.fs.split_retentate[i].product,
                    destination=getattr(m.fs.mix_product_retentate, f"inlet_{i}"),
                )
                propagate_state(
                    source=m.fs.split_permeate[i].product,
                    destination=getattr(m.fs.mix_product_permeate, f"inlet_{i}"),
                )
            mixer_initializer.initialize(m.fs.mix_product_retentate)
            mixer_initializer.initialize(m.fs.mix_product_permeate)

            for prod in ["retentate", "permeate"]:
                propagate_state(
                    source=getattr(m.fs, f"mix_product_{prod}").outlet,
                    destination=m.fs.precipitator[prod].inlet,
                )
                split_initializer.initialize(m.fs.precipitator[prod])

                propagate_state(
                    source=m.fs.precipitator[prod].recycle,
                    destination=getattr(m.fs.mix_precipitate_recycle, f"{prod}"),
                )
            mixer_initializer.initialize(m.fs.mix_precipitate_recycle)
            propagate_state(
                source=m.fs.mix_precipitate_recycle.outlet,
                destination=m.fs.split_precipitate_recycle.inlet,
            )
            split_initializer.initialize(m.fs.split_precipitate_recycle)
            propagate_state(
                source=m.fs.split_precipitate_recycle.recycle,
                destination=m.fs.split_diafiltrate.inlet,
            )

    def check_precipitates(self, m, precipitate):
        """Check to see if precipitator exists."""
        if not precipitate:
            if hasattr(m.fs, "precipitator"):
                raise NotImplementedError(
                    "Precipitator was found in model. "
                    "Try setting its initialize option to `True`."
                )

    def initialize_feeds(self, m, mixing):
        """Set inlets and and initialize feed."""
        if mixing == "tube":
            inlets = self.ns * self.nt
        elif mixing == "stage":
            inlets = self.ns
        for i in pyo.RangeSet(inlets - 1):
            m.fs.split_feed.split_fraction[0, f"outlet_{i}"].fix(1 / inlets)
            m.fs.split_diafiltrate.split_fraction[0, f"outlet_{i}"].fix(1 / inlets)
        split_initializer = SeparatorInitializer()
        split_initializer.initialize(m.fs.split_feed)

    def initialize(self, m, mixing="tube", precipitate=False):
        """Initialize all IDAES unit models."""
        # initialize feed and diafiltrate splitters
        # 1. we put a bit of the feed streams into all inlet locations
        # 2. then we use fixed-point iteration to converge recycles

        # print initialization header
        # TODO could incorporate initialization info into idaes logger
        print("*" * 50)
        print("* Diafiltration Flowsheet Initialization")
        print("*" * 50)

        # Change logger level to avoid flood of initialization INFO
        logger = logging.getLogger("idaes.init")
        logger.setLevel(idaeslog.WARNING)

        self.initialize_feeds(m, mixing)
        split_initializer = SeparatorInitializer()

        # outer loop for precipitator initialization
        diaf_old = 0
        diaf_new = 1
        itr = 0
        check_sol = self.solutes[0]
        while not np.isclose(diaf_old, diaf_new) and itr < 100:
            diaf_old = (
                m.fs.split_diafiltrate.mixed_state[0].mass_solute[check_sol].value
            )
            split_initializer.initialize(m.fs.split_diafiltrate)

            # inner loop for membranes initialization
            check_old = 0
            check_new = 1
            iteration = 0
            # check first mixer
            check_loc = m.fs.inlet_mixers.index_set().first()
            while not np.isclose(check_old, check_new) and iteration < 100:
                print(f"  Membranes Iteration: {iteration}")
                # propagate feed and diafiltrate state to mixers
                for i in m.fs.inlet_mixers:
                    if type(i) is tuple:
                        num = (i[0] - 1) * self.nt + i[1]
                    else:
                        num = i
                    propagate_state(
                        source=getattr(m.fs.split_feed, f"outlet_{num}"),
                        destination=m.fs.inlet_mixers[i].feed,
                    )
                    propagate_state(
                        source=getattr(m.fs.split_diafiltrate, f"outlet_{num}"),
                        destination=m.fs.inlet_mixers[i].diafiltrate,
                    )

                    if iteration == 0:
                        # fix recycle streams as 0
                        if m.fs.inlet_mixers[i].recycle.flow_vol[0].value != 1e-8:
                            m.fs.inlet_mixers[i].recycle.flow_vol[0].fix(0)
                            for sol in self.solutes:
                                m.fs.inlet_mixers[i].recycle.mass_solute[0, sol].fix(0)
                    else:
                        if i == check_loc:
                            check_old = (
                                m.fs.inlet_mixers[i].mixed_state[0].flow_vol.value
                            )

                    mixer_initializer = MixerInitializer()
                    mixer_initializer.initialize(m.fs.inlet_mixers[i])
                    if i == check_loc:
                        check_new = m.fs.inlet_mixers[i].mixed_state[0].flow_vol.value

                    # propagate state if using stage mixing to the splitter
                    if mixing == "stage":
                        propagate_state(
                            source=m.fs.inlet_mixers[i].outlet,
                            destination=m.fs.splitters[i].inlet,
                        )
                        for j in pyo.RangeSet(self.nt - 1):
                            m.fs.splitters[i].split_fraction[0, f"outlet_{j}"].fix(
                                1 / self.nt
                            )
                        split_initializer.initialize(m.fs.splitters[i])

                # initialize stages and relevant splitters
                for i in pyo.RangeSet(self.ns):
                    self.initialize_single_stage(m, stage=i, mixing=mixing)

                print(f"  =>Recycle Error: {check_new - check_old}")
                iteration += 1

            ###############################################################
            # check if precipitators are included
            self.check_precipitates(m, precipitate)

            if not precipitate:
                break

            print(f"Precipitator Iteration {itr}")

            self.initialize_precipitators(m, mixing, precipitate)

            diaf_new = (
                m.fs.split_diafiltrate.mixed_state[0].mass_solute[check_sol].value
            )
            print(f"=>Recycle Error: {diaf_new - diaf_old}\n")
            itr += 1

    def num_inlets(self, mixing):
        """Find the number of inlets."""
        if mixing == "tube":
            inlets = self.ns * self.nt
        elif mixing == "stage":
            inlets = self.ns
        return inlets

    def unfix_dof(self, m, mixing="tube", precipitate=False):
        """Unfix model DoF."""
        inlets = self.num_inlets(mixing)
        # go through all split feed/diafiltrate streams and unfix them
        for i in pyo.RangeSet(inlets - 1):
            m.fs.split_feed.split_fraction[0, f"outlet_{i}"].unfix()
            m.fs.split_diafiltrate.split_fraction[0, f"outlet_{i}"].unfix()

        # Unfix stage recycle inlets
        for i in m.fs.inlet_mixers:
            if m.fs.inlet_mixers[i].recycle.flow_vol[0].value != 1e-8:
                m.fs.inlet_mixers[i].recycle.flow_vol[0].unfix()
                for sol in self.solutes:
                    m.fs.inlet_mixers[i].recycle.mass_solute[0, sol].unfix()
            if mixing == "stage":
                for j in pyo.RangeSet(self.nt - 1):
                    m.fs.splitters[i].split_fraction[0, f"outlet_{j}"].unfix()

        # go through each stage and unfix relevant membrane units and splitters
        for i in pyo.RangeSet(self.ns):
            for j in pyo.RangeSet(self.nt):
                m.fs.stage[i].retentate_side_stream_state[0, j].flow_vol.unfix()
                for sol in self.solutes:
                    m.fs.stage[i].retentate_side_stream_state[0, j].mass_solute[
                        sol
                    ].unfix()
            m.fs.split_retentate[i].split_fraction[0, "product"].unfix()
            m.fs.split_permeate[i].split_fraction[0, "product"].unfix()

            # unfix recycle splitter
            if i != 1:
                if mixing == "tube":
                    for j in pyo.RangeSet(self.nt - 1):
                        m.fs.recycle_splitters[i].split_fraction[
                            0, f"outlet_{j}"
                        ].unfix()

            # unfix stage lengths
            m.fs.stage[i].length.unfix()

        # the recycle and forward stream of respective first/last stages are
        # not used. Fix these to 0.
        m.fs.split_retentate[1].split_fraction[0, "recycle"].fix(0)
        m.fs.split_permeate[self.ns].split_fraction[0, "forward"].fix(0)

        # unfix recycle split fraction and set diafiltrate upper bound
        # the upper bound is set to the saved diafiltrate flow rate
        self.check_precipitates(m, precipitate)
        if precipitate:
            m.fs.split_precipitate_recycle.split_fraction[0, "recycle"].unfix()
            m.fs.split_diafiltrate.inlet.flow_vol.setub(self.diaf["solvent"])
            m.fs.precipitator["retentate"].V.unfix()
            m.fs.precipitator["permeate"].V.unfix()
