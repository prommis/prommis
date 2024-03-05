from pyomo.environ import ConcreteModel, Constraint, Set, SolverFactory, TransformationFactory
from pyomo.network import Arc, SequentialDecomposition

from idaes.core import FlowDirection, FlowsheetBlock, MaterialBalanceType, MomentumBalanceType
from idaes.core.initialization import InitializationStatus
from idaes.core.util.initialization import propagate_state
from idaes.core.initialization.block_triangularization import (
    BlockTriangularizationInitializer,
)
from idaes.core.util.model_statistics import degrees_of_freedom as dof
from idaes.models.unit_models.separator import (
    EnergySplittingType,
    Separator,
    SplittingType,
)
from idaes.models.unit_models.feed import Feed, FeedInitializer
from idaes.models.unit_models.mixer import Mixer, MixingType, MomentumMixingType

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.solvent_extraction import SolventExtraction
from prommis.precipitate.precipitate_liquid_properties import AqueousParameter
from idaes.core.util.model_diagnostics import DiagnosticsToolbox


m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.prop_o = REESolExOgParameters()
m.fs.leach_soln = LeachSolutionParameters()
m.fs.properties_aq = AqueousParameter()

m.fs.solex_rougher = SolventExtraction(
    number_of_finite_elements=3,
    dynamic=False,
    aqueous_stream={
        "property_package": m.fs.leach_soln,
        "flow_direction": FlowDirection.forward,
        "has_energy_balance": False,
        "has_pressure_balance": False,
    },
    organic_stream={
        "property_package": m.fs.prop_o,
        "flow_direction": FlowDirection.backward,
        "has_energy_balance": False,
        "has_pressure_balance": False,
    },
)

m.fs.aqueous_feed = Feed(property_package=m.fs.leach_soln)

m.fs.solex_cleaner = SolventExtraction(
    number_of_finite_elements=3,
    dynamic=False,
    aqueous_stream={
        "property_package": m.fs.properties_aq,
        "flow_direction": FlowDirection.forward,
        "has_energy_balance": False,
        "has_pressure_balance": False,
    },
    organic_stream={
        "property_package": m.fs.prop_o,
        "flow_direction": FlowDirection.backward,
        "has_energy_balance": False,
        "has_pressure_balance": False,
    },
)

m.fs.sep1 = Separator(
    property_package=m.fs.prop_o,
    outlet_list=["recycle", "purge"],
    split_basis=SplittingType.totalFlow,
    material_balance_type=MaterialBalanceType.componentTotal,
    momentum_balance_type=MomentumBalanceType.none,
    energy_split_basis=EnergySplittingType.none,
)

m.fs.sep2 = Separator(
    property_package=m.fs.prop_o,
    outlet_list=["make_up", "purge"],
    split_basis=SplittingType.totalFlow,
    material_balance_type=MaterialBalanceType.componentTotal,
    momentum_balance_type=MomentumBalanceType.none,
    energy_split_basis=EnergySplittingType.none,
)

m.fs.mixer = Mixer(
    property_package=m.fs.prop_o,
    num_inlets=2,
    inlet_list=["make_up", "recycle"],
    material_balance_type=MaterialBalanceType.componentTotal,
    energy_mixing_type=MixingType.none,
    momentum_mixing_type=MomentumMixingType.none,
)

m.fs.feed = Arc(
    source=m.fs.aqueous_feed.outlet,
    destination=m.fs.solex_rougher.mscontactor.aqueous_inlet,
)
m.fs.recycle = Arc(
    source=m.fs.mixer.outlet,
    destination=m.fs.solex_rougher.mscontactor.organic_inlet,
)
m.fs.organic_stream1 = Arc(
    source=m.fs.solex_rougher.mscontactor.organic_outlet,
    destination=m.fs.solex_cleaner.mscontactor.organic_inlet,
)
m.fs.organic_stream2 = Arc(
    source=m.fs.solex_cleaner.mscontactor.organic_outlet,
    destination=m.fs.sep1.inlet,
)
m.fs.organic_recycle = Arc(
    source=m.fs.sep1.recycle,
    destination=m.fs.mixer.recycle,
)
m.fs.organic_purge = Arc(
    source=m.fs.sep1.purge,
    destination=m.fs.sep2.inlet,
)
m.fs.organic_make_up = Arc(
    source=m.fs.sep2.make_up,
    destination=m.fs.mixer.make_up,
)


TransformationFactory("network.expand_arcs").apply_to(m)


# m.fs.aqueous_feed.conc_mass_comp.fix(1e-5)
m.fs.aqueous_feed.conc_mass_comp[0, "H2O"].fix(1e-9)
m.fs.aqueous_feed.conc_mass_comp[0, "H"].fix(1e-9)
m.fs.aqueous_feed.conc_mass_comp[0, "SO4"].fix(1e-9)
m.fs.aqueous_feed.conc_mass_comp[0, "HSO4"].fix(1e-9)

m.fs.aqueous_feed.conc_mass_comp[0, "Al"].fix(820)
m.fs.aqueous_feed.conc_mass_comp[0, "Ca"].fix(5230)
m.fs.aqueous_feed.conc_mass_comp[0, "Fe"].fix(270)
m.fs.aqueous_feed.conc_mass_comp[0, "Sc"].fix(209.31)
m.fs.aqueous_feed.conc_mass_comp[0, "Y"].fix(637.74)
m.fs.aqueous_feed.conc_mass_comp[0, "La"].fix(2032.77)
m.fs.aqueous_feed.conc_mass_comp[0, "Ce"].fix(4516.13)
m.fs.aqueous_feed.conc_mass_comp[0, "Pr"].fix(756.64)
m.fs.aqueous_feed.conc_mass_comp[0, "Nd"].fix(2047.85)
m.fs.aqueous_feed.conc_mass_comp[0, "Sm"].fix(369.1)
m.fs.aqueous_feed.conc_mass_comp[0, "Gd"].fix(174.38)
m.fs.aqueous_feed.conc_mass_comp[0, "Dy"].fix(101.12)

m.fs.aqueous_feed.flow_vol.fix(4.4)

m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Al"].fix(7.54e-10)
m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Ca"].fix(4.955e-9)
m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Fe"].fix(1.491e-7)
m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Sc"].fix(321.34)
m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Y"].fix(5.67e-6)
m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["La"].fix(1.78e-05)
m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Ce"].fix(4.019e-5)
m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Pr"].fix(6.73e-6)
m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Nd"].fix(1.82e-5)
m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Sm"].fix(3.285e-6)
m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Gd"].fix(1.55e-6)
m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Dy"].fix(9e-7)

m.fs.solex_rougher.mscontactor.organic_inlet_state[0].flow_vol.fix(62.01)

recycle_fraction = 0.95
purge_fraction = 0.50

m.fs.sep1.split_fraction[:, "recycle"].fix(recycle_fraction)

# TODO: Make sure make up is pure organic - use componentFlow?
m.fs.sep2.split_fraction[:, "recycle"].fix(purge_fraction)

# m.fs.sep2.split_fraction[0, "purge", "DEHPA"].fix(purge_fraction)
# m.fs.sep2.split_fraction[0, "purge", "Al"].fix(0.99)
# m.fs.sep2.split_fraction[0, "purge", "Ca"].fix(0.99)
# m.fs.sep2.split_fraction[0, "purge", "Fe"].fix(0.99)
# m.fs.sep2.split_fraction[0, "purge", "Sc"].fix(0.99)
# m.fs.sep2.split_fraction[0, "purge", "Y"].fix(0.99)
# m.fs.sep2.split_fraction[0, "purge", "La"].fix(0.99)
# m.fs.sep2.split_fraction[0, "purge", "Ce"].fix(0.99)
# m.fs.sep2.split_fraction[0, "purge", "Pr"].fix(0.99)
# m.fs.sep2.split_fraction[0, "purge", "Nd"].fix(0.99)
# m.fs.sep2.split_fraction[0, "purge", "Sm"].fix(0.99)
# m.fs.sep2.split_fraction[0, "purge", "Gd"].fix(0.99)
# m.fs.sep2.split_fraction[0, "purge", "Dy"].fix(0.99)

# print(dof(m))
dt = DiagnosticsToolbox(model=m)
dt.report_structural_issues()
dt.display_underconstrained_set()
dt.display_overconstrained_set()

# Initializing of the model

seq = SequentialDecomposition()
seq.options.tear_method = "Direct"
seq.options.iterLim = 1
seq.options.tear_set = [m.fs.recycle]

# Solving of the model

solver = SolverFactory("ipopt")
solver.solve(m, tee=True)

# m.fs.solex_rougher.mscontactor.organic_inlet_state[0].flow_vol.unfix()
#
# m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Al"].unfix()
# m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Ca"].unfix()
# m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Fe"].unfix()
# m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Sc"].unfix()
# m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Y"].unfix()
# m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["La"].unfix()
# m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Ce"].unfix()
# m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Pr"].unfix()
# m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Nd"].unfix()
# m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Sm"].unfix()
# m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Gd"].unfix()
# m.fs.solex_rougher.mscontactor.organic_inlet_state[0].conc_mass_comp["Dy"].unfix()
#
#
# dissolved_elements = Set(
#     initialize=[
#         "Al",
#         "Ca",
#         "Fe",
#         "Sc",
#         "Y",
#         "La",
#         "Ce",
#         "Pr",
#         "Nd",
#         "Sm",
#         "Gd",
#         "Dy",
#     ]
# )
#
#
# def rule_organic_mass_balance(b, t, j):
#     return (
#                 b.organic_inlet_state[t].flow_vol * b.organic_inlet_state[t].conc_mass_comp[j]
#                 == b.organic[t, 1].flow_vol * b.organic[t, 1].conc_mass_comp[j] * recycle_fraction
#         )
#
#
# m.fs.solex_rougher.mscontactor.organic_mass_balance_constraint = Constraint(
#     m.fs.time,
#     dissolved_elements,
#     rule=rule_organic_mass_balance,
# )
#
# m.fs.solex_cleaner.mscontactor.organic_mass_balance_constraint = Constraint(
#     m.fs.time,
#     dissolved_elements,
#     rule=rule_organic_mass_balance,
# )
#
# solver = SolverFactory("ipopt")
# solver.solve(m, tee=True)

# Final organic outlet display
# m.fs.solex_rougher.display()
# m.fs.solex_cleaner.display()
