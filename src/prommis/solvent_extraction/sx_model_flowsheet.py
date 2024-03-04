from pyomo.environ import ConcreteModel, Constraint, Set, SolverFactory
from pyomo.network import Arc

from idaes.core import FlowDirection, FlowsheetBlock, MaterialBalanceType, MomentumBalanceType
from idaes.core.initialization import InitializationStatus
from idaes.core.initialization.block_triangularization import (
    BlockTriangularizationInitializer,
)
from idaes.core.util.model_statistics import degrees_of_freedom as dof
from idaes.models.unit_models.separator import (
    EnergySplittingType,
    Separator,
    SplittingType,
)

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.solvent_extraction import SolventExtraction
from prommis.precipitate.precipitate_liquid_properties import AqueousParameter

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
    outlet_list=["make-up", "purge"],
    split_basis=SplittingType.componentFlow,
    material_balance_type=MaterialBalanceType.componentTotal,
    momentum_balance_type=MomentumBalanceType.none,
    energy_split_basis=EnergySplittingType.none,
)

dissolved_elements = Set(
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
    ]
)

recycle_fraction = 0.95
purge_fraction = 0.50


def rule_organic_mass_balance(b, t, j):
    return (
                b.organic_inlet_state[t].flow_vol * b.organic_inlet_state[t].conc_mass_comp[j] * recycle_fraction
                == b.organic[t, 1].flow_vol * b.organic[t, 1].conc_mass_comp[j]
        )


m.fs.solex_rougher.mscontactor.organic_mass_balance_constraint = Constraint(
    m.fs.time,
    dissolved_elements,
    rule=rule_organic_mass_balance,
)

m.fs.solex_cleaner.mscontactor.organic_mass_balance_constraint = Constraint(
    m.fs.time,
    dissolved_elements,
    rule=rule_organic_mass_balance,
)

m.fs.organic_stream1 = Arc(
    source=m.fs.solex_rougher.mscontactor.organic_outlet,
    destination=m.fs.solex_cleaner.mscontactor.organic_inlet,
)
m.fs.organic_stream2 = Arc(
    source=m.fs.solex_cleaner.mscontactor.organic_outlet,
    destination=m.fs.solex_rougher.mscontactor.organic_inlet,
)
m.fs.organic_purge = Arc(
    source=m.fs.solex_cleaner.mscontactor.organic_outlet,
    destination=m.fs.solex_rougher.mscontactor.organic_inlet,
)
m.fs.liq_feed = Arc(
    source=m.fs.leach_liquid_feed.outlet, destination=m.fs.leach_mixer.feed
)

m.fs.solex_rougher.mscontactor.aqueous_inlet_state[0].conc_mass_comp["H2O"].fix(1e-9)
m.fs.solex_rougher.mscontactor.aqueous_inlet_state[0].conc_mass_comp["H"].fix(1e-9)
m.fs.solex_rougher.mscontactor.aqueous_inlet_state[0].conc_mass_comp["SO4"].fix(1e-9)
m.fs.solex_rougher.mscontactor.aqueous_inlet_state[0].conc_mass_comp["HSO4"].fix(1e-9)
m.fs.solex_rougher.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Al"].fix(820)
m.fs.solex_rougher.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Ca"].fix(5230)
m.fs.solex_rougher.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Fe"].fix(270)
m.fs.solex_rougher.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Sc"].fix(209.31)
m.fs.solex_rougher.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Y"].fix(637.74)
m.fs.solex_rougher.mscontactor.aqueous_inlet_state[0].conc_mass_comp["La"].fix(2032.77)
m.fs.solex_rougher.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Ce"].fix(4516.13)
m.fs.solex_rougher.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Pr"].fix(756.64)
m.fs.solex_rougher.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Nd"].fix(2047.85)
m.fs.solex_rougher.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Sm"].fix(369.1)
m.fs.solex_rougher.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Gd"].fix(174.38)
m.fs.solex_rougher.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Dy"].fix(101.12)

m.fs.solex_rougher.mscontactor.aqueous_inlet_state[0].flow_vol.fix(4.4)

m.fs.sep1.split_fraction[:, "recycle"].fix(recycle_fraction)

m.fs.sep2.split_fraction[0, "recycle", "DEHPA"].fix(purge_fraction)
m.fs.sep2.split_fraction[0, "recycle", "Al"].fix(0.99)
m.fs.sep2.split_fraction[0, "recycle", "Ca"].fix(0.99)
m.fs.sep2.split_fraction[0, "recycle", "Fe"].fix(0.99)
m.fs.sep2.split_fraction[0, "recycle", "Sc"].fix(0.99)
m.fs.sep2.split_fraction[0, "recycle", "Y"].fix(0.99)
m.fs.sep2.split_fraction[0, "recycle", "La"].fix(0.99)
m.fs.sep2.split_fraction[0, "recycle", "Ce"].fix(0.99)
m.fs.sep2.split_fraction[0, "recycle", "Pr"].fix(0.99)
m.fs.sep2.split_fraction[0, "recycle", "Nd"].fix(0.99)
m.fs.sep2.split_fraction[0, "recycle", "Sm"].fix(0.99)
m.fs.sep2.split_fraction[0, "recycle", "Gd"].fix(0.99)
m.fs.sep2.split_fraction[0, "recycle", "Dy"].fix(0.99)

# m.fs.solex.mscontactor.organic_inlet_state[0].flow_vol.fix(62.01)

print(dof(m))

# Initializing of the model

initializer = BlockTriangularizationInitializer()
initializer.initialize(m.fs.solex)
assert initializer.summary[m.fs.solex]["status"] == InitializationStatus.Ok

# Solving of the model

solver = SolverFactory("ipopt")
solver.options["bound_push"] = 1e-8
solver.options["mu_init"] = 1e-8
solver.solve(m, tee=True)

# Final organic outlet display
m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp.display()
m.fs.solex.mscontactor.organic[0, 1].conc_mol_comp.display()

# Final aqueous outlets display
m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp.display()
m.fs.solex.mscontactor.aqueous[0, 3].conc_mol_comp.display()
