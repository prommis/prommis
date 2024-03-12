from pyomo.environ import ConcreteModel, SolverFactory, TransformationFactory, Suffix
from pyomo.network import Arc, SequentialDecomposition

from idaes.core import FlowDirection, FlowsheetBlock, MaterialBalanceType, MomentumBalanceType
from idaes.core.initialization.block_triangularization import (
    BlockTriangularizationInitializer,
)
from idaes.models.unit_models.separator import (
    EnergySplittingType,
    Separator,
    SplittingType,
)
from idaes.models.unit_models.mixer import Mixer, MixingType, MomentumMixingType
from idaes.models.unit_models.product import Product
from idaes.models.unit_models.feed import Feed, FeedInitializer

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.solvent_extraction import SolventExtraction

from idaes.core.util.initialization import propagate_state


def main():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.prop_o = REESolExOgParameters()
    m.fs.leach_soln = LeachSolutionParameters()

    m.fs.aqueous_feed = Feed(property_package=m.fs.leach_soln)

    m.fs.organic_make_up = Feed(property_package=m.fs.prop_o)

    m.fs.solex_rougher1 = SolventExtraction(
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
        aqueous_to_organic=True,
    )

    m.fs.solex_rougher1.partition_coefficient[1, "aqueous", "organic", "Al"] = 5.2 / 100
    m.fs.solex_rougher1.partition_coefficient[1, "aqueous", "organic", "Ca"] = 3.0 / 100
    m.fs.solex_rougher1.partition_coefficient[1, "aqueous", "organic", "Fe"] = 24.7 / 100
    m.fs.solex_rougher1.partition_coefficient[1, "aqueous", "organic", "Sc"] = 99.9 / 100
    m.fs.solex_rougher1.partition_coefficient[1, "aqueous", "organic", "Y"] = 99.9 / 100
    m.fs.solex_rougher1.partition_coefficient[1, "aqueous", "organic", "La"] = 32.4 / 100
    m.fs.solex_rougher1.partition_coefficient[1, "aqueous", "organic", "Ce"] = 58.2 / 100
    m.fs.solex_rougher1.partition_coefficient[1, "aqueous", "organic", "Pr"] = 58.2 / 100
    m.fs.solex_rougher1.partition_coefficient[1, "aqueous", "organic", "Nd"] = 87.6 / 100
    m.fs.solex_rougher1.partition_coefficient[1, "aqueous", "organic", "Sm"] = 99.9 / 100
    m.fs.solex_rougher1.partition_coefficient[1, "aqueous", "organic", "Gd"] = 69.8 / 100
    m.fs.solex_rougher1.partition_coefficient[1, "aqueous", "organic", "Dy"] = 96.6 / 100

    m.fs.solex_rougher1.partition_coefficient[2, "aqueous", "organic", "Al"] = 4.9 / 100
    m.fs.solex_rougher1.partition_coefficient[2, "aqueous", "organic", "Ca"] = 12.3 / 100
    m.fs.solex_rougher1.partition_coefficient[2, "aqueous", "organic", "Fe"] = 6.4 / 100
    m.fs.solex_rougher1.partition_coefficient[2, "aqueous", "organic", "Sc"] = 16.7 / 100
    m.fs.solex_rougher1.partition_coefficient[2, "aqueous", "organic", "Y"] = 99.9 / 100
    m.fs.solex_rougher1.partition_coefficient[2, "aqueous", "organic", "La"] = 23.2 / 100
    m.fs.solex_rougher1.partition_coefficient[2, "aqueous", "organic", "Ce"] = 24.9 / 100
    m.fs.solex_rougher1.partition_coefficient[2, "aqueous", "organic", "Pr"] = 15.1 / 100
    m.fs.solex_rougher1.partition_coefficient[2, "aqueous", "organic", "Nd"] = 99.9 / 100
    m.fs.solex_rougher1.partition_coefficient[2, "aqueous", "organic", "Sm"] = 99.9 / 100
    m.fs.solex_rougher1.partition_coefficient[2, "aqueous", "organic", "Gd"] = 7.6 / 100
    m.fs.solex_rougher1.partition_coefficient[2, "aqueous", "organic", "Dy"] = 5.0 / 100

    m.fs.solex_rougher1.partition_coefficient[3, "aqueous", "organic", "Al"] = 4.9 / 100
    m.fs.solex_rougher1.partition_coefficient[3, "aqueous", "organic", "Ca"] = 12.3 / 100
    m.fs.solex_rougher1.partition_coefficient[3, "aqueous", "organic", "Fe"] = 6.4 / 100
    m.fs.solex_rougher1.partition_coefficient[3, "aqueous", "organic", "Sc"] = 16.7 / 100
    m.fs.solex_rougher1.partition_coefficient[3, "aqueous", "organic", "Y"] = 99.9 / 100
    m.fs.solex_rougher1.partition_coefficient[3, "aqueous", "organic", "La"] = 23.2 / 100
    m.fs.solex_rougher1.partition_coefficient[3, "aqueous", "organic", "Ce"] = 24.9 / 100
    m.fs.solex_rougher1.partition_coefficient[3, "aqueous", "organic", "Pr"] = 15.1 / 100
    m.fs.solex_rougher1.partition_coefficient[3, "aqueous", "organic", "Nd"] = 99.9 / 100
    m.fs.solex_rougher1.partition_coefficient[3, "aqueous", "organic", "Sm"] = 99.9 / 100
    m.fs.solex_rougher1.partition_coefficient[3, "aqueous", "organic", "Gd"] = 7.6 / 100
    m.fs.solex_rougher1.partition_coefficient[3, "aqueous", "organic", "Dy"] = 5.0 / 100

    m.fs.leach_recycle = Product(property_package=m.fs.leach_soln)

    m.fs.acid_feed1 = Feed(property_package=m.fs.leach_soln)

    m.fs.solex_rougher2 = SolventExtraction(
        number_of_finite_elements=1,
        dynamic=False,
        aqueous_stream={
            "property_package": m.fs.leach_soln,
            "flow_direction": FlowDirection.backward,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        organic_stream={
            "property_package": m.fs.prop_o,
            "flow_direction": FlowDirection.forward,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        aqueous_to_organic=False,
    )
    m.fs.solex_rougher2.partition_coefficient[1, "aqueous", "organic", "Al"] = (100 - 0.12) / 100
    m.fs.solex_rougher2.partition_coefficient[1, "aqueous", "organic", "Ca"] = (100 - 0.55) / 100
    m.fs.solex_rougher2.partition_coefficient[1, "aqueous", "organic", "Fe"] = (100 - 0.007) / 100
    m.fs.solex_rougher2.partition_coefficient[1, "aqueous", "organic", "Sc"] = (100 - 99.9) / 100
    m.fs.solex_rougher2.partition_coefficient[1, "aqueous", "organic", "Y"] = (100 - 99.9) / 100
    m.fs.solex_rougher2.partition_coefficient[1, "aqueous", "organic", "La"] = (100 - 99.8) / 100
    m.fs.solex_rougher2.partition_coefficient[1, "aqueous", "organic", "Ce"] = (100 - 99.9) / 100
    m.fs.solex_rougher2.partition_coefficient[1, "aqueous", "organic", "Pr"] = (100 - 99.9) / 100
    m.fs.solex_rougher2.partition_coefficient[1, "aqueous", "organic", "Nd"] = (100 - 99.9) / 100
    m.fs.solex_rougher2.partition_coefficient[1, "aqueous", "organic", "Sm"] = (100 - 99.9) / 100
    m.fs.solex_rougher2.partition_coefficient[1, "aqueous", "organic", "Gd"] = (100 - 99.9) / 100
    m.fs.solex_rougher2.partition_coefficient[1, "aqueous", "organic", "Dy"] = (100 - 99.9) / 100

    m.fs.leach_recycle2 = Product(property_package=m.fs.leach_soln)

    m.fs.acid_feed2 = Feed(property_package=m.fs.leach_soln)

    m.fs.solex_rougher3 = SolventExtraction(
        number_of_finite_elements=2,
        dynamic=False,
        aqueous_stream={
            "property_package": m.fs.leach_soln,
            "flow_direction": FlowDirection.backward,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        organic_stream={
            "property_package": m.fs.prop_o,
            "flow_direction": FlowDirection.forward,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        aqueous_to_organic=False,
    )
    m.fs.solex_rougher3.partition_coefficient[:, "aqueous", "organic", "Al"] = (100 - 0.5) / 100
    m.fs.solex_rougher3.partition_coefficient[:, "aqueous", "organic", "Ca"] = (100 - 0.5) / 100
    m.fs.solex_rougher3.partition_coefficient[:, "aqueous", "organic", "Fe"] = (100 - 0.5) / 100
    m.fs.solex_rougher3.partition_coefficient[:, "aqueous", "organic", "Sc"] = (100 - 98.5) / 100
    m.fs.solex_rougher3.partition_coefficient[:, "aqueous", "organic", "Y"] = (100 - 0.5) / 100
    m.fs.solex_rougher3.partition_coefficient[:, "aqueous", "organic", "La"] = (100 - 0.5) / 100
    m.fs.solex_rougher3.partition_coefficient[:, "aqueous", "organic", "Ce"] = (100 - 0.5) / 100
    m.fs.solex_rougher3.partition_coefficient[:, "aqueous", "organic", "Pr"] = (100 - 0.5) / 100
    m.fs.solex_rougher3.partition_coefficient[:, "aqueous", "organic", "Nd"] = (100 - 0.5) / 100
    m.fs.solex_rougher3.partition_coefficient[:, "aqueous", "organic", "Sm"] = (100 - 0.5) / 100
    m.fs.solex_rougher3.partition_coefficient[:, "aqueous", "organic", "Gd"] = (100 - 0.5) / 100
    m.fs.solex_rougher3.partition_coefficient[:, "aqueous", "organic", "Dy"] = (100 - 0.5) / 100

    m.fs.sep = Separator(
        property_package=m.fs.prop_o,
        outlet_list=["recycle", "purge"],
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

    m.fs.sc_circuit_purge = Product(property_package=m.fs.prop_o)

    m.fs.solex_cleaner1 = SolventExtraction(
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
        aqueous_to_organic=True,
    )
    m.fs.solex_cleaner1.partition_coefficient[:, "aqueous", "organic", "Al"] = 3.6 / 100
    m.fs.solex_cleaner1.partition_coefficient[:, "aqueous", "organic", "Ca"] = 3.7 / 100
    m.fs.solex_cleaner1.partition_coefficient[:, "aqueous", "organic", "Fe"] = 2.1 / 100
    m.fs.solex_cleaner1.partition_coefficient[:, "aqueous", "organic", "Sc"] = 99.9 / 100
    m.fs.solex_cleaner1.partition_coefficient[:, "aqueous", "organic", "Y"] = 99.9 / 100
    m.fs.solex_cleaner1.partition_coefficient[:, "aqueous", "organic", "La"] = 75.2 / 100
    m.fs.solex_cleaner1.partition_coefficient[:, "aqueous", "organic", "Ce"] = 95.7 / 100
    m.fs.solex_cleaner1.partition_coefficient[:, "aqueous", "organic", "Pr"] = 96.5 / 100
    m.fs.solex_cleaner1.partition_coefficient[:, "aqueous", "organic", "Nd"] = 99.2 / 100
    m.fs.solex_cleaner1.partition_coefficient[:, "aqueous", "organic", "Sm"] = 99.9 / 100
    m.fs.solex_cleaner1.partition_coefficient[:, "aqueous", "organic", "Gd"] = 98.6 / 100
    m.fs.solex_cleaner1.partition_coefficient[:, "aqueous", "organic", "Dy"] = 99.9 / 100

    m.fs.solex_cleaner2 = SolventExtraction(
        number_of_finite_elements=3,
        dynamic=False,
        aqueous_stream={
            "property_package": m.fs.leach_soln,
            "flow_direction": FlowDirection.backward,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        organic_stream={
            "property_package": m.fs.prop_o,
            "flow_direction": FlowDirection.forward,
            "has_energy_balance": False,
            "has_pressure_balance": False,
        },
        aqueous_to_organic=False,
    )
    m.fs.solex_cleaner2.partition_coefficient[:, "aqueous", "organic", "Al"] = (100 - 0.5) / 100
    m.fs.solex_cleaner2.partition_coefficient[:, "aqueous", "organic", "Ca"] = (100 - 0.5) / 100
    m.fs.solex_cleaner2.partition_coefficient[:, "aqueous", "organic", "Fe"] = (100 - 5) / 100
    m.fs.solex_cleaner2.partition_coefficient[:, "aqueous", "organic", "Sc"] = (100 - 98.5) / 100
    m.fs.solex_cleaner2.partition_coefficient[:, "aqueous", "organic", "Y"] = (100 - 0.5) / 100
    m.fs.solex_cleaner2.partition_coefficient[:, "aqueous", "organic", "La"] = (100 - 0.5) / 100
    m.fs.solex_cleaner2.partition_coefficient[:, "aqueous", "organic", "Ce"] = (100 - 0.5) / 100
    m.fs.solex_cleaner2.partition_coefficient[:, "aqueous", "organic", "Pr"] = (100 - 0.5) / 100
    m.fs.solex_cleaner2.partition_coefficient[:, "aqueous", "organic", "Nd"] = (100 - 0.5) / 100
    m.fs.solex_cleaner2.partition_coefficient[:, "aqueous", "organic", "Sm"] = (100 - 0.5) / 100
    m.fs.solex_cleaner2.partition_coefficient[:, "aqueous", "organic", "Gd"] = (100 - 0.5) / 100
    m.fs.solex_cleaner2.partition_coefficient[:, "aqueous", "organic", "Dy"] = (100 - 0.5) / 100

    m.fs.organic_make_up2 = Feed(property_package=m.fs.prop_o)

    m.fs.mixer2 = Mixer(
        property_package=m.fs.prop_o,
        num_inlets=2,
        inlet_list=["make_up", "recycle"],
        material_balance_type=MaterialBalanceType.componentTotal,
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )
    m.fs.mixer3 = Mixer(
        property_package=m.fs.leach_soln,
        num_inlets=2,
        inlet_list=["leach", "cleaner"],
        material_balance_type=MaterialBalanceType.componentTotal,
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
    )

    m.fs.sep2 = Separator(
        property_package=m.fs.prop_o,
        outlet_list=["recycle", "purge"],
        split_basis=SplittingType.totalFlow,
        material_balance_type=MaterialBalanceType.componentTotal,
        momentum_balance_type=MomentumBalanceType.none,
        energy_split_basis=EnergySplittingType.none,
    )

    m.fs.acid_feed3 = Feed(property_package=m.fs.leach_soln)

    m.fs.precipitation_recycle = Product(property_package=m.fs.leach_soln)
    m.fs.cleaner_purge = Product(property_package=m.fs.prop_o)

    m.fs.aq_feed = Arc(
        source=m.fs.aqueous_feed.outlet, destination=m.fs.mixer3.leach
    )
    m.fs.mixed_aq_feed = Arc(
        source=m.fs.mixer3.outlet, destination=m.fs.solex_rougher1.mscontactor.aqueous_inlet
    )
    m.fs.org_feed = Arc(
        source=m.fs.organic_make_up.outlet, destination=m.fs.mixer.make_up
    )
    m.fs.mixed_org_feed = Arc(
        source=m.fs.mixer.outlet, destination=m.fs.solex_rougher1.mscontactor.organic_inlet
    )
    m.fs.s01 = Arc(
        source=m.fs.solex_rougher1.mscontactor.aqueous_outlet, destination=m.fs.leach_recycle.inlet
    )
    m.fs.s02 = Arc(
        source=m.fs.solex_rougher1.mscontactor.organic_outlet, destination=m.fs.solex_rougher2.mscontactor.organic_inlet
    )
    m.fs.s03 = Arc(
        source=m.fs.acid_feed1.outlet, destination=m.fs.solex_rougher2.mscontactor.aqueous_inlet
    )
    m.fs.s04 = Arc(
        source=m.fs.solex_rougher2.mscontactor.aqueous_outlet, destination=m.fs.leach_recycle2.inlet
    )
    m.fs.s05 = Arc(
        source=m.fs.solex_rougher2.mscontactor.organic_outlet, destination=m.fs.solex_rougher3.mscontactor.organic_inlet
    )
    m.fs.s06 = Arc(
        source=m.fs.acid_feed2.outlet, destination=m.fs.solex_rougher3.mscontactor.aqueous_inlet
    )
    m.fs.s07 = Arc(
        source=m.fs.solex_rougher3.mscontactor.organic_outlet, destination=m.fs.sep.inlet
    )
    m.fs.s08 = Arc(
        source=m.fs.sep.purge, destination=m.fs.sc_circuit_purge.inlet
    )
    m.fs.s09 = Arc(
        source=m.fs.sep.recycle, destination=m.fs.mixer.recycle
    )
    m.fs.s10 = Arc(
        source=m.fs.solex_rougher3.mscontactor.aqueous_outlet, destination=m.fs.solex_cleaner1.mscontactor.aqueous_inlet
    )
    m.fs.s11 = Arc(
        source=m.fs.organic_make_up2.outlet, destination=m.fs.mixer2.make_up
    )
    m.fs.s12 = Arc(
        source=m.fs.mixer2.outlet, destination=m.fs.solex_cleaner1.mscontactor.organic_inlet
    )
    m.fs.s13 = Arc(
        source=m.fs.solex_cleaner1.mscontactor.aqueous_outlet, destination=m.fs.mixer3.cleaner
    )
    m.fs.s14 = Arc(
        source=m.fs.acid_feed3.outlet, destination=m.fs.solex_cleaner2.mscontactor.aqueous_inlet
    )
    m.fs.s15 = Arc(
        source=m.fs.solex_cleaner1.mscontactor.organic_outlet, destination=m.fs.solex_cleaner2.mscontactor.organic_inlet
    )
    m.fs.s16 = Arc(
        source=m.fs.solex_cleaner2.mscontactor.organic_outlet, destination=m.fs.sep2.inlet
    )
    m.fs.s17 = Arc(
        source=m.fs.sep2.purge, destination=m.fs.cleaner_purge.inlet
    )
    m.fs.s18 = Arc(
        source=m.fs.sep2.recycle, destination=m.fs.mixer2.recycle
    )
    m.fs.s19 = Arc(
        source=m.fs.solex_cleaner2.mscontactor.aqueous_outlet, destination=m.fs.precipitation_recycle.inlet
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    m.scaling_factor = Suffix(direction=Suffix.EXPORT)

    component_set1 = [
        "H2O",
        "H",
        "HSO4",
        "SO4",
        "Cl",
        "Sc",
        "Y",
        "La",
        "Ce",
        "Pr",
        "Nd",
        "Sm",
        "Gd",
        "Dy",
        "Al",
        "Ca",
        "Fe",
    ]

    component_set2 = [
        "Sc",
        "Y",
        "La",
        "Ce",
        "Pr",
        "Nd",
        "Sm",
        "Gd",
        "Dy",
        "Al",
        "Ca",
        "Fe",
    ]

    for component in component_set1:
        m.scaling_factor[m.fs.aqueous_feed.properties[0.0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.solex_rougher1.mscontactor.aqueous[0, 1].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.solex_rougher1.mscontactor.aqueous[0, 2].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.solex_rougher1.mscontactor.aqueous[0, 3].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.solex_rougher1.mscontactor.aqueous_inlet_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.leach_recycle.properties[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.acid_feed1.properties[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.solex_rougher2.mscontactor.aqueous[0, 1].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.solex_rougher2.mscontactor.aqueous_inlet_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.leach_recycle2.properties[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.acid_feed2.properties[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.solex_rougher3.mscontactor.aqueous[0, 1].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.solex_rougher3.mscontactor.aqueous[0, 2].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.solex_rougher3.mscontactor.aqueous_inlet_state[0].conc_mol_comp[component]] = 1e5

    for component in component_set2:
        m.scaling_factor[m.fs.organic_make_up.properties[0.0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.solex_rougher1.mscontactor.organic[0, 1].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.solex_rougher1.mscontactor.organic[0, 2].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.solex_rougher1.mscontactor.organic[0, 3].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.solex_rougher1.mscontactor.organic_inlet_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.solex_rougher2.mscontactor.organic[0, 1].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.solex_rougher2.mscontactor.organic_inlet_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.solex_rougher3.mscontactor.organic[0, 1].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.solex_rougher3.mscontactor.organic[0, 2].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.solex_rougher3.mscontactor.organic_inlet_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.sep.mixed_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.sep.recycle_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.sep.purge_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.mixer.make_up_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.mixer.recycle_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.mixer.mixed_state[0].conc_mol_comp[component]] = 1e5
        m.scaling_factor[m.fs.sc_circuit_purge.properties[0].conc_mol_comp[component]] = 1e5

    scaling = TransformationFactory("core.scale_model")
    m = scaling.create_using(m, rename=False)

    m.fs.aqueous_feed.flow_vol.fix(4.4)

    m.fs.aqueous_feed.conc_mass_comp[0, "H2O"].fix(1e-9)
    m.fs.aqueous_feed.conc_mass_comp[0, "H"].fix(1e-9)
    m.fs.aqueous_feed.conc_mass_comp[0, "SO4"].fix(1e-9)
    m.fs.aqueous_feed.conc_mass_comp[0, "HSO4"].fix(1e-9)
    m.fs.aqueous_feed.conc_mass_comp[0, "Cl"].fix(1e-9)
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

    # Note: This stream + m.fs.s09 = 62.01 L/hr
    m.fs.organic_make_up.flow_vol.fix(6.201)

    m.fs.organic_make_up.conc_mass_comp[0, "Al"].fix(7.54e-10)
    m.fs.organic_make_up.conc_mass_comp[0, "Ca"].fix(4.955e-9)
    m.fs.organic_make_up.conc_mass_comp[0, "Fe"].fix(1.491e-7)
    m.fs.organic_make_up.conc_mass_comp[0, "Sc"].fix(321.34)
    m.fs.organic_make_up.conc_mass_comp[0, "Y"].fix(5.67e-6)
    m.fs.organic_make_up.conc_mass_comp[0, "La"].fix(1.78e-05)
    m.fs.organic_make_up.conc_mass_comp[0, "Ce"].fix(4.019e-5)
    m.fs.organic_make_up.conc_mass_comp[0, "Pr"].fix(6.73e-6)
    m.fs.organic_make_up.conc_mass_comp[0, "Nd"].fix(1.82e-5)
    m.fs.organic_make_up.conc_mass_comp[0, "Sm"].fix(3.285e-6)
    m.fs.organic_make_up.conc_mass_comp[0, "Gd"].fix(1.55e-6)
    m.fs.organic_make_up.conc_mass_comp[0, "Dy"].fix(9e-7)

    eps = 1e-7

    m.fs.acid_feed1.flow_vol.fix(0.09)
    m.fs.acid_feed1.conc_mass_comp[0, "H2O"].fix(1000000)
    m.fs.acid_feed1.conc_mass_comp[0, "H"].fix(10.36)
    m.fs.acid_feed1.conc_mass_comp[0, "SO4"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "HSO4"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "Cl"].fix(359.64)
    m.fs.acid_feed1.conc_mass_comp[0, "Al"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "Ca"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "Fe"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "Sc"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "Y"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "La"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "Ce"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "Pr"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "Nd"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "Sm"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "Gd"].fix(eps)
    m.fs.acid_feed1.conc_mass_comp[0, "Dy"].fix(eps)

    m.fs.acid_feed2.flow_vol.fix(0.09)
    m.fs.acid_feed2.conc_mass_comp[0, "H2O"].fix(1000000)
    m.fs.acid_feed2.conc_mass_comp[0, "H"].fix(10.36 * 4)  # Arbitrarily choose 4x the dilute solution
    m.fs.acid_feed2.conc_mass_comp[0, "SO4"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "HSO4"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "Cl"].fix(359.64 * 4)
    m.fs.acid_feed2.conc_mass_comp[0, "Al"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "Ca"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "Fe"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "Sc"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "Y"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "La"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "Ce"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "Pr"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "Nd"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "Sm"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "Gd"].fix(eps)
    m.fs.acid_feed2.conc_mass_comp[0, "Dy"].fix(eps)

    m.fs.sep.split_fraction[:, "recycle"].fix(0.9)

    # Note: This stream + m.fs.s18 = 62.01 L/hr
    m.fs.organic_make_up2.flow_vol.fix(6.201)

    m.fs.organic_make_up2.conc_mass_comp[0, "Al"].fix(7.54e-10)
    m.fs.organic_make_up2.conc_mass_comp[0, "Ca"].fix(4.955e-9)
    m.fs.organic_make_up2.conc_mass_comp[0, "Fe"].fix(1.491e-7)
    m.fs.organic_make_up2.conc_mass_comp[0, "Sc"].fix(321.34)
    m.fs.organic_make_up2.conc_mass_comp[0, "Y"].fix(5.67e-6)
    m.fs.organic_make_up2.conc_mass_comp[0, "La"].fix(1.78e-05)
    m.fs.organic_make_up2.conc_mass_comp[0, "Ce"].fix(4.019e-5)
    m.fs.organic_make_up2.conc_mass_comp[0, "Pr"].fix(6.73e-6)
    m.fs.organic_make_up2.conc_mass_comp[0, "Nd"].fix(1.82e-5)
    m.fs.organic_make_up2.conc_mass_comp[0, "Sm"].fix(3.285e-6)
    m.fs.organic_make_up2.conc_mass_comp[0, "Gd"].fix(1.55e-6)
    m.fs.organic_make_up2.conc_mass_comp[0, "Dy"].fix(9e-7)

    m.fs.acid_feed3.flow_vol.fix(0.09)
    m.fs.acid_feed3.conc_mass_comp[0, "H2O"].fix(1000000)
    m.fs.acid_feed3.conc_mass_comp[0, "H"].fix(10.36 * 4)  # Arbitrarily choose 4x the dilute solution
    m.fs.acid_feed3.conc_mass_comp[0, "SO4"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "HSO4"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "Cl"].fix(359.64 * 4)
    m.fs.acid_feed3.conc_mass_comp[0, "Al"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "Ca"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "Fe"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "Sc"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "Y"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "La"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "Ce"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "Pr"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "Nd"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "Sm"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "Gd"].fix(eps)
    m.fs.acid_feed3.conc_mass_comp[0, "Dy"].fix(eps)

    m.fs.sep2.split_fraction[:, "recycle"].fix(0.9)

    seq = SequentialDecomposition()
    seq.options.tear_method = "Direct"
    seq.options.iterLim = 1
    seq.options.tear_set = [m.fs.mixed_org_feed, m.fs.mixed_aq_feed, m.fs.s12]

    G = seq.create_graph(m)
    order = seq.calculation_order(G)
    print("Initialization Order")
    for o in order:
        print(o[0].name)

    tear_guesses1 = {
        "flow_vol": {0: 64},
        "conc_mass_comp": {
            (0, "Al"): 1e-9,
            (0, "Ca"): 1e-9,
            (0, "Ce"): 1e-5,
            (0, "Dy"): 1e-7,
            (0, "Fe"): 1e-7,
            (0, "Gd"): 1e-6,
            (0, "La"): 1e-5,
            (0, "Nd"): 1e-5,
            (0, "Pr"): 1e-6,
            (0, "Sc"): 321.34,
            (0, "Sm"): 1e-6,
            (0, "Y"): 1e-6,
        },
    }
    tear_guesses2 = {
        "flow_vol": {0: 40},
        "conc_mass_comp": {
            (0, "Al"): 820,
            (0, "Ca"): 5230,
            (0, "Ce"): 4516,
            (0, "Dy"): 101,
            (0, "Fe"): 270,
            (0, "Gd"): 174,
            (0, "H"): 2.274,
            (0, "H2O"): 1000000,
            (0, "HSO4"): 881.19,
            (0, "Cl"): 881.19,
            (0, "La"): 2032,
            (0, "Nd"): 2047,
            (0, "Pr"): 756,
            (0, "SO4"): 3924.07,
            (0, "Sc"): 209.31,
            (0, "Sm"): 369,
            (0, "Y"): 637,
        },
    }
    tear_guesses3 = {
        "flow_vol": {0: 64},
        "conc_mass_comp": {
            (0, "Al"): 1e-9,
            (0, "Ca"): 1e-9,
            (0, "Ce"): 1e-5,
            (0, "Dy"): 1e-7,
            (0, "Fe"): 1e-7,
            (0, "Gd"): 1e-6,
            (0, "La"): 1e-5,
            (0, "Nd"): 1e-5,
            (0, "Pr"): 1e-6,
            (0, "Sc"): 321.34,
            (0, "Sm"): 1e-6,
            (0, "Y"): 1e-6,
        },
    }

    # Pass the tear_guess to the SD tool
    seq.set_guesses_for(m.fs.solex_rougher1.mscontactor.organic_inlet, tear_guesses1)
    seq.set_guesses_for(m.fs.solex_rougher1.mscontactor.aqueous_inlet, tear_guesses2)
    seq.set_guesses_for(m.fs.solex_cleaner1.mscontactor.organic_inlet, tear_guesses3)

    def function(stream):
        initializer_feed = FeedInitializer()
        initializer = BlockTriangularizationInitializer()

        propagate_state(m.fs.aq_feed)
        propagate_state(m.fs.org_feed)

        if stream == m.fs.organic_make_up:
            initializer_feed.initialize(m.fs.organic_make_up)
        elif stream == m.fs.organic_make_up2:
            initializer_feed.initialize(m.fs.organic_make_up2)
        elif stream == m.fs.acid_feed1:
            initializer_feed.initialize(m.fs.acid_feed1)
        elif stream == m.fs.acid_feed2:
            initializer_feed.initialize(m.fs.acid_feed2)
        elif stream == m.fs.acid_feed3:
            initializer_feed.initialize(m.fs.acid_feed3)
        elif stream == m.fs.solex_rougher1.mscontactor:
            print(f"Initializing {stream}")
            initializer.initialize(m.fs.solex_rougher1)
        elif stream == m.fs.solex_rougher2.mscontactor:
            print(f"Initializing {stream}")
            try:
                initializer.initialize(m.fs.solex_rougher2)
            except:
                # Fix feed states
                m.fs.solex_rougher2.mscontactor.organic_inlet_state[0].flow_vol.fix()
                m.fs.solex_rougher2.mscontactor.aqueous_inlet_state[0].flow_vol.fix()
                m.fs.solex_rougher2.mscontactor.organic_inlet_state[0].conc_mass_comp.fix()
                m.fs.solex_rougher2.mscontactor.aqueous_inlet_state[0].conc_mass_comp.fix()
                # Re-solve leach unit
                solver = SolverFactory("ipopt")
                solver.solve(m.fs.solex_rougher2, tee=True)
                # Unfix feed states
                m.fs.solex_rougher2.mscontactor.organic_inlet_state[0].flow_vol.unfix()
                m.fs.solex_rougher2.mscontactor.aqueous_inlet_state[0].flow_vol.unfix()
                m.fs.solex_rougher2.mscontactor.organic_inlet_state[0].conc_mass_comp.unfix()
                m.fs.solex_rougher2.mscontactor.aqueous_inlet_state[0].conc_mass_comp.unfix()
        elif stream == m.fs.solex_rougher3.mscontactor:
            print(f"Initializing {stream}")
            try:
                initializer.initialize(m.fs.solex_rougher3)
            except:
                # Fix feed states
                m.fs.solex_rougher3.mscontactor.organic_inlet_state[0].flow_vol.fix()
                m.fs.solex_rougher3.mscontactor.aqueous_inlet_state[0].flow_vol.fix()
                m.fs.solex_rougher3.mscontactor.organic_inlet_state[0].conc_mass_comp.fix()
                m.fs.solex_rougher3.mscontactor.aqueous_inlet_state[0].conc_mass_comp.fix()
                # Re-solve leach unit
                solver = SolverFactory("ipopt")
                solver.solve(m.fs.solex_rougher3, tee=True)
                # Unfix feed states
                m.fs.solex_rougher3.mscontactor.organic_inlet_state[0].flow_vol.unfix()
                m.fs.solex_rougher3.mscontactor.aqueous_inlet_state[0].flow_vol.unfix()
                m.fs.solex_rougher3.mscontactor.organic_inlet_state[0].conc_mass_comp.unfix()
                m.fs.solex_rougher3.mscontactor.aqueous_inlet_state[0].conc_mass_comp.unfix()
        elif stream == m.fs.solex_cleaner1.mscontactor:
            print(f"Initializing {stream}")
            try:
                initializer.initialize(m.fs.solex_cleaner1)
            except:
                # Fix feed states
                m.fs.solex_cleaner1.mscontactor.organic_inlet_state[0].flow_vol.fix()
                m.fs.solex_cleaner1.mscontactor.aqueous_inlet_state[0].flow_vol.fix()
                m.fs.solex_cleaner1.mscontactor.organic_inlet_state[0].conc_mass_comp.fix()
                m.fs.solex_cleaner1.mscontactor.aqueous_inlet_state[0].conc_mass_comp.fix()
                # Re-solve leach unit
                solver = SolverFactory("ipopt")
                solver.solve(m.fs.solex_cleaner1, tee=True)
                # Unfix feed states
                m.fs.solex_cleaner1.mscontactor.organic_inlet_state[0].flow_vol.unfix()
                m.fs.solex_cleaner1.mscontactor.aqueous_inlet_state[0].flow_vol.unfix()
                m.fs.solex_cleaner1.mscontactor.organic_inlet_state[0].conc_mass_comp.unfix()
                m.fs.solex_cleaner1.mscontactor.aqueous_inlet_state[0].conc_mass_comp.unfix()
        elif stream == m.fs.solex_cleaner2.mscontactor:
            print(f"Initializing {stream}")
            try:
                initializer.initialize(m.fs.solex_cleaner2)
            except:
                # Fix feed states
                m.fs.solex_cleaner2.mscontactor.organic_inlet_state[0].flow_vol.fix()
                m.fs.solex_cleaner2.mscontactor.aqueous_inlet_state[0].flow_vol.fix()
                m.fs.solex_cleaner2.mscontactor.organic_inlet_state[0].conc_mass_comp.fix()
                m.fs.solex_cleaner2.mscontactor.aqueous_inlet_state[0].conc_mass_comp.fix()
                # Re-solve leach unit
                solver = SolverFactory("ipopt")
                solver.solve(m.fs.solex_cleaner2, tee=True)
                # Unfix feed states
                m.fs.solex_cleaner2.mscontactor.organic_inlet_state[0].flow_vol.unfix()
                m.fs.solex_cleaner2.mscontactor.aqueous_inlet_state[0].flow_vol.unfix()
                m.fs.solex_cleaner2.mscontactor.organic_inlet_state[0].conc_mass_comp.unfix()
                m.fs.solex_cleaner2.mscontactor.aqueous_inlet_state[0].conc_mass_comp.unfix()

        else:
            print(f"Initializing {stream}")
            initializer.initialize(stream)

    seq.run(m, function)

    solver = SolverFactory("ipopt")
    # solver.options["bound_push"] = 1e-8
    # solver.options["mu_init"] = 1e-8
    results = solver.solve(m, tee=True)

    m.fs.solex_rougher1.display()
    m.fs.solex_rougher2.display()
    m.fs.solex_rougher3.display()

    m.fs.solex_cleaner1.display()
    m.fs.solex_cleaner2.display()
    return results


if __name__ == "__main__":
    results = main()
