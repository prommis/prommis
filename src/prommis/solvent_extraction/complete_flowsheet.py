from pyomo.environ import ConcreteModel, SolverFactory, units, TransformationFactory
from pyomo.network import Arc

import numpy as np

from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
    MaterialBalanceType,
    MomentumBalanceType,
)
from idaes.core.initialization.block_triangularization import (
    BlockTriangularizationInitializer,
)
from idaes.models.unit_models.mixer import (
    Mixer,
    MixingType,
    MomentumMixingType,
    MixerInitializer,
)
from idaes.models.unit_models.feed import Feed, FeedInitializer
from idaes.models.unit_models.product import Product, ProductInitializer
from idaes.models.unit_models.separator import (
    EnergySplittingType,
    Separator,
    SplittingType,
    SeparatorInitializer,
)
from idaes.core.util.model_statistics import degrees_of_freedom as dof

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.solvent_extraction import (
    SolventExtraction,
    SolventExtractionInitializer,
)


# Model development and flowsheet creation

m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)


# Property packages

m.fs.leach_soln = LeachSolutionParameters()
m.fs.prop_o = REESolExOgParameters()


# Rougher circuit

m.fs.rougher_org_make_up = Feed(property_package=m.fs.prop_o)

m.fs.solex_rougher_load = SolventExtraction(
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

stage_number = np.arange(1, 4)

for s in stage_number:
    if s == 1:
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "Al"] = (
            5.2 / 100
        )
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "Ca"] = (
            3 / 100
        )
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "Fe"] = (
            24.7 / 100
        )
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "Sc"] = (
            99.1 / 100
        )
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "Y"] = (
            99.9 / 100
        )
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "La"] = (
            32.4 / 100
        )
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "Ce"] = (
            58.2 / 100
        )
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "Pr"] = (
            58.2 / 100
        )
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "Nd"] = (
            87.6 / 100
        )
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "Sm"] = (
            99.9 / 100
        )
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "Gd"] = (
            69.8 / 100
        )
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "Dy"] = (
            96.6 / 100
        )
    else:
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "Al"] = (
            4.9 / 100
        )
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "Ca"] = (
            12.3 / 100
        )
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "Fe"] = (
            6.4 / 100
        )
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "Sc"] = (
            16.7 / 100
        )
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "Y"] = (
            99.9 / 100
        )
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "La"] = (
            23.2 / 100
        )
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "Ce"] = (
            24.9 / 100
        )
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "Pr"] = (
            15.1 / 100
        )
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "Nd"] = (
            99.9 / 100
        )
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "Sm"] = (
            99.9 / 100
        )
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "Gd"] = (
            7.6 / 100
        )
        m.fs.solex_rougher_load.partition_coefficient[s, "aqueous", "organic", "Dy"] = (
            5 / 100
        )


m.fs.solex_rougher_scrub = SolventExtraction(
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

m.fs.acid_feed_to_rough_scrub = Feed(property_package=m.fs.leach_soln)

m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "Al"] = (
    100 - 0.12
) / 100
m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "Ca"] = (
    100 - 0.55
) / 100
m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "Fe"] = (
    100 - 0.007
) / 100
m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "Sc"] = (
    100 - 99.9
) / 100
m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "Y"] = (
    100 - 99.9
) / 100
m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "La"] = (
    100 - 99.8
) / 100
m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "Ce"] = (
    100 - 99.9
) / 100
m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "Pr"] = (
    100 - 99.9
) / 100
m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "Nd"] = (
    100 - 99.9
) / 100
m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "Sm"] = (
    100 - 99.9
) / 100
m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "Gd"] = (
    100 - 99.9
) / 100
m.fs.solex_rougher_scrub.partition_coefficient[1, "aqueous", "organic", "Dy"] = (
    100 - 99.9
) / 100

m.fs.acid_feed_to_rough_strip = Feed(property_package=m.fs.leach_soln)

m.fs.solex_rougher_strip = SolventExtraction(
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

m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "Al"] = (
    100 - 0.5
) / 100
m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "Ca"] = (
    100 - 0.5
) / 100
m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "Fe"] = (
    100 - 0.5
) / 100
m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "Sc"] = (
    100 - 98.5
) / 100
m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "Y"] = (
    100 - 0.5
) / 100
m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "La"] = (
    100 - 0.5
) / 100
m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "Ce"] = (
    100 - 0.5
) / 100
m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "Pr"] = (
    100 - 0.5
) / 100
m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "Nd"] = (
    100 - 0.5
) / 100
m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "Sm"] = (
    100 - 0.5
) / 100
m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "Gd"] = (
    100 - 0.5
) / 100
m.fs.solex_rougher_strip.partition_coefficient[:, "aqueous", "organic", "Dy"] = (
    100 - 0.5
) / 100

m.fs.rougher_sep = Separator(
    property_package=m.fs.prop_o,
    outlet_list=["recycle", "purge"],
    split_basis=SplittingType.totalFlow,
    material_balance_type=MaterialBalanceType.componentTotal,
    momentum_balance_type=MomentumBalanceType.none,
    energy_split_basis=EnergySplittingType.none,
)

m.fs.rougher_mixer = Mixer(
    property_package=m.fs.prop_o,
    num_inlets=2,
    inlet_list=["make_up", "recycle"],
    material_balance_type=MaterialBalanceType.componentTotal,
    energy_mixing_type=MixingType.none,
    momentum_mixing_type=MomentumMixingType.none,
)

m.fs.sc_circuit_purge = Product(property_package=m.fs.prop_o)


# Cleaner Circuit

m.fs.solex_cleaner_load = SolventExtraction(
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

m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "Al"] = 3.6 / 100
m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "Ca"] = 3.7 / 100
m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "Fe"] = 2.1 / 100
m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "Sc"] = (
    99.9 / 100
)
m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "Y"] = 99.9 / 100
m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "La"] = (
    75.2 / 100
)
m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "Ce"] = (
    95.7 / 100
)
m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "Pr"] = (
    96.5 / 100
)
m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "Nd"] = (
    99.2 / 100
)
m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "Sm"] = (
    99.9 / 100
)
m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "Gd"] = (
    98.6 / 100
)
m.fs.solex_cleaner_load.partition_coefficient[:, "aqueous", "organic", "Dy"] = (
    99.9 / 100
)

m.fs.solex_cleaner_strip = SolventExtraction(
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

m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "Al"] = (
    100 - 0.5
) / 100
m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "Ca"] = (
    100 - 0.5
) / 100
m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "Fe"] = (
    100 - 5
) / 100
m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "Sc"] = (
    100 - 98.5
) / 100
m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "Y"] = (
    100 - 0.5
) / 100
m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "La"] = (
    100 - 0.5
) / 100
m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "Ce"] = (
    100 - 0.5
) / 100
m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "Pr"] = (
    100 - 0.5
) / 100
m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "Nd"] = (
    100 - 0.5
) / 100
m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "Sm"] = (
    100 - 0.5
) / 100
m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "Gd"] = (
    100 - 0.5
) / 100
m.fs.solex_cleaner_strip.partition_coefficient[:, "aqueous", "organic", "Dy"] = (
    100 - 0.5
) / 100

m.fs.cleaner_org_make_up = Feed(property_package=m.fs.prop_o)

m.fs.cleaner_sep = Separator(
    property_package=m.fs.prop_o,
    outlet_list=["recycle", "purge"],
    split_basis=SplittingType.totalFlow,
    material_balance_type=MaterialBalanceType.componentTotal,
    momentum_balance_type=MomentumBalanceType.none,
    energy_split_basis=EnergySplittingType.none,
)

m.fs.cleaner_mixer = Mixer(
    property_package=m.fs.prop_o,
    num_inlets=2,
    inlet_list=["make_up", "recycle"],
    material_balance_type=MaterialBalanceType.componentTotal,
    energy_mixing_type=MixingType.none,
    momentum_mixing_type=MomentumMixingType.none,
)

m.fs.cleaner_purge = Product(property_package=m.fs.prop_o)
m.fs.acid_feed_to_clean_strip = Feed(property_package=m.fs.leach_soln)

m.fs.leach_sx_mixer = Mixer(
    property_package=m.fs.leach_soln,
    num_inlets=2,
    inlet_list=["leach", "cleaner"],
    material_balance_type=MaterialBalanceType.componentTotal,
    energy_mixing_type=MixingType.none,
    momentum_mixing_type=MomentumMixingType.none,
)

m.fs.acid_feed_to_leach_mixer = Feed(property_package=m.fs.leach_soln)

m.fs.precip_sx_mixer = Mixer(
    property_package=m.fs.leach_soln,
    num_inlets=2,
    inlet_list=["precip", "rougher"],
    material_balance_type=MaterialBalanceType.componentTotal,
    energy_mixing_type=MixingType.none,
    momentum_mixing_type=MomentumMixingType.none,
)

# Define outlets from the flowsheet
m.fs.precip_purge = Feed(property_package=m.fs.leach_soln)


# Connections

m.fs.rough_org_load_to_scrub = Arc(
    source=m.fs.solex_rougher_load.mscontactor.organic_outlet,
    destination=m.fs.solex_rougher_scrub.mscontactor.organic_inlet,
)
m.fs.rough_org_scrub_to_strip = Arc(
    source=m.fs.solex_rougher_scrub.mscontactor.organic_outlet,
    destination=m.fs.solex_rougher_strip.mscontactor.organic_inlet,
)
m.fs.rough_org_strip_to_sep = Arc(
    source=m.fs.solex_rougher_strip.mscontactor.organic_outlet,
    destination=m.fs.rougher_sep.inlet,
)
m.fs.rough_org_sep_to_purge = Arc(
    source=m.fs.rougher_sep.purge, destination=m.fs.sc_circuit_purge.inlet
)
m.fs.rough_org_sep_to_mixer = Arc(
    source=m.fs.rougher_sep.recycle, destination=m.fs.rougher_mixer.recycle
)
m.fs.rough_org_makeup_to_mixer = Arc(
    source=m.fs.rougher_org_make_up.outlet, destination=m.fs.rougher_mixer.make_up
)
m.fs.rough_org_mixer_to_load = Arc(
    source=m.fs.rougher_mixer.outlet,
    destination=m.fs.solex_rougher_load.mscontactor.organic_inlet,
)

m.fs.rough_aq_feed_to_mixer = Arc(
    source=m.fs.acid_feed_to_leach_mixer.outlet, destination=m.fs.leach_sx_mixer.leach
)
m.fs.rough_aq_clean_to_mixer = Arc(
    source=m.fs.solex_cleaner_load.mscontactor.aqueous_outlet,
    destination=m.fs.leach_sx_mixer.cleaner,
)
m.fs.rough_aq_mixer_to_load = Arc(
    source=m.fs.leach_sx_mixer.outlet,
    destination=m.fs.solex_rougher_load.mscontactor.aqueous_inlet,
)

m.fs.rough_aq_feed_to_scrub = Arc(
    source=m.fs.acid_feed_to_rough_scrub.outlet,
    destination=m.fs.solex_rougher_scrub.mscontactor.aqueous_inlet,
)

m.fs.rough_aq_feed_to_strip = Arc(
    source=m.fs.acid_feed_to_rough_strip.outlet,
    destination=m.fs.solex_rougher_strip.mscontactor.aqueous_inlet,
)

m.fs.rough_aq_strip_to_mixer = Arc(
    source=m.fs.solex_rougher_strip.mscontactor.aqueous_outlet,
    destination=m.fs.precip_sx_mixer.rougher,
)
m.fs.rough_precip_feed_to_mixer = Arc(
    source=m.fs.precip_purge.outlet,
    destination=m.fs.precip_sx_mixer.precip,
)
m.fs.rough_mixer_to_clean_load = Arc(
    source=m.fs.precip_sx_mixer.outlet,
    destination=m.fs.solex_cleaner_load.mscontactor.aqueous_inlet,
)
m.fs.clean_org_load_to_strip = Arc(
    source=m.fs.solex_cleaner_load.mscontactor.organic_outlet,
    destination=m.fs.solex_cleaner_strip.mscontactor.organic_inlet,
)
m.fs.clean_org_strip_to_sep = Arc(
    source=m.fs.solex_cleaner_strip.mscontactor.organic_outlet,
    destination=m.fs.cleaner_sep.inlet,
)
m.fs.clean_org_sep_to_purge = Arc(
    source=m.fs.cleaner_sep.purge, destination=m.fs.cleaner_purge.inlet
)
m.fs.clean_org_sep_to_mixer = Arc(
    source=m.fs.cleaner_sep.recycle, destination=m.fs.cleaner_mixer.recycle
)
m.fs.clean_org_makeup_to_mixer = Arc(
    source=m.fs.cleaner_org_make_up.outlet, destination=m.fs.cleaner_mixer.make_up
)
m.fs.clean_org_mixer_to_load = Arc(
    source=m.fs.cleaner_mixer.outlet,
    destination=m.fs.solex_cleaner_load.mscontactor.organic_inlet,
)

m.fs.clean_aq_feed_to_strip = Arc(
    source=m.fs.acid_feed_to_clean_strip.outlet,
    destination=m.fs.solex_cleaner_strip.mscontactor.aqueous_inlet,
)

TransformationFactory("network.expand_arcs").apply_to(m)

# Fix the inlets

eps = 1e-7 * units.mg / units.L

m.fs.rougher_org_make_up.flow_vol.fix(6.201)
m.fs.rougher_org_make_up.conc_mass_comp[0, "Al"].fix(eps)
m.fs.rougher_org_make_up.conc_mass_comp[0, "Ca"].fix(eps)
m.fs.rougher_org_make_up.conc_mass_comp[0, "Fe"].fix(eps)
m.fs.rougher_org_make_up.conc_mass_comp[0, "Sc"].fix(eps)
m.fs.rougher_org_make_up.conc_mass_comp[0, "Y"].fix(eps)
m.fs.rougher_org_make_up.conc_mass_comp[0, "La"].fix(eps)
m.fs.rougher_org_make_up.conc_mass_comp[0, "Ce"].fix(eps)
m.fs.rougher_org_make_up.conc_mass_comp[0, "Pr"].fix(eps)
m.fs.rougher_org_make_up.conc_mass_comp[0, "Nd"].fix(eps)
m.fs.rougher_org_make_up.conc_mass_comp[0, "Sm"].fix(eps)
m.fs.rougher_org_make_up.conc_mass_comp[0, "Gd"].fix(eps)
m.fs.rougher_org_make_up.conc_mass_comp[0, "Dy"].fix(eps)

m.fs.acid_feed_to_rough_scrub.flow_vol.fix(0.09)
m.fs.acid_feed_to_rough_scrub.conc_mass_comp[0, "H2O"].fix(1000000)
m.fs.acid_feed_to_rough_scrub.conc_mass_comp[0, "H"].fix(10.36)
m.fs.acid_feed_to_rough_scrub.conc_mass_comp[0, "SO4"].fix(eps)
m.fs.acid_feed_to_rough_scrub.conc_mass_comp[0, "HSO4"].fix(eps)
m.fs.acid_feed_to_rough_scrub.conc_mass_comp[0, "Al"].fix(eps)
m.fs.acid_feed_to_rough_scrub.conc_mass_comp[0, "Ca"].fix(eps)
m.fs.acid_feed_to_rough_scrub.conc_mass_comp[0, "Cl"].fix(eps)
m.fs.acid_feed_to_rough_scrub.conc_mass_comp[0, "Fe"].fix(eps)
m.fs.acid_feed_to_rough_scrub.conc_mass_comp[0, "Sc"].fix(eps)
m.fs.acid_feed_to_rough_scrub.conc_mass_comp[0, "Y"].fix(eps)
m.fs.acid_feed_to_rough_scrub.conc_mass_comp[0, "La"].fix(eps)
m.fs.acid_feed_to_rough_scrub.conc_mass_comp[0, "Ce"].fix(eps)
m.fs.acid_feed_to_rough_scrub.conc_mass_comp[0, "Pr"].fix(eps)
m.fs.acid_feed_to_rough_scrub.conc_mass_comp[0, "Nd"].fix(eps)
m.fs.acid_feed_to_rough_scrub.conc_mass_comp[0, "Sm"].fix(eps)
m.fs.acid_feed_to_rough_scrub.conc_mass_comp[0, "Gd"].fix(eps)
m.fs.acid_feed_to_rough_scrub.conc_mass_comp[0, "Dy"].fix(eps)

m.fs.acid_feed_to_rough_strip.flow_vol.fix(0.09)
m.fs.acid_feed_to_rough_strip.conc_mass_comp[0, "H2O"].fix(1000000)
m.fs.acid_feed_to_rough_strip.conc_mass_comp[0, "H"].fix(10.36 * 4)
m.fs.acid_feed_to_rough_strip.conc_mass_comp[0, "SO4"].fix(eps)
m.fs.acid_feed_to_rough_strip.conc_mass_comp[0, "HSO4"].fix(eps)
m.fs.acid_feed_to_rough_strip.conc_mass_comp[0, "Al"].fix(eps)
m.fs.acid_feed_to_rough_strip.conc_mass_comp[0, "Ca"].fix(eps)
m.fs.acid_feed_to_rough_strip.conc_mass_comp[0, "Cl"].fix(eps)
m.fs.acid_feed_to_rough_strip.conc_mass_comp[0, "Fe"].fix(eps)
m.fs.acid_feed_to_rough_strip.conc_mass_comp[0, "Sc"].fix(eps)
m.fs.acid_feed_to_rough_strip.conc_mass_comp[0, "Y"].fix(eps)
m.fs.acid_feed_to_rough_strip.conc_mass_comp[0, "La"].fix(eps)
m.fs.acid_feed_to_rough_strip.conc_mass_comp[0, "Ce"].fix(eps)
m.fs.acid_feed_to_rough_strip.conc_mass_comp[0, "Pr"].fix(eps)
m.fs.acid_feed_to_rough_strip.conc_mass_comp[0, "Nd"].fix(eps)
m.fs.acid_feed_to_rough_strip.conc_mass_comp[0, "Sm"].fix(eps)
m.fs.acid_feed_to_rough_strip.conc_mass_comp[0, "Gd"].fix(eps)
m.fs.acid_feed_to_rough_strip.conc_mass_comp[0, "Dy"].fix(eps)

m.fs.acid_feed_to_clean_strip.flow_vol.fix(9)
m.fs.acid_feed_to_clean_strip.conc_mass_comp[0, "H2O"].fix(1000000)
m.fs.acid_feed_to_clean_strip.conc_mass_comp[0, "H"].fix(10.36 * 4)
m.fs.acid_feed_to_clean_strip.conc_mass_comp[0, "SO4"].fix(eps)
m.fs.acid_feed_to_clean_strip.conc_mass_comp[0, "HSO4"].fix(eps)
m.fs.acid_feed_to_clean_strip.conc_mass_comp[0, "Al"].fix(eps)
m.fs.acid_feed_to_clean_strip.conc_mass_comp[0, "Ca"].fix(eps)
m.fs.acid_feed_to_clean_strip.conc_mass_comp[0, "Cl"].fix(eps)
m.fs.acid_feed_to_clean_strip.conc_mass_comp[0, "Fe"].fix(eps)
m.fs.acid_feed_to_clean_strip.conc_mass_comp[0, "Sc"].fix(eps)
m.fs.acid_feed_to_clean_strip.conc_mass_comp[0, "Y"].fix(eps)
m.fs.acid_feed_to_clean_strip.conc_mass_comp[0, "La"].fix(eps)
m.fs.acid_feed_to_clean_strip.conc_mass_comp[0, "Ce"].fix(eps)
m.fs.acid_feed_to_clean_strip.conc_mass_comp[0, "Pr"].fix(eps)
m.fs.acid_feed_to_clean_strip.conc_mass_comp[0, "Nd"].fix(eps)
m.fs.acid_feed_to_clean_strip.conc_mass_comp[0, "Sm"].fix(eps)
m.fs.acid_feed_to_clean_strip.conc_mass_comp[0, "Gd"].fix(eps)
m.fs.acid_feed_to_clean_strip.conc_mass_comp[0, "Dy"].fix(eps)

m.fs.cleaner_org_make_up.flow_vol.fix(6.201)
m.fs.cleaner_org_make_up.conc_mass_comp[0, "Al"].fix(eps)
m.fs.cleaner_org_make_up.conc_mass_comp[0, "Ca"].fix(eps)
m.fs.cleaner_org_make_up.conc_mass_comp[0, "Fe"].fix(eps)
m.fs.cleaner_org_make_up.conc_mass_comp[0, "Sc"].fix(eps)
m.fs.cleaner_org_make_up.conc_mass_comp[0, "Y"].fix(eps)
m.fs.cleaner_org_make_up.conc_mass_comp[0, "La"].fix(eps)
m.fs.cleaner_org_make_up.conc_mass_comp[0, "Ce"].fix(eps)
m.fs.cleaner_org_make_up.conc_mass_comp[0, "Pr"].fix(eps)
m.fs.cleaner_org_make_up.conc_mass_comp[0, "Nd"].fix(eps)
m.fs.cleaner_org_make_up.conc_mass_comp[0, "Sm"].fix(eps)
m.fs.cleaner_org_make_up.conc_mass_comp[0, "Gd"].fix(eps)
m.fs.cleaner_org_make_up.conc_mass_comp[0, "Dy"].fix(eps)

m.fs.acid_feed_to_leach_mixer.flow_vol.fix(62.01 - 0.09 - 0.09)
m.fs.acid_feed_to_leach_mixer.conc_mass_comp[0, "H2O"].fix(1000000)
m.fs.acid_feed_to_leach_mixer.conc_mass_comp[0, "H"].fix(10.36)
m.fs.acid_feed_to_leach_mixer.conc_mass_comp[0, "SO4"].fix(3999.818)
m.fs.acid_feed_to_leach_mixer.conc_mass_comp[0, "HSO4"].fix(693.459)
m.fs.acid_feed_to_leach_mixer.conc_mass_comp[0, "Al"].fix(422.375)
m.fs.acid_feed_to_leach_mixer.conc_mass_comp[0, "Ca"].fix(109.542)
m.fs.acid_feed_to_leach_mixer.conc_mass_comp[0, "Cl"].fix(1e-7)
m.fs.acid_feed_to_leach_mixer.conc_mass_comp[0, "Fe"].fix(688.266)
m.fs.acid_feed_to_leach_mixer.conc_mass_comp[0, "Sc"].fix(0.032)
m.fs.acid_feed_to_leach_mixer.conc_mass_comp[0, "Y"].fix(0.124)
m.fs.acid_feed_to_leach_mixer.conc_mass_comp[0, "La"].fix(0.986)
m.fs.acid_feed_to_leach_mixer.conc_mass_comp[0, "Ce"].fix(2.277)
m.fs.acid_feed_to_leach_mixer.conc_mass_comp[0, "Pr"].fix(0.303)
m.fs.acid_feed_to_leach_mixer.conc_mass_comp[0, "Nd"].fix(0.946)
m.fs.acid_feed_to_leach_mixer.conc_mass_comp[0, "Sm"].fix(0.097)
m.fs.acid_feed_to_leach_mixer.conc_mass_comp[0, "Gd"].fix(0.2584)
m.fs.acid_feed_to_leach_mixer.conc_mass_comp[0, "Dy"].fix(0.047)

m.fs.precip_purge.flow_vol.fix(0.09)
m.fs.precip_purge.conc_mass_comp[0, "H2O"].fix(1000000)
m.fs.precip_purge.conc_mass_comp[0, "H"].fix(eps)
m.fs.precip_purge.conc_mass_comp[0, "SO4"].fix(eps)
m.fs.precip_purge.conc_mass_comp[0, "HSO4"].fix(eps)
m.fs.precip_purge.conc_mass_comp[0, "Al"].fix(eps)
m.fs.precip_purge.conc_mass_comp[0, "Ca"].fix(eps)
m.fs.precip_purge.conc_mass_comp[0, "Cl"].fix(eps)
m.fs.precip_purge.conc_mass_comp[0, "Fe"].fix(eps)
m.fs.precip_purge.conc_mass_comp[0, "Sc"].fix(eps)
m.fs.precip_purge.conc_mass_comp[0, "Y"].fix(eps)
m.fs.precip_purge.conc_mass_comp[0, "La"].fix(eps)
m.fs.precip_purge.conc_mass_comp[0, "Ce"].fix(eps)
m.fs.precip_purge.conc_mass_comp[0, "Pr"].fix(eps)
m.fs.precip_purge.conc_mass_comp[0, "Nd"].fix(eps)
m.fs.precip_purge.conc_mass_comp[0, "Sm"].fix(eps)
m.fs.precip_purge.conc_mass_comp[0, "Gd"].fix(eps)
m.fs.precip_purge.conc_mass_comp[0, "Dy"].fix(eps)

m.fs.rougher_sep.split_fraction[:, "recycle"].fix(0.9)
m.fs.cleaner_sep.split_fraction[:, "recycle"].fix(0.9)

print(dof(m))

initializer_feed = FeedInitializer()
feed_units = [
    m.fs.rougher_org_make_up,
    m.fs.acid_feed_to_rough_scrub,
    m.fs.acid_feed_to_rough_strip,
    m.fs.acid_feed_to_clean_strip,
    m.fs.cleaner_org_make_up,
    m.fs.acid_feed_to_leach_mixer,
    m.fs.precip_purge,
]

initializer_product = ProductInitializer()
product_units = [m.fs.sc_circuit_purge, m.fs.cleaner_purge]

initializer_mixer = MixerInitializer()
mixer_units = [
    m.fs.rougher_mixer,
    m.fs.cleaner_mixer,
    m.fs.leach_sx_mixer,
    m.fs.precip_sx_mixer,
]

initializer_sep = SeparatorInitializer()
sep_units = [m.fs.rougher_sep, m.fs.cleaner_sep]

initializer_sx = SolventExtractionInitializer()
sx_units = [
    m.fs.solex_rougher_load,
    m.fs.solex_rougher_scrub,
    m.fs.solex_rougher_strip,
    m.fs.solex_cleaner_load,
    m.fs.solex_cleaner_strip,
]

for feed in feed_units:
    initializer_feed.initialize(feed)

for product in product_units:
    initializer_product.initialize(product)

for mixer in mixer_units:
    initializer_mixer.initialize(mixer)

for sep in sep_units:
    initializer_sep.initialize(sep)

for sx in sx_units:
    initializer_sx.initialize(sx)

solver = SolverFactory("ipopt")
results = solver.solve(m, tee=True)
