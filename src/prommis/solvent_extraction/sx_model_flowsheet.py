from pyomo.environ import ConcreteModel, SolverFactory, TransformationFactory, Suffix
from pyomo.network import Arc, SequentialDecomposition

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
from idaes.models.unit_models.mixer import Mixer, MixingType, MomentumMixingType
from idaes.models.unit_models.product import Product, ProductInitializer
from idaes.models.unit_models.feed import Feed, FeedInitializer
from idaes.models.unit_models.mscontactor import MSContactor, MSContactorInitializer

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.solvent_extraction import SolventExtraction

from idaes.core.util.model_diagnostics import DiagnosticsToolbox
import idaes.core.util.scaling as iscale
from idaes.core.util.initialization import propagate_state

def main():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.prop_o = REESolExOgParameters()
    m.fs.leach_soln = LeachSolutionParameters()

    # TODO: May need to adjust feed such that mixed flow matches REESim
    m.fs.organic_make_up = Feed(property_package=m.fs.prop_o)

    # TODO: Double check configuration and partition coefficients
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
    )

    m.fs.leach_recycle = Product(property_package=m.fs.leach_soln)
    # TODO: Double check what this HCl stream should be doing
    m.fs.acid_feed1 = Feed(property_package=m.fs.leach_soln)
    # TODO: Double check configuration and partition coefficients
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
    )

    m.fs.leach_recycle2 = Product(property_package=m.fs.leach_soln)
    # TODO: Double check what this HCl stream should be doing
    m.fs.acid_feed2 = Feed(property_package=m.fs.leach_soln)
    # TODO: Double check configuration and partition coefficients
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
    )
    # TODO: Double check split fraction
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
    #
    m.fs.sc_circuit_purge = Product(property_package=m.fs.prop_o)
    m.fs.solex_cleaner_recycle = Product(property_package=m.fs.leach_soln)

    m.fs.feed = Arc(
        source=m.fs.organic_make_up.outlet, destination=m.fs.mixer.make_up
    )
    m.fs.mixed_feed = Arc(
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
        source=m.fs.solex_rougher3.mscontactor.aqueous_outlet, destination=m.fs.solex_cleaner_recycle.inlet
    )
    m.fs.s08 = Arc(
        source=m.fs.solex_rougher3.mscontactor.organic_outlet, destination=m.fs.sep.inlet
    )
    m.fs.s09 = Arc(
        source=m.fs.sep.purge, destination=m.fs.sc_circuit_purge.inlet
    )
    m.fs.s10 = Arc(
        source=m.fs.sep.recycle, destination=m.fs.mixer.recycle
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
        m.scaling_factor[m.fs.solex_cleaner_recycle.properties[0].conc_mol_comp[component]] = 1e5

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

    # TODO: verify where Arko got these values from
    m.fs.solex_rougher1.mscontactor.aqueous_inlet_state[0].flow_vol.fix(4.4)

    m.fs.solex_rougher1.mscontactor.aqueous_inlet_state[0].conc_mass_comp["H2O"].fix(1e-9)
    m.fs.solex_rougher1.mscontactor.aqueous_inlet_state[0].conc_mass_comp["H"].fix(1e-9)
    m.fs.solex_rougher1.mscontactor.aqueous_inlet_state[0].conc_mass_comp["SO4"].fix(1e-9)
    m.fs.solex_rougher1.mscontactor.aqueous_inlet_state[0].conc_mass_comp["HSO4"].fix(1e-9)
    m.fs.solex_rougher1.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Cl"].fix(1e-9)
    m.fs.solex_rougher1.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Al"].fix(820)
    m.fs.solex_rougher1.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Ca"].fix(5230)
    m.fs.solex_rougher1.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Fe"].fix(270)
    m.fs.solex_rougher1.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Sc"].fix(209.31)
    m.fs.solex_rougher1.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Y"].fix(637.74)
    m.fs.solex_rougher1.mscontactor.aqueous_inlet_state[0].conc_mass_comp["La"].fix(2032.77)
    m.fs.solex_rougher1.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Ce"].fix(4516.13)
    m.fs.solex_rougher1.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Pr"].fix(756.64)
    m.fs.solex_rougher1.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Nd"].fix(2047.85)
    m.fs.solex_rougher1.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Sm"].fix(369.1)
    m.fs.solex_rougher1.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Gd"].fix(174.38)
    m.fs.solex_rougher1.mscontactor.aqueous_inlet_state[0].conc_mass_comp["Dy"].fix(101.12)

    # TODO: This make-up + recycle stream ~ 62.01
    m.fs.organic_make_up.flow_vol.fix(62.01)

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

    # print(dof(m))
    dt = DiagnosticsToolbox(model=m)
    dt.report_structural_issues()

    badly_scaled_var_list = iscale.badly_scaled_var_generator(m, large=1e2, small=1e-2)
    print("----------------   badly_scaled_var_list   ----------------")
    for x in badly_scaled_var_list:
        print(f"{x[0].name}\t{x[0].value}\tsf: {iscale.get_scaling_factor(x[0])}")

    seq = SequentialDecomposition()
    seq.options.tear_method = "Direct"
    seq.options.iterLim = 1
    seq.options.tear_set = [m.fs.mixed_feed]

    G = seq.create_graph(m)
    order = seq.calculation_order(G)
    print("Initialization Order")
    for o in order:
        print(o[0].name)

    tear_guesses1 = {
        "flow_mass": {0.007},
        "conc_mass_comp": {
            "Al": 1493939.39,
            "Ca": 501864.01,
            "Ce": 49698.79,
            "Dy": 466.86,
            "Fe": 1685228.86,
            "Gd": 1624.81,
            "La": 32143.22,
            "Nd": 12552.13,
            "Pr": 4084.40,
            "Sc": 17310.17,
            "Sm": 931.35,
            "Y": 2666.95,
        },
        "flow_mol_comp": {
            "Al": 577.62,
            "Ca": 64.65,
            "Ce": 1.17,
            "Dy": 0.0086,
            "Fe": 231.93,
            "Gd": 0.054,
            "La": 0.63,
            "Nd": 0.34,
            "Pr": 0.099,
            "Sc": 303.33,
            "Sm": 0.021,
            "Y": 0.090,
        },
    }

    # Pass the tear_guess to the SD tool
    seq.set_guesses_for(m.fs.precipitator.cv_aqueous.properties_out[0], tear_guesses1)

    def function(stream):
        initializer_feed = FeedInitializer()
        initializer_product = ProductInitializer()
        initializer1 = MSContactorInitializer()
        initializer2 = BlockTriangularizationInitializer()

        propagate_state(m.fs.liq_feed)
        propagate_state(m.fs.sol_feed)

        if stream == m.fs.leach_liquid_feed:
            initializer_feed.initialize(m.fs.leach_liquid_feed)
        elif stream == m.fs.leach_solid_feed:
            initializer_feed.initialize(m.fs.leach_solid_feed)
        elif stream == m.fs.leach_filter_cake:
            print(f"Initializing {stream}")
            initializer_product.initialize(m.fs.leach_filter_cake)
        elif stream == m.fs.leach:
            print(f"Initializing {stream}")
            try:
                initializer1.initialize(m.fs.leach)
            except:
                # Fix feed states
                m.fs.leach.liquid_inlet.flow_vol.fix()
                m.fs.leach.liquid_inlet.conc_mass_comp.fix()
                m.fs.leach.solid_inlet.flow_mass.fix()
                m.fs.leach.solid_inlet.mass_frac_comp.fix()
                # Re-solve leach unit
                solver = SolverFactory("ipopt")
                solver.solve(m.fs.leach, tee=True)
                # Unfix feed states
                m.fs.leach_liquid_feed.flow_vol.unfix()
                m.fs.leach.liquid_inlet.conc_mass_comp.unfix()
                m.fs.leach.solid_inlet.flow_mass.unfix()
                m.fs.leach.solid_inlet.mass_frac_comp.unfix()
        elif stream == m.fs.leach_mixer:
            initializer2.initialize(m.fs.leach_mixer)
        elif stream == m.fs.solex_rougher.mscontactor:
            print(f"Initializing {stream}")
            initializer2.initialize(m.fs.solex_rougher)
        elif stream == m.fs.solex_cleaner.mscontactor:
            print(f"Initializing {stream}")
            initializer2.initialize(m.fs.solex_cleaner)
        elif stream == m.fs.precipitator:
            print(f"Initializing {stream}")
            try:
                initializer2.initialize(m.fs.precipitator)
            except:
                # Fix feed states
                m.fs.precipitator.cv_aqueous.properties_in[0].flow_vol.fix()
                m.fs.precipitator.cv_aqueous.properties_in[0].conc_mass_comp.fix()
                m.fs.precipitator.cv_precipitate.properties_in[0].flow_mol_comp.fix()
                # Re-solve leach unit
                solver = SolverFactory("ipopt")
                solver.solve(m.fs.precipitator, tee=True)
                # Unfix feed states
                m.fs.precipitator.cv_aqueous.properties_in[0].flow_vol.unfix()
                m.fs.precipitator.cv_aqueous.properties_in[0].conc_mass_comp.unfix()
                m.fs.precipitator.cv_precipitate.properties_in[0].flow_mol_comp.unfix()
        else:
            print(f"Initializing {stream}")
            initializer2.initialize(stream)

    seq.run(m, function)


    print("Numerical issues after initialization")
    dt.report_numerical_issues()

    # Solving of the model

    solver = SolverFactory("ipopt")
    # solver.options["bound_push"] = 1e-8
    # solver.options["mu_init"] = 1e-8
    results = solver.solve(m, tee=True)

    print("Numerical issues after solve")
    dt.report_numerical_issues()

    # Final organic outlet display
    # m.fs.solex.mscontactor.organic[0, 1].conc_mass_comp.display()
    # m.fs.solex.mscontactor.organic[0, 1].conc_mol_comp.display()
    #
    # # Final aqueous outlets display
    # m.fs.solex.mscontactor.aqueous[0, 3].conc_mass_comp.display()
    # m.fs.solex.mscontactor.aqueous[0, 3].conc_mol_comp.display()
    # m.fs.solex_rougher1.display()
    # m.fs.solex_rougher2.display()
    return results


if __name__ == "__main__":
    results = main()
