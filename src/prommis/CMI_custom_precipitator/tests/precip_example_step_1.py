from pyomo.environ import (
    Constraint,
    Var,
    ConcreteModel,
    Expression,
    Objective,
    TransformationFactory,
    value,
    units as pyunits,
    SolverFactory,
)
from pyomo.network import Arc

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)

from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    InitializationStatus,
)

from idaes.core.util import DiagnosticsToolbox


from idaes.models.unit_models import Feed, StoichiometricReactor, Product
from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock,
)
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)

from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state
from idaes.core.util.scaling import report_scaling_issues
import idaes.logger as idaeslog

import AqueousProperties as aq_thermo_prop_pack
import PrecipitateProperties as precip_thermo_prop_pack
from opt_based_precipitator import Precipitator

# property package
aq_comp_list = ["HNO3", "H^+", "OH^-", "NO3^-", "Fe^3+"]
# precip_comp_list = ["FeOH3", "Fe2O3"]
precip_comp_list = ["FeOH3"]


# reaction package
aq_log_keq_dict = {
    "E1": 1.084, 
    "E2": -14,
    }

precip_log_keq_dict = {
    "E3": -33.498,
    # "E4": -1.763
    }

eq_stoich_dict = {
    "E1": {"HNO3": -1, "H^+": 1, "NO3^-": 1, "OH^-": 0, "Fe^3+": 0},
    "E2": {"H^+": 1, "OH^-": 1, "HNO3": 0, "NO3^-": 0, "Fe^3+": 0},
    "E3": {"OH^-": 3, "Fe^3+": 1, "HNO3": 0, "NO3^-": 0, "H^+": 0},
    # "E4": {"OH^-": 0, "Fe^3+": 2, "HNO3": 0, "NO3^-": 0, "H^+": -6},
    }
precip_eq_stoich_dict = {
    "E1": {"FeOH3": 0},
    "E2": {"FeOH3": 0},
    "E3": {"FeOH3": -1},
    # "E4": {"FeOH3": 0, "Fe2O3": -1},
    }


m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.aq_properties = aq_thermo_prop_pack.AqueousParameter(
    aq_comp_list=aq_comp_list,
    eq_rxn_logkeq_dict=aq_log_keq_dict, 
    eq_rxn_stoich_dict=eq_stoich_dict,
    )
m.fs.precip_properties = precip_thermo_prop_pack.PrecipitateParameter(
    precip_comp_list=precip_comp_list,
    precip_eq_rxn_logkeq_dict=precip_log_keq_dict,
    precip_eq_rxn_stoich_dict=precip_eq_stoich_dict,
)

m.fs.unit = Precipitator(
    property_package_aqueous=m.fs.aq_properties,
    property_package_precipitate=m.fs.precip_properties,
)

# initial concentrations
pH = 7
Hp_init = 10**(-pH)
OHm_init = 10**(-14 + pH)

Feppp_init = 1.05e-1
NO3m_init = 1.8e-1
HNO3m_init = 1e-20
FeOH3_init = 1e-20
Fe2O3_init = 1e-20

m.fs.unit.aqueous_inlet.molality_aq_comp[0, "HNO3"].fix(HNO3m_init)
m.fs.unit.aqueous_inlet.molality_aq_comp[0, "H^+"].fix(Hp_init)
m.fs.unit.aqueous_inlet.molality_aq_comp[0, "OH^-"].fix(OHm_init)
m.fs.unit.aqueous_inlet.molality_aq_comp[0, "NO3^-"].fix(NO3m_init)
m.fs.unit.aqueous_inlet.molality_aq_comp[0, "Fe^3+"].fix(Feppp_init)
m.fs.unit.precipitate_inlet.molality_precip_comp[0, "FeOH3"].fix(FeOH3_init)

m.pprint()

solver_obj = SolverFactory('ipopt', options={"max_iter": 10000, "tol": 1e-8,})


results = solver_obj.solve(m, tee=True)

m.fs.unit.display()
# m.fs.unit.pprint()


# dt = DiagnosticsToolbox(m)
# dt.report_numerical_issues()
# dt.report_structural_issues()
# dt.display_underconstrained_set()
# m.fs.unit.rxn_extent.pprint()
# m.fs.unit.log_equil_rxn_eqns.pprint()
# m.fs.unit.mole_balance_eqns.pprint()

# report_scaling_issues(m.fs.unit)



# print(degrees_of_freedom(m))

# initializer = BlockTriangularizationInitializer(constraint_tolerance=2e-5)
# initializer.initialize(m.fs.unit)

# assert initializer.summary[m.fs.unit]["status"] == InitializationStatus.Ok

# solver_obj = SolverFactory("ipopt", options={"max_iter": 10000, "tol": 1e-8,}) # "linear_solver": "ma57",},)

# solver_obj.solve(
#         m,
#         tee=True,
#     )


# results = solver_obj.solve(m, tee=True)

# m.display()

# m.fs.unit.rxn_extent.pprint()
# # m.fs.unit.mole_balance_eqns.pprint()

# dt.report_structural_issues()
# dt.display_components_with_inconsistent_units()