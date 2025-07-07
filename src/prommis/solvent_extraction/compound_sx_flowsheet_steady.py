#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################

from pyomo.environ import ConcreteModel, units

from idaes.core import (
    FlowDirection,
    FlowsheetBlock,
)
from idaes.core.util import to_json


from idaes.core.solvers import get_solver

from prommis.leaching.leach_solution_properties import LeachSolutionParameters
from prommis.solvent_extraction.ree_og_distribution import REESolExOgParameters
from prommis.solvent_extraction.compound_solvent_extraction import (
    CompoundSolventExtraction,
)
from prommis.solvent_extraction.solvent_extraction_reaction_package import (
    SolventExtractionReactions,
)


def build_model(dosage, number_of_stages):
    """
    Build the model
    """

    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.prop_o = REESolExOgParameters()
    m.fs.leach_soln = LeachSolutionParameters()
    m.fs.reaxn = SolventExtractionReactions()

    m.fs.reaxn.extractant_dosage = dosage

    m.fs.compound_solex = CompoundSolventExtraction(
        number_of_finite_elements=number_of_stages,
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
        heterogeneous_reaction_package=m.fs.reaxn,
        has_holdup=True,
        settler_transformation_method="dae.finite_difference",
        settler_transformation_scheme="BACKWARD",
        settler_finite_elements=4,
    )

    return m


def set_inputs(m, dosage):
    """
    Set inputs to the streams
    """

    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "H2O"].fix(1e6)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "H"].fix(10.75)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "SO4"].fix(100)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "HSO4"].fix(1e4)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Al"].fix(422.375)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Ca"].fix(109.542)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Cl"].fix(1e-7)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Fe"].fix(688.266)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Sc"].fix(0.032)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Y"].fix(0.124)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "La"].fix(0.986)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Ce"].fix(2.277)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Pr"].fix(0.303)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Nd"].fix(0.946)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Sm"].fix(0.097)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Gd"].fix(0.2584)
    m.fs.compound_solex.aqueous_inlet.conc_mass_comp[0, "Dy"].fix(0.047)

    m.fs.compound_solex.aqueous_inlet.flow_vol.fix(62.01)

    m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Kerosene"].fix(820e3)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "DEHPA"].fix(
        975.8e3 * dosage / 100
    )
    m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Al_o"].fix(1.267e-5)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Ca_o"].fix(2.684e-5)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Fe_o"].fix(2.873e-6)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Sc_o"].fix(1.734)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Y_o"].fix(2.179e-5)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "La_o"].fix(0.000105)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Ce_o"].fix(0.00031)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Pr_o"].fix(3.711e-5)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Nd_o"].fix(0.000165)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Sm_o"].fix(1.701e-5)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Gd_o"].fix(3.357e-5)
    m.fs.compound_solex.organic_inlet.conc_mass_comp[0, "Dy_o"].fix(8.008e-6)

    m.fs.compound_solex.organic_inlet.flow_vol.fix(62.01)

    m.fs.compound_solex.mixer[:].unit.mscontactor.aqueous[:, :].temperature.fix(
        305.15 * units.K
    )
    m.fs.compound_solex.mixer[:].unit.mscontactor.organic[:, :].temperature.fix(
        305.15 * units.K
    )

    m.fs.compound_solex.mixer[:].unit.mscontactor.volume[:].fix(0.4 * units.m**3)

    m.fs.compound_solex.organic_settler[:].unit.area.fix(1)
    m.fs.compound_solex.aqueous_settler[:].unit.area.fix(1)
    m.fs.compound_solex.aqueous_settler[:].unit.length.fix(1)
    m.fs.compound_solex.organic_settler[:].unit.length.fix(1)


def main(dosage, number_of_stages):
    """
    Main function to run the model.
    """
    m = build_model(dosage, number_of_stages)
    set_inputs(m, dosage)

    return m


dosage = 5
number_of_stages = 3

if __name__ == "__main__":

    m = main(dosage, number_of_stages)

    m.fs.compound_solex.default_initializer().initialize(m.fs.compound_solex)

    solver = get_solver("ipopt_v2")
    results = solver.solve(m, tee=True)

    to_json(m, fname="compound_solvent_extraction.json", human_read=True)
