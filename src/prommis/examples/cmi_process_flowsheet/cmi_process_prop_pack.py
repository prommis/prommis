#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Property package for the Critical Minerals Innovation Hub process.
"""

# Import Pyomo units
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import (
    AqueousPhase,
    Component,
    LiquidPhase,
    SolidPhase,
    Solute,
    Solvent,
    VaporPhase,
)
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.pure.ConstantProperties import Constant
from idaes.models.properties.modular_properties.pure.NIST import NIST
from idaes.models.properties.modular_properties.pure.Perrys import Perrys
from idaes.models.properties.modular_properties.pure.RPP4 import RPP4
from idaes.models.properties.modular_properties.state_definitions import FpcTP

# all ion heat capacities approximated as 75000 kJ/kmol/K.

# ---------------------------------------------------------------------

# Data Sources:
# [1] The Properties of Gases and Liquids (1987)
#     4th edition, Chemical Engineering Series - Robert C. Reid
# [2] Perry's Chemical Engineers' Handbook 7th Ed.
# [3] NIST Chemistry WebBook, https://webbook.nist.gov/chemistry/
#     Retrieved 23rd September, 2021
# [4] B. Ruscic and D. H. Bross, Active Thermochemical Tables (ATcT)
#     values based on ver. 1.122r of the Thermochemical Network (2021);
#     available at ATcT.anl.gov
# [5] CRC Handbook of Chemistry and Physics, 97th Ed., W.M. Haynes
# [6] Journal of Physical and Chemical Reference Data 20, 1157
#     (1991); https:// doi.org/10.1063/1.555899
# [7] https://www.sciencedirect.com/science/article/pii/S0040603120301696
# [8] https://e-magnetsuk.com/introduction-to-neodymium-magnets/characteristics-of-ndfeb-magnets/
# [9] https://ecampusontario.pressbooks.pub/genchemforgeegees/chapter/appendix-g-standard-enthalpies-of-formation-for-selected-substances/
# [10] https://sistemas.eel.usp.br/docentes/arquivos/5817712/TDQ%20I/R-standard_enthalpy_of_formation.pdf
# [11] https://www.engineeringtoolbox.com/standard-state-enthalpy-formation-definition-value-Gibbs-free-energy-entropy-molar-heat-capacity-d_1978.html
# [12] https://www.osti.gov/servlets/purl/
# [13] https://en.wikipedia.org/wiki/Copper(II)_borate
# [14] https://en.wikipedia.org/wiki/Table_of_specific_heat_capacities
# [15] https://www.americanelements.com/copper-i-oxide-1317-39-1
# [16] https://www.boyiprototyping.com/materials-guide/density-of-copper/
# [17] https://winter.group.shef.ac.uk/webelements/copper/thermochemistry.html
# [18] https://chemister.ru/Databases/Chemdatabase/properties-en.php?dbid=1&id=450
# [19] https://srd.nist.gov/jpcrdreprint/1.555964.pdf
# [20] https://www.osti.gov/servlets/purl/1530415
# [21] https://www.chemicalbook.com/ChemicalProductProperty_EN_CB0323998.htm
# [22] https://pubmed.ncbi.nlm.nih.gov/30901227/
# [23] https://www.chemicalbook.com/ChemicalProductProperty_US_CB9262164.aspx
# [24] https://deepblue.lib.umich.edu/bitstream/handle/2027.42/23403/0000348.pdf;sequence=1
# [25] https://www.osti.gov/servlets/purl/10111965


thermo_config = {
    "components": {
        ### Liquids
        "H2O": {
            "type": Solvent,
            "valid_phase_types": [1],  # Liquid
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (18.0153, pyunits.g / pyunits.mol),
                # Parameters here come from Perry's Handbook:  p. 2-98
                "dens_mol_liq_comp_coeff": {
                    "eqn_type": 1,
                    "1": (5.459, pyunits.kmol * pyunits.m**-3),
                    "2": (0.30542, pyunits.dimensionless),
                    "3": (647.13, pyunits.K),
                    "4": (0.081, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-285.830, pyunits.kJ / pyunits.mol),
                "enth_mol_form_vap_comp_ref": (0, pyunits.kJ / pyunits.mol),
                # Parameters here come Perry's Handbook:  p. 2-174
                "cp_mol_liq_comp_coeff": {
                    "1": (2.7637e5, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (-2.0901e3, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (8.125, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (-1.4116e-2, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (9.3701e-6, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "cp_mol_ig_comp_coeff": {
                    "A": (30.09200, pyunits.J / pyunits.mol / pyunits.K),
                    "B": (
                        6.832514,
                        pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-1,
                    ),
                    "C": (
                        6.793435,
                        pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-2,
                    ),
                    "D": (
                        -2.534480,
                        pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**-3,
                    ),
                    "E": (
                        0.082139,
                        pyunits.J * pyunits.mol**-1 * pyunits.K**-1 * pyunits.kiloK**2,
                    ),
                    "F": (-250.8810, pyunits.kJ / pyunits.mol),
                    "G": (223.3967, pyunits.J / pyunits.mol / pyunits.K),
                    "H": (0, pyunits.kJ / pyunits.mol),
                },
                "entr_mol_form_liq_comp_ref": (
                    69.95,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
        },
        ### Aqueous
        "H2C2O4": {
            "type": Solute,
            "valid_phase_types": [4],  # Aqueous
            "dens_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "parameter_data": {
                "mw": (0.09003, pyunits.kg / pyunits.mol),
                "dens_mol_liq_comp_coeff": (
                    11,
                    pyunits.kmol * pyunits.m**-3,
                ),  # [21]
                "cp_mol_liq_comp_coeff": (
                    146.0,
                    pyunits.J / pyunits.mol / pyunits.K,
                ),  # [3]
                "enth_mol_form_liq_comp_ref": (-829.7e3, pyunits.J / pyunits.mol),
            },  # [22]
        },
        "OH_-": {
            "type": Solute,
            "valid_phase_types": [4],  # Aqueous
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            "parameter_data": {
                "mw": (17.008, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55, pyunits.kmol * pyunits.m**-3),
                "enth_mol_form_liq_comp_ref": (-230.000, pyunits.kJ / pyunits.mol),
                "cp_mol_liq_comp_coeff": (75000, pyunits.J / pyunits.kmol / pyunits.K),
                "entr_mol_form_liq_comp_ref": (
                    -10.75,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
        },
        "Cu_2+": {
            "type": Solute,
            "valid_phase_types": [4],  # Aqueous
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            "parameter_data": {
                "mw": (0.06355, pyunits.kg / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55, pyunits.kmol * pyunits.m**-3),
                "enth_mol_form_liq_comp_ref": (64.4, pyunits.kJ / pyunits.mol),  # [10]
                "cp_mol_liq_comp_coeff": (75000, pyunits.J / pyunits.kmol / pyunits.K),
                "entr_mol_form_liq_comp_ref": (
                    -98.7,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),  # [11]
            },
        },
        "NO3_-": {
            "type": Solute,
            "valid_phase_types": [4],  # Aqueous
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            "parameter_data": {
                "mw": (0.06201, pyunits.kg / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55, pyunits.kmol * pyunits.m**-3),
                "enth_mol_form_liq_comp_ref": (-207.4, pyunits.kJ / pyunits.mol),  # [3]
                "cp_mol_liq_comp_coeff": (75000, pyunits.J / pyunits.kmol / pyunits.K),
                "entr_mol_form_liq_comp_ref": (
                    146.4,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),  # [3]
            },
        },
        "Nd_3+": {
            "type": Solute,
            "valid_phase_types": [4],  # Aqueous
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            "parameter_data": {
                "mw": (0.14424, pyunits.kg / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55, pyunits.kmol * pyunits.m**-3),
                "enth_mol_form_liq_comp_ref": (-672, pyunits.kJ / pyunits.mol),  # [12]
                "cp_mol_liq_comp_coeff": (75000, pyunits.J / pyunits.kmol / pyunits.K),
                "entr_mol_form_liq_comp_ref": (
                    -315,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),  # [12]
            },
        },
        "Fe_2+": {
            "type": Solute,
            "valid_phase_types": [4],  # Aqueous
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            "parameter_data": {
                "mw": (0.05585, pyunits.kg / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55, pyunits.kmol * pyunits.m**-3),
                "enth_mol_form_liq_comp_ref": (-87.9, pyunits.kJ / pyunits.mol),  # [3]
                "cp_mol_liq_comp_coeff": (75000, pyunits.J / pyunits.kmol / pyunits.K),
                "entr_mol_form_liq_comp_ref": (
                    -137.7,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),  # [3]
            },
        },
        "Fe_3+": {
            "type": Solute,
            "valid_phase_types": [4],  # Aqueous
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            "parameter_data": {
                "mw": (0.05585, pyunits.kg / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55, pyunits.kmol * pyunits.m**-3),
                "enth_mol_form_liq_comp_ref": (-49, pyunits.kJ / pyunits.mol),  # [19]
                "cp_mol_liq_comp_coeff": (75000, pyunits.J / pyunits.kmol / pyunits.K),
                "entr_mol_form_liq_comp_ref": (
                    -278.4,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),  # [3]
            },
        },
        "NH4_+": {
            "type": Solute,
            "valid_phase_types": [4],  # Aqueous
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            "parameter_data": {
                "mw": (0.018039, pyunits.kg / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55, pyunits.kmol * pyunits.m**-3),
                "enth_mol_form_liq_comp_ref": (-132.5, pyunits.kJ / pyunits.mol),  # [5]
                "cp_mol_liq_comp_coeff": (75000, pyunits.J / pyunits.kmol / pyunits.K),
                "entr_mol_form_liq_comp_ref": (
                    -113.4,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),  # [5]
            },
        },
        "C2O4_2-": {
            "type": Solute,
            "valid_phase_types": [4],  # Aqueous
            "dens_mol_liq_comp": Constant,
            "enth_mol_liq_comp": Constant,
            "cp_mol_liq_comp": Constant,
            "entr_mol_liq_comp": Constant,
            "parameter_data": {
                "mw": (0.08802, pyunits.kg / pyunits.mol),
                "dens_mol_liq_comp_coeff": (55, pyunits.kmol * pyunits.m**-3),
                "enth_mol_form_liq_comp_ref": (
                    -825.1,
                    pyunits.kJ / pyunits.mol,
                ),  # [20]
                "cp_mol_liq_comp_coeff": (75000, pyunits.J / pyunits.kmol / pyunits.K),
                "entr_mol_form_liq_comp_ref": (
                    146.4,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),  # [3]
            },
        },
        ### Gases
        "O2": {
            "type": Component,
            "valid_phase_types": [2],  # Vapor
            "enth_mol_ig_comp": RPP4,
            "entr_mol_ig_comp": RPP4,
            "pressure_sat_comp": NIST,
            "parameter_data": {
                "mw": (31.999e-3, pyunits.kg / pyunits.mol),  # [1]
                "pressure_crit": (50.43e5, pyunits.Pa),  # [1]
                "temperature_crit": (154.58, pyunits.K),  # [1]
                "omega": 0.025,  # [1]
                "cp_mol_ig_comp_coeff": {
                    "A": (2.811e1, pyunits.J / pyunits.mol / pyunits.K),
                    "B": (-3.680e-6, pyunits.J / pyunits.mol / pyunits.K**2),
                    "C": (1.746e-5, pyunits.J / pyunits.mol / pyunits.K**3),
                    "D": (-1.065e-8, pyunits.J / pyunits.mol / pyunits.K**4),
                },
                "enth_mol_form_vap_comp_ref": (0.0, pyunits.J / pyunits.mol),  # [3]
                "entr_mol_form_vap_comp_ref": (
                    205.152,
                    pyunits.J / pyunits.mol / pyunits.K,
                ),  # [3]
                "pressure_sat_comp_coeff": {
                    "A": (3.85845, None),  # [3]
                    "B": (325.675, pyunits.K),
                    "C": (-5.667, pyunits.K),
                },
            },
        },
        ### Solids
        "Nd2Fe14B": {
            "type": Component,
            "valid_phase_types": [3],  # Solid
            "dens_mol_sol_comp": Constant,
            "cp_mol_sol_comp": Constant,
            "enth_mol_sol_comp": Constant,
            "parameter_data": {
                "mw": (1.08112, pyunits.kg / pyunits.mol),
                "dens_mol_sol_comp_coeff": (6.84, pyunits.kmol * pyunits.m**-3),
                "cp_mol_sol_comp_coeff": (
                    542.2,
                    pyunits.J / pyunits.mol / pyunits.K,
                ),  # [8]
                "enth_mol_form_sol_comp_ref": (-93.39e3, pyunits.J / pyunits.mol),
            },
        },  # [7]
        "Cu3(BO3)2": {
            "type": Component,
            "valid_phase_types": [3],  # Solid
            "dens_mol_sol_comp": Constant,
            "cp_mol_sol_comp": Constant,
            "enth_mol_sol_comp": Constant,
            "parameter_data": {
                "mw": (0.3083, pyunits.kg / pyunits.mol),
                "dens_mol_sol_comp_coeff": (
                    14.72,
                    pyunits.kmol * pyunits.m**-3,
                ),  # [13]
                "cp_mol_sol_comp_coeff": (
                    24.47,
                    pyunits.J / pyunits.mol / pyunits.K,
                ),  # approximated as Cu [14]
                "enth_mol_form_sol_comp_ref": (-170, pyunits.kJ / pyunits.mol),
            },  # approximated as Cu2O [15]
        },
        "Cu2O": {
            "type": Component,
            "valid_phase_types": [3],  # Solid
            "dens_mol_sol_comp": Constant,
            "cp_mol_sol_comp": Constant,
            "enth_mol_sol_comp": Constant,
            "parameter_data": {
                "mw": (0.14309, pyunits.kg / pyunits.mol),
                "dens_mol_sol_comp_coeff": (
                    41.93,
                    pyunits.kmol * pyunits.m**-3,
                ),  # [15]
                "cp_mol_sol_comp_coeff": (
                    24.47,
                    pyunits.J / pyunits.mol / pyunits.K,
                ),  # approximated as Cu [14]
                "enth_mol_form_sol_comp_ref": (-170, pyunits.kJ / pyunits.mol),
            },  # [15]
        },
        "Cu": {
            "type": Component,
            "valid_phase_types": [3],  # Solid
            "dens_mol_sol_comp": Constant,
            "cp_mol_sol_comp": Constant,
            "enth_mol_sol_comp": Constant,
            "parameter_data": {
                "mw": (0.063546, pyunits.kg / pyunits.mol),
                "dens_mol_sol_comp_coeff": (
                    141,
                    pyunits.kmol * pyunits.m**-3,
                ),  # [16]
                "cp_mol_sol_comp_coeff": (
                    23.43,
                    pyunits.J / pyunits.mol / pyunits.K,
                ),  # [17]
                "enth_mol_form_sol_comp_ref": (0, pyunits.J / pyunits.mol),
            },  # [17]
        },
        "Fe(OH)3": {
            "type": Component,
            "valid_phase_types": [3],  # Solid
            "dens_mol_sol_comp": Constant,
            "cp_mol_sol_comp": Constant,
            "enth_mol_sol_comp": Constant,
            "parameter_data": {
                "mw": (0.10687, pyunits.kg / pyunits.mol),
                "dens_mol_sol_comp_coeff": (
                    34.6,
                    pyunits.kmol * pyunits.m**-3,
                ),  # [18]
                "cp_mol_sol_comp_coeff": (
                    101.7,
                    pyunits.J / pyunits.mol / pyunits.K,
                ),  # [18]
                "enth_mol_form_sol_comp_ref": (-826, pyunits.kJ / pyunits.mol),
            },  # [18]
        },
        "Nd(OH)3": {
            "type": Component,
            "valid_phase_types": [3],  # Solid
            "dens_mol_sol_comp": Constant,
            "cp_mol_sol_comp": Constant,
            "enth_mol_sol_comp": Constant,
            "parameter_data": {
                "mw": (0.19526, pyunits.kg / pyunits.mol),
                "dens_mol_sol_comp_coeff": (
                    23.88,
                    pyunits.kmol * pyunits.m**-3,
                ),  # [23]
                "cp_mol_sol_comp_coeff": (
                    542.2,
                    pyunits.J / pyunits.mol / pyunits.K,
                ),  # [24]
                "enth_mol_form_sol_comp_ref": (-1400, pyunits.kJ / pyunits.mol),
            },  # [25]
        },
        "Nd2(C2O4)3 * 10H2O": {
            "type": Component,
            "valid_phase_types": [3],  # Solid
            "dens_mol_sol_comp": Constant,
            "cp_mol_sol_comp": Constant,
            "enth_mol_sol_comp": Constant,
            "parameter_data": {
                "mw": (0.7326938, pyunits.kg / pyunits.mol),
                "dens_mol_sol_comp_coeff": (
                    23.88,
                    pyunits.kmol * pyunits.m**-3,
                ),  # Estimated as Nd(OH)3 [23]
                "cp_mol_sol_comp_coeff": (
                    542.2,
                    pyunits.J / pyunits.mol / pyunits.K,
                ),  # Estimated as Nd(OH)3 [24]
                "enth_mol_form_sol_comp_ref": (-1400, pyunits.kJ / pyunits.mol),
            },  # Estimated as Nd(OH)3 [25]
        },
    },
    ### Specifying phases
    "phases": {
        "Liq": {"type": LiquidPhase, "equation_of_state": Ideal},
        "Sol": {"type": SolidPhase, "equation_of_state": Ideal},
        "Aq": {"type": AqueousPhase, "equation_of_state": Ideal},
        "Vap": {"type": VaporPhase, "equation_of_state": Ideal},
    },
    ### Set base units of measurement
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    ### Specifying state definition
    "state_definition": FpcTP,
    "state_bounds": {
        "flow_mol_phase_comp": (0, 100, 100000, pyunits.mol / pyunits.s),
        "temperature": (273.15, 298.15, 1500, pyunits.K),
        "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
    },
    "pressure_ref": (1e5, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),
}
