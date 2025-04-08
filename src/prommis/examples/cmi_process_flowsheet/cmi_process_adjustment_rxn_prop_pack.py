#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Reaction package for pH Adjustment in the Acid-Free Dissolution Process of
Neodymium Magnets.
"""

from pyomo.environ import units as pyunits

config_dict = {
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    "rate_reactions": {
        "R1": {
            "stoichiometry": {
                ("Sol", "Nd2Fe14B"): 0,
                ("Aq", "Cu_2+"): 0,
                ("Vap", "O2"): 0,
                ("Aq", "Nd_3+"): 0,
                ("Aq", "Fe_2+"): 0,
                ("Sol", "Cu3(BO3)2"): 0,
                ("Sol", "Cu2O"): 0,
                ("Sol", "Cu"): 0,
                ("Aq", "Fe_3+"): -1,
                ("Aq", "NO3_-"): 0,
                ("Aq", "H2C2O4"): 0,
                ("Aq", "C2O4_2-"): 0,
                ("Aq", "NH4_+"): 0,
                ("Aq", "OH_-"): -3,
                ("Liq", "H2O"): 0,
                ("Sol", "Fe(OH)3"): 1,
                ("Sol", "Nd(OH)3"): 0,
                ("Sol", "Nd2(C2O4)3 * 10H2O"): 0,
            },
        },
        "R2": {
            "stoichiometry": {
                ("Sol", "Nd2Fe14B"): 0,
                ("Aq", "Cu_2+"): 0,
                ("Vap", "O2"): 0,
                ("Aq", "Nd_3+"): -1,
                ("Aq", "Fe_2+"): 0,
                ("Sol", "Cu3(BO3)2"): 0,
                ("Sol", "Cu2O"): 0,
                ("Sol", "Cu"): 0,
                ("Aq", "Fe_3+"): 0,
                ("Aq", "NO3_-"): 0,
                ("Aq", "H2C2O4"): 0,
                ("Aq", "C2O4_2-"): 0,
                ("Aq", "NH4_+"): 0,
                ("Aq", "OH_-"): -3,
                ("Liq", "H2O"): 0,
                ("Sol", "Fe(OH)3"): 0,
                ("Sol", "Nd(OH)3"): 1,
                ("Sol", "Nd2(C2O4)3 * 10H2O"): 0,
            },
        },
    },
}
