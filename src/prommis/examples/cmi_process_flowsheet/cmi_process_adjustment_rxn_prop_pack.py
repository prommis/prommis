#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Reaction package for pH Adjustment in the Acid-Free Dissolution Process of 
Neodymium Magnets.
"""

from pyomo.environ import units as pyunits

from idaes.models.properties.modular_properties.base.generic_reaction import (
    ConcentrationForm,
)
from idaes.models.properties.modular_properties.reactions.dh_rxn import constant_dh_rxn

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
                ("Aq", "Fe(NO3)3"): -1,
                ("Aq", "Nd(NO3)3"): 0,
                ("Aq", "NH4OH"): -3,
                ("Sol", "Fe(OH)3"): 1,
                ("Sol", "Nd(OH)3"): 0,
                ("Aq", "NH4NO3"): 3,
                ("Aq", "Fe(NO3)2"): 0,
                ("Liq", "H2O"): 0,
                ("Vap", "O2"): 0,
                ("Aq", "H2C2O4"): 0,
                ("Sol", "Nd2(C2O4)3 * 10H2O"): 0,
                ("Sol", "Nd2Fe14B"): 0,
                ("Aq", "Cu(NO3)2"): 0,
                ("Sol", "Cu3(BO3)2"): 0,
                ("Sol", "Cu2O"): 0,
                ("Sol", "Cu"): 0,
                ("Aq", "(NH4)3[Fe(C2O4)3]"): 0,
            },
            "heat_of_reaction": constant_dh_rxn,
            "concentration_form": ConcentrationForm.molarity,
            "parameter_data": {"dh_rxn_ref": (-1, pyunits.J / pyunits.mol)},
        },  # Placeholder Value
        "R2": {
            "stoichiometry": {
                ("Aq", "Fe(NO3)3"): 0,
                ("Aq", "Nd(NO3)3"): -1,
                ("Aq", "NH4OH"): -3,
                ("Sol", "Fe(OH)3"): 0,
                ("Sol", "Nd(OH)3"): 1,
                ("Aq", "NH4NO3"): 3,
                ("Aq", "Fe(NO3)2"): 0,
                ("Liq", "H2O"): 0,
                ("Vap", "O2"): 0,
                ("Aq", "H2C2O4"): 0,
                ("Sol", "Nd2(C2O4)3 * 10H2O"): 0,
                ("Sol", "Nd2Fe14B"): 0,
                ("Aq", "Cu(NO3)2"): 0,
                ("Sol", "Cu3(BO3)2"): 0,
                ("Sol", "Cu2O"): 0,
                ("Sol", "Cu"): 0,
                ("Aq", "(NH4)3[Fe(C2O4)3]"): 0,
            },
            "heat_of_reaction": constant_dh_rxn,
            "concentration_form": ConcentrationForm.molarity,
            "parameter_data": {"dh_rxn_ref": (-1, pyunits.J / pyunits.mol)},
        },  # Placeholder Value
    },
}
