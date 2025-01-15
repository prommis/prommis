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
Property package for CMI Process for recovery of REEs from EOL Neodymium Magnets.
"""

# Import Python libraries
import logging

# Import Pyomo units
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import (
    Anion,
    AqueousPhase,
    Cation,
    Component,
    LiquidPhase,
    SolidPhase,
    Solute,
    Solvent,
    VaporPhase,
)
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.phase_equil.forms import fugacity
from idaes.models.properties.modular_properties.pure.ConstantProperties import Constant
from idaes.models.properties.modular_properties.pure.NIST import NIST
from idaes.models.properties.modular_properties.pure.Perrys import Perrys
from idaes.models.properties.modular_properties.pure.RPP4 import RPP4
from idaes.models.properties.modular_properties.state_definitions import FpcTP

# Set up logger
_log = logging.getLogger(__name__)


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

config_dict = {
    # Specifying components
    "components": {
        'Nd2Fe14B':
            {"type": Component,
             "valid_phase_types": [3], # Solid
             "dens_mol_sol_comp": Constant,
             "cp_mol_sol_comp": Constant,
             "enth_mol_sol_comp": Constant,
             "parameter_data": {
                 "mw": (1.08112, pyunits.kg/pyunits.mol),
                 "dens_mol_sol_comp_coeff": (6.84, pyunits.kmol*pyunits.m**-3),
                 "cp_mol_sol_comp_coeff": (542.2, pyunits.J/pyunits.mol/pyunits.K), # [8]
                 "enth_mol_form_sol_comp_ref":(-93.39E3, pyunits.J/pyunits.mol)}}, # [7]
        'Cu(NO3)2':
            {"type": Solute,
             "valid_phase_types": [4], # Aqueous
             "dens_mol_liq_comp": Constant,
             "cp_mol_liq_comp": Constant,
             "enth_mol_liq_comp": Constant,
             "parameter_data": {
                 "mw": (0.18756, pyunits.kg/pyunits.mol),
                 "dens_mol_liq_comp_coeff": (1, pyunits.kmol*pyunits.m**-3), # PLACEHOLDER
                 "cp_mol_liq_comp_coeff": (1, pyunits.J/pyunits.mol/pyunits.K), # PLACEHOLDER
                 "enth_mol_form_liq_comp_ref":(-1, pyunits.J/pyunits.mol)}}, # PLACEHOLDER
        'O2':
            {"type": Component,
             "valid_phase_types": [2], # Vapor
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
                     "C": (-5.667, pyunits.K)}}},
        'Nd(NO3)3':
            {"type": Solute,
             "valid_phase_types": [4], # Aqueous
             "dens_mol_liq_comp": Constant,
             "cp_mol_liq_comp": Constant,
             "enth_mol_liq_comp": Constant,
             "parameter_data": {
                 "mw": (0.33025, pyunits.kg/pyunits.mol),
                 "dens_mol_liq_comp_coeff": (1, pyunits.kmol*pyunits.m**-3), # PLACEHOLDER
                 "cp_mol_liq_comp_coeff": (1, pyunits.J/pyunits.mol/pyunits.K), # PLACEHOLDER
                 "enth_mol_form_liq_comp_ref":(1, pyunits.J/pyunits.mol)}}, # PLACEHOLDER
        'Fe(NO3)2':
            {"type": Solute,
             "valid_phase_types": [4], # Aqueous
             "dens_mol_liq_comp": Constant,
             "cp_mol_liq_comp": Constant,
             "enth_mol_liq_comp": Constant,
             "parameter_data": {
                 "dens_mol_liq_comp_coeff": (1, pyunits.kmol*pyunits.m**-3), # PLACEHOLDER
                 "cp_mol_liq_comp_coeff": (1, pyunits.J/pyunits.mol/pyunits.K), # PLACEHOLDER
                 "enth_mol_form_liq_comp_ref":(-1, pyunits.J/pyunits.mol)}}, # PLACEHOLDER
        'Cu3(BO3)2':
            {"type": Component,
             "valid_phase_types": [3], # Solid
             "dens_mol_sol_comp": Constant,
             "cp_mol_sol_comp": Constant,
             "enth_mol_sol_comp": Constant,
             "parameter_data": {
                 "mw": (0.29134, pyunits.kg/pyunits.mol),
                 "dens_mol_sol_comp_coeff": (6.84, pyunits.kmol*pyunits.m**-3), # PLACEHOLDER
                 "cp_mol_sol_comp_coeff": (542.2, pyunits.J/pyunits.mol/pyunits.K), # PLACEHOLDER
                 "enth_mol_form_sol_comp_ref":(0, pyunits.J/pyunits.mol)}}, # PLACEHOLDER
        'Cu2O':
            {"type": Component,
             "valid_phase_types": [3], # Solid
             "dens_mol_sol_comp": Constant,
             "cp_mol_sol_comp": Constant,
             "enth_mol_sol_comp": Constant,
             "parameter_data": {
                 "mw": (0.14309, pyunits.kg/pyunits.mol),
                 "dens_mol_sol_comp_coeff": (6.84, pyunits.kmol*pyunits.m**-3), # PLACEHOLDER
                 "cp_mol_sol_comp_coeff": (542.2, pyunits.J/pyunits.mol/pyunits.K), # PLACEHOLDER
                 "enth_mol_form_sol_comp_ref":(-180E3, pyunits.J/pyunits.mol)}},  # PLACEHOLDER
        'Cu':
            {"type": Component,
             "valid_phase_types": [3], # Solid
             "dens_mol_sol_comp": Constant,
             "cp_mol_sol_comp": Constant,
             "enth_mol_sol_comp": Constant,
             "parameter_data": {
                 "mw": (0.063546, pyunits.kg/pyunits.mol),
                 "dens_mol_sol_comp_coeff": (6.84, pyunits.kmol*pyunits.m**-3), # PLACEHOLDER
                 "cp_mol_sol_comp_coeff": (542.2, pyunits.J/pyunits.mol/pyunits.K), # PLACEHOLDER
                 "enth_mol_form_sol_comp_ref":(0, pyunits.J/pyunits.mol)}}, # PLACEHOLDER
        'H2O':
            {"type": Component,
             "valid_phase_types": [1], # Liquid
             "dens_mol_liq_comp": Perrys,
             "enth_mol_liq_comp": Perrys,
             "enth_mol_ig_comp": RPP4,
             "pressure_sat_comp": RPP4,
             "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
             "parameter_data": {
                 "mw": (18.015e-3, pyunits.kg / pyunits.mol),  # [1] pg. 667
                 "pressure_crit": (221.2e5, pyunits.Pa),  # [1] pg. 667
                 "temperature_crit": (647.3, pyunits.K),  # [1] pg. 667
                 "dens_mol_liq_comp_coeff": {  # [2] pg. 2-98
                     "eqn_type": 1,
                     "1": (5.459, pyunits.kmol * pyunits.m**-3),
                     "2": (0.30542, None),
                     "3": (647.13, pyunits.K),
                     "4": (0.081, None),
                 },
                 "cp_mol_ig_comp_coeff": {  # [1] pg. 668
                     "A": (3.224e1, pyunits.J / pyunits.mol / pyunits.K),
                     "B": (1.924e-3, pyunits.J / pyunits.mol / pyunits.K**2),
                     "C": (1.055e-5, pyunits.J / pyunits.mol / pyunits.K**3),
                     "D": (-3.596e-9, pyunits.J / pyunits.mol / pyunits.K**4),
                 },
                 "cp_mol_liq_comp_coeff": {  # [2] pg. 2-174
                     "1": (2.7637e5, pyunits.J / pyunits.kmol / pyunits.K),
                     "2": (-2.0901e3, pyunits.J / pyunits.kmol / pyunits.K**2),
                     "3": (8.1250, pyunits.J / pyunits.kmol / pyunits.K**3),
                     "4": (-1.4116e-2, pyunits.J / pyunits.kmol / pyunits.K**4),
                     "5": (9.3701e-6, pyunits.J / pyunits.kmol / pyunits.K**5),
                 },
                 "enth_mol_form_liq_comp_ref": (
                     -285.830e3,
                     pyunits.J / pyunits.mol,
                 ),  # [3] updated 5/10/24
                 "enth_mol_form_vap_comp_ref": (
                     -241.826e3,
                     pyunits.J / pyunits.mol,
                 ),  # [3] updated 5/10/24
                 "pressure_sat_comp_coeff": {
                     "A": (-7.76451, None),  # [1] pg. 669
                     "B": (1.45838, None),
                     "C": (-2.77580, None),
                     "D": (-1.23303, None)}}},
        'Fe(NO3)3':
            {"type": Solute,
             "valid_phase_types": [4], # Aqueous
             "dens_mol_liq_comp": Constant,
             "cp_mol_liq_comp": Constant,
             "enth_mol_liq_comp": Constant,
             "parameter_data": {
                 "mw": (0.24186, pyunits.kg/pyunits.mol),
                 "dens_mol_liq_comp_coeff": (1, pyunits.kmol*pyunits.m**-3), # PLACEHOLDER
                 "cp_mol_liq_comp_coeff": (1, pyunits.J/pyunits.mol/pyunits.K), # PLACEHOLDER
                 "enth_mol_form_liq_comp_ref":(-1, pyunits.J/pyunits.mol)}}, # PLACEHOLDER
        'Fe(OH)3':
            {"type": Component,
             "valid_phase_types": [3], # Solid
             "dens_mol_sol_comp": Constant,
             "cp_mol_sol_comp": Constant,
             "enth_mol_sol_comp": Constant,
             "parameter_data": {
                 "mw": (0.10687, pyunits.kg/pyunits.mol),
                 "dens_mol_sol_comp_coeff": (6.84, pyunits.kmol*pyunits.m**-3), # PLACEHOLDER
                 "cp_mol_sol_comp_coeff": (542.2, pyunits.J/pyunits.mol/pyunits.K), # PLACEHOLDER
                 "enth_mol_form_sol_comp_ref":(-826E3, pyunits.J/pyunits.mol)}}, # PLACEHOLDER
        'NH4OH':
            {"type": Solute,
             "valid_phase_types": [4], # Aqueous
             "dens_mol_liq_comp": Constant,
             "cp_mol_liq_comp": Constant,
             "enth_mol_liq_comp": Constant,
             "parameter_data": {
                 "mw": (0.03505, pyunits.kg/pyunits.mol),
                 "dens_mol_liq_comp_coeff": (1, pyunits.kmol*pyunits.m**-3), # PLACEHOLDER
                 "cp_mol_liq_comp_coeff": (1, pyunits.J/pyunits.mol/pyunits.K), # PLACEHOLDER
                 "enth_mol_form_liq_comp_ref":(-1, pyunits.J/pyunits.mol)}}, # PLACEHOLDER
        'Nd(OH)3':
            {"type": Component,
             "valid_phase_types": [3], # Solid
             "dens_mol_sol_comp": Constant,
             "cp_mol_sol_comp": Constant,
             "enth_mol_sol_comp": Constant,
             "parameter_data": {
                 "mw": (0.19526, pyunits.kg/pyunits.mol),
                 "dens_mol_sol_comp_coeff": (6.84, pyunits.kmol*pyunits.m**-3), # PLACEHOLDER
                 "cp_mol_sol_comp_coeff": (542.2, pyunits.J/pyunits.mol/pyunits.K), # PLACEHOLDER
                 "enth_mol_form_sol_comp_ref":(-1400E3, pyunits.J/pyunits.mol)}}, # PLACEHOLDER
        'NH4NO3': 
            {"type": Solute,
             "valid_phase_types": [4], # Aqueous
             "dens_mol_liq_comp": Constant,
             "cp_mol_liq_comp": Constant,
             "enth_mol_liq_comp": Constant,
             "parameter_data": {
                 "mw": (0.080043, pyunits.kg/pyunits.mol),
                 "dens_mol_liq_comp_coeff": (1, pyunits.kmol*pyunits.m**-3), # PLACEHOLDER
                 "cp_mol_liq_comp_coeff": (1, pyunits.J/pyunits.mol/pyunits.K), # PLACEHOLDER
                 "enth_mol_form_liq_comp_ref":(-1, pyunits.J/pyunits.mol)}}, # PLACEHOLDER
        'H2C2O4':
            {"type": Solute,
             "valid_phase_types": [4], # Aqueous
             "dens_mol_liq_comp": Constant,
             "cp_mol_liq_comp": Constant,
             "enth_mol_liq_comp": Constant,
             "parameter_data": {
                 "mw": (0.09003, pyunits.kg/pyunits.mol),
                 "dens_mol_liq_comp_coeff": (1, pyunits.kmol*pyunits.m**-3), # PLACEHOLDER
                 "cp_mol_liq_comp_coeff": (1, pyunits.J/pyunits.mol/pyunits.K), # PLACEHOLDER
                 "enth_mol_form_liq_comp_ref":(-1, pyunits.J/pyunits.mol)}}, # PLACEHOLDER
        '(NH4)3[Fe(C2O4)3]':
            {"type": Solute,
             "valid_phase_types": [4], # Aqueous
             "dens_mol_liq_comp": Constant,
             "cp_mol_liq_comp": Constant,
             "enth_mol_liq_comp": Constant,
             "parameter_data": {
                 "mw": (0.197983, pyunits.kg/pyunits.mol),
                 "dens_mol_liq_comp_coeff": (1, pyunits.kmol*pyunits.m**-3), # PLACEHOLDER
                 "cp_mol_liq_comp_coeff": (1, pyunits.J/pyunits.mol/pyunits.K), # PLACEHOLDER
                 "enth_mol_form_liq_comp_ref":(-1, pyunits.J/pyunits.mol)}}, # PLACEHOLDER
        'Nd2(C2O4)3 * 10H2O':
            {"type": Component,
             "valid_phase_types": [3], # Solid
             "dens_mol_sol_comp": Constant,
             "cp_mol_sol_comp": Constant,
             "enth_mol_sol_comp": Constant,
             "parameter_data": {
                 "mw": (0.7326938, pyunits.kg/pyunits.mol), 
                 "dens_mol_sol_comp_coeff": (6.84, pyunits.kmol*pyunits.m**-3), # PLACEHOLDER
                 "cp_mol_sol_comp_coeff": (542.2, pyunits.J/pyunits.mol/pyunits.K), # PLACEHOLDER
                 "enth_mol_form_sol_comp_ref":(0, pyunits.J/pyunits.mol)}}}, # PLACEHOLDER

    # Specifying phases
    "phases":  {'Liq': {"type": LiquidPhase,
                        "equation_of_state": Ideal},
                'Sol': {"type": SolidPhase,
                        "equation_of_state": Ideal},
                'Aq': {"type": AqueousPhase,
                        "equation_of_state": Ideal},
                'Vap': {"type": VaporPhase,
                        "equation_of_state": Ideal}},

    # Set base units of measurement
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},

    # Specifying state definition
    "state_definition": FpcTP,
    "state_bounds": {"flow_mol_phase_comp": (0, 100, 100000,
                                             pyunits.mol/pyunits.s),
                     "temperature": (273.15, 298.15, 1500, pyunits.K),
                     "pressure": (5e4, 1e5, 1e6, pyunits.Pa)},
    "pressure_ref": (1e5, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K)}