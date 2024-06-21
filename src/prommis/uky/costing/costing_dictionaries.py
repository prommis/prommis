#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Python script to read costing components
This script reads the library of costing components (scaled cost, reference
parameters, costing exponents, etc.) from the json files.
First, open json file, then create a python dictionary that gets imported into
power_plant_capcost.py

Python dictionaries that are loaded:
* REE_costing_params

"""

__author__ = "Costing Team (B. Paul, A. Fritz, A. Ojo, A. Dasgupta, and M. Zamarripa)"
__version__ = "1.0.0"

import json
import os

from pyomo.common.fileutils import this_file_dir

directory = this_file_dir()


def load_REE_costing_dictionary():
    """
    The costing data dictionary contains information from the University
    of Kentucky pilot study "Pilot-Scale Testing of an Integrated Circuit
    for the Extraction of Rare Earth Minerals and Elements from Coal and
    Coal Byproducts Using Advanced Separation Technologies" (2021) and
    from the NETL Quality Guidelines for Energy Systems Studies (Feb 2021.
    Specifically it includes scaling exponents, valid ranges for the
    scaled parameter, and units for those ranges. It is important to note
    the units only apply to the ranges and are not necessarily the units
    that the reference parameter value will be given in.. It includes the
    total plant cost (TPC), reference parameter value, and units for that
    value.

    This dictionary is nested with the following structure:
    source --> account --> property name --> property values
    """
    with open(os.path.join(directory, "REE_costing_parameters.json"), "r") as file:
        REE_costing_parameters = json.load(file)
    return REE_costing_parameters
