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

__author__ = (
    "Costing Team (B. Paul, A. Fritz, A. Ojo, A. Dasgupta, L. Deng, and M. Zamarripa)"
)
__version__ = "1.0.0"

import json
import os

from pyomo.common.fileutils import this_file_dir
from pyomo.environ import units as pyunits

from idaes.core import register_idaes_currency_units

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


def load_default_sale_prices():
    """
    Dictionary of default prices
    MUSD: the currency units are millions of USD, so its price need a 1e-6 multiplier to get USD
    """
    register_idaes_currency_units()
    CE_index_units = pyunits.MUSD_2023
    default_sale_prices = {
        # pure elements. 33 pure elements prices are cited from https://www.usgs.gov/centers/national-minerals-information-center/minerals-yearbook-metals-and-minerals, year 2023, unless otherwise noticed.
        # minimum purity of 99.9% for pure and oxides is assumed, unless mentioned otherwise.
        "Al": 2.78 * 1e-6 * CE_index_units / pyunits.kg,
        "Sb": 12.1 * 1e-6 * CE_index_units / pyunits.kg,
        "As": 4.52 * 1e-6 * CE_index_units / pyunits.kg,
        "Ba": 0.22 * 1e-6 * CE_index_units / pyunits.kg,
        "Be": 1400 * 1e-6 * CE_index_units / pyunits.kg,
        "Bi": 8.99 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Ce price not available in usgs, but USGS has CeO2's price.
        "Ce": 100
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.5%.1kg price. https://www.luciteria.com/elements-for-sale/buy-cerium.
        "Cs": 91600 * 1e-6 * CE_index_units / pyunits.kg,
        "Cr": 11.13 * 1e-6 * CE_index_units / pyunits.kg,
        "Co": 34.13 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure element Dysprosium (Dy) price is not available in USGS, but USGS has Dy2O3's price.
        "Dy": 1900
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.99%. https://www.luciteria.com/elements-for-sale/dysprosium-metal-9999-dendritic. 1kg price.
        # Pure Erbium's (Er) price is not available in USGS, nor its oxides' price.
        "Er": 39
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.99%. https://www.luciteria.com/elements-for-sale/erbium-metal-9999-dendritic. 1kg price.
        # Pure Europium (Eu) price is not available in USGS, but USGS has Eu2O3's price.
        "Eu": 1550
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.9%. https://www.luciteria.com/elements-for-sale/buy-europium. 1kg price.
        "CaF2": 0.296
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Fluorspar, alternate Names: Fluorite, Calcium Fluoride.
        # Pure Gadolinium's (Gd) price is not available in USGS, nor its oxides' price.
        "Gd": 850
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.95%. https://www.luciteria.com/elements-for-sale/gadolinium-metal-9995-dendritic. 1kg price.
        "Ga": 450 * 1e-6 * CE_index_units / pyunits.kg,
        "Ge": 1392 * 1e-6 * CE_index_units / pyunits.kg,
        "C": 1.08 * 1e-6 * CE_index_units / pyunits.kg,
        "Ha": 6150 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Holmium's (Ho) price is not available in USGS, nor its oxides' price.
        "Ho": 1600
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.95%. https://www.luciteria.com/elements-for-sale/holmium-metal-9995-dendritic. 1kg price.
        "In": 244 * 1e-6 * CE_index_units / pyunits.kg,
        "Ir": 150233.37 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Lathanum's (Ho) price is not available in USGS, but USGS has La2O3's price.
        "La": 110
        * 1e-6
        * CE_index_units
        / pyunits.kg,  #  https://www.luciteria.com/elements-for-sale/buy-lanthanum. Purity 99.5%, 1kg price.
        "Li": 41.3 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Lutetium's (Lu) price is not available in USGS, nor its oxides' price.
        "Lu": 3600
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.9%. https://www.luciteria.com/elements-for-sale/buy-lutetium. 1kg price.
        "Mg": 11.02 * 1e-6 * CE_index_units / pyunits.kg,
        "Mn": 0.0048 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Neodymium's (Lu) price is not available in USGS, but USGS has Nd2O3's price.
        "Nd": 425
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.5%. https://www.luciteria.com/elements-for-sale/buy-neodymium-metal. 1kg price.
        "Ni": 21.495 * 1e-6 * CE_index_units / pyunits.kg,
        "Nb": 25 * 1e-6 * CE_index_units / pyunits.kg,
        "Pd": 43456.88 * 1e-6 * CE_index_units / pyunits.kg,
        "Pt": 31282.68 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Praseodymium's (Pr) price is not available in USGS, nor its oxides' price.
        "Pr": 500
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.6%. https://www.luciteria.com/elements-for-sale/buy-praseodymium. 1kg price.
        "Rh": 214142.62 * 1e-6 * CE_index_units / pyunits.kg,
        "Rb": 121000 * 1e-6 * CE_index_units / pyunits.kg,
        "Ru": 14998 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Samarium's (Sm) price is not available in USGS, nor its oxides' price.
        "Sm": 140
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.95%. https://www.luciteria.com/elements-for-sale/samarium-metal-999-dendritic. 1kg price.
        "Sc": 153000 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Tantalum's (Ta) price is not available in USGS, but USGS has Ta2O3's price.
        "Ta": 810
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.95%.https://www.luciteria.com/elements-for-sale/buy-tantalum
        "Te": 79.09 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Terbium's (Tb) price is not available in USGS, but USGS has Tb4O7's price.
        "Tb": 2850
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.9%. https://www.luciteria.com/elements-for-sale/terbium-metal-999-pieces. 1kg price.
        # Pure Thulium's (Tm) price is not available in USGS, nor its oxides' price.
        "Tm": 1500
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.95%. https://www.luciteria.com/elements-for-sale/buy-thulium?srsltid=AfmBOoogD1Sf1Nl5FvCGmo3jYTov6EfHWAWUZ8eWkhdSlKb-jKr8csYF
        "Sn": 27.69 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Titanium's (Ti) price is not available in USGS, but USGS has TiO2's price.
        "Ti": 62
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.9%. https://www.luciteria.com/elements-for-sale/buy-titanium
        # Pure Tungsten's (W) price is not available in USGS, but USGS has WO3's price.
        "W": 155
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.95%. https://www.luciteria.com/elements-for-sale/buy-tungsten
        "V": 16.53 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Ytterbium's (Yb) price is not available in USGS, nor its oxides' price.
        "Yb": 375
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.99%. https://www.luciteria.com/elements-for-sale/buy-ytterbium. 1kg price.
        "Y": 33 * 1e-6 * CE_index_units / pyunits.kg,
        "Zn": 3.34 * 1e-6 * CE_index_units / pyunits.kg,
        "Zr": 28 * 1e-6 * CE_index_units / pyunits.kg,
        # oxides. 11 oxides are cited from https://www.usgs.gov/centers/national-minerals-information-center/minerals-yearbook-metals-and-minerals, year 2023, unless otherwise noticed.
        "CeO2": 1 * 1e-6 * CE_index_units / pyunits.kg,  # Purity 99.5%.
        "Dy2O3": 330 * 1e-6 * CE_index_units / pyunits.kg,  # Purity 99.5%.
        "Eu2O3": 27 * 1e-6 * CE_index_units / pyunits.kg,  # Purity 99.99%.
        "La2O3": 1 * 1e-6 * CE_index_units / pyunits.kg,  # Purity 99.5%.
        "Nd2O3": 78 * 1e-6 * CE_index_units / pyunits.kg,  # Purity 99.5%.
        "Sc2O3": 2100 * 1e-6 * CE_index_units / pyunits.kg,  # Purity 99.99%. year 2022.
        "Ta2O5": 170 * 1e-6 * CE_index_units / pyunits.kg,  # Purity 99.9%.
        "Tb4O7": 1298 * 1e-6 * CE_index_units / pyunits.kg,  # Purity 99.9%.
        "TiO2": 1.46 * 1e-6 * CE_index_units / pyunits.kg,  # Purity 99.9%.
        "WO3": 0.26 * 1e-6 * CE_index_units / pyunits.kg,  # Purity 99.9%.
        "Y2O3": 8 * 1e-6 * CE_index_units / pyunits.kg,  # purity 99.9%.
        "Er2O3": 1775.95
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.99%. https://www.msesupplies.com/products/mse-pro-erbium-iii-oxide-er-sub-2-sub-o-sub-3-sub-99-99-4n-powder?variant=31401261629498. 1kg price.
        "Ho2O3": 1995.95
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.99%. https://www.msesupplies.com/products/mse-pro-holmium-iii-oxide-ho-sub-2-sub-o-sub-3-sub-99-99-4n-powder?variant=31405247856698. 1kg price.
        "Gd2O3": 795.95
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.99%. https://www.msesupplies.com/products/mse-pro-gadolinium-iii-oxide-gd-sub-2-sub-o-sub-3-sub-99-999-5n-powder?variant=31402238345274. 1kg price.
        "Lu2O3": 797
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.995%. https://www.msesupplies.com/products/mse-pro-lutetium-iii-oxide-lu-sub-2-sub-o-sub-3-sub-99-995-4n5-powder?variant=31401378644026. 1kg price.
        "Pr6O11": 995.95
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.9%. https://www.msesupplies.com/products/mse-pro-praseodymium-iii-iv-oxide-pr-sub-6-sub-o-sub-11-sub-99-9-3n-powder?variant=31797802729530. 1kg price.
        "Sm2O3": 445.95
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.99%. https://www.msesupplies.com/products/mse-pro-samarium-iii-oxide-sm-sub-2-sub-o-sub-3-sub-99-99-4n-powder?variant=31799234101306. 1kg price.
        "Tm2O3": 1495.95
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.99%. https://www.msesupplies.com/products/mse-pro-thulium-oxide-tm-sub-2-sub-o-sub-3-sub-99-99-4n-powder?variant=31799741677626. 1kg price.
        "Yb2O3": 2319.5
        * 1e-6
        * CE_index_units
        / pyunits.kg,  # Purity 99.99%. https://www.msesupplies.com/products/mse-pro-ytterbium-oxide-yb-sub-2-sub-o-sub-3-sub-powder-99-99-4n-high-purity. 1kg price
    }
    return default_sale_prices


def load_default_resource_prices():
    """
    Dictionary of default prices
    MUSD: the currency units are millions of USD, so its price need a 1e-6 multiplier to get USD
    """
    register_idaes_currency_units()
    CE_index_units = pyunits.MUSD_2023

    default_resource_prices = {
        "power": 0.0804
        * pyunits.USD_2023
        / pyunits.kWh,  # Average industrial electricity rates (2023). https://www.eia.gov/electricity/monthly/epm_table_grapher.php?t=epmt_5_3
        "water": 3.86e-3
        * pyunits.USD_2021
        / pyunits.gallon,  # Average industrial water rates (2021). https://www.osti.gov/servlets/purl/1975260.
        "diesel": 4.214
        * pyunits.USD_2023
        / pyunits.gallon,  # https://www.eia.gov/dnav/pet/pet_pri_gnd_dcus_nus_a.htm. Diesel price annual average in U.S (2023).
        "bioleaching_solution": 0.008 * 1e-6 * CE_index_units / pyunits.L,
        "H2SO4": 128.00
        * pyunits.USD_2023
        / pyunits.tonne,  # Average price of year 2023. https://businessanalytiq.com/procurementanalytics/index/sulfuric-acid-price-index/. Accessed 1/16/2025
        "natural_gas": 4.53
        * 1e-3
        * pyunits.USD_2023
        / pyunits.ft
        ** 3,  # U.S. Annual Industrial price. https://www.eia.gov/dnav/ng/ng_pri_sum_dcu_nus_a.htm
        "polymer": 33.61 * 1e-6 * CE_index_units / pyunits.kg,
        "NAOH": 350.00
        * pyunits.USD_2020
        / pyunits.tonne,  # (price year 2020) https://www.intratec.us/chemical-markets/caustic-soda-price. Accessed 1/16/2025
        "CACO3": 1030.00
        * pyunits.USD_2020
        / pyunits.tonne,  # (price year 2020) https://www.intratec.us/chemical-markets/calcium-carbonate-price. Accessed 1/16/2025
        "coal_calcite": 0.50 * 1e-6 * CE_index_units / pyunits.tonne,
        "HCL": 250.00 * 1e-6 * CE_index_units / pyunits.tonne,
        "oxalic_acid": 1.00 * 1e-6 * CE_index_units / pyunits.kg,
        "ascorbic_acid": 2.00 * 1e-6 * CE_index_units / pyunits.kg,
        "kerosene": 2.699
        * pyunits.USD_2023
        / pyunits.gallon,  # Annual average Kerosene price (2023). https://www.eia.gov/dnav/pet/hist/LeafHandler.ashx?n=PET&s=EER_EPJK_PF4_RGC_DPG&f=A
        "D2EHPA": 15.00
        * pyunits.USD_2023
        / pyunits.kg,  # industry grade. https://kemcore.com/products/d2ehpa-95. Accessed 1/16/2025
        "NA2S": 655.00
        * pyunits.USD_2020
        / pyunits.tonne,  # (price year 2020) https://www.intratec.us/chemical-markets/sodium-sulfides-price. Accessed 1/16/2025
        "nonhazardous_solid_waste": 1.00 * 1e-6 * CE_index_units / pyunits.ton,
        "nonhazardous_precipitate_waste": 5.00 * 1e-6 * CE_index_units / pyunits.ton,
        "dust_and_volatiles": 1.00 * 1e-6 * CE_index_units / pyunits.ton,
    }
    return default_resource_prices
