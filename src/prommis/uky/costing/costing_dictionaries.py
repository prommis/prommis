#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
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
    # pyunits.MUSD_2023
    CE_index_units = pyunits.convert(pyunits.MUSD_2023, to_units=pyunits.MUSD_2021)

    default_sale_prices = {
        # Pure elements. 33 pure elements prices are cited from year 2023, unless otherwise noticed.
        # https://www.usgs.gov/centers/national-minerals-information-center/
        # minerals-yearbook-metals-and-minerals
        # minimum purity of 99.9% for pure and oxides is assumed, unless mentioned otherwise.
        "Al": 2.78 * 1e-6 * CE_index_units / pyunits.kg,
        "Sb": 12.1 * 1e-6 * CE_index_units / pyunits.kg,
        "As": 4.52 * 1e-6 * CE_index_units / pyunits.kg,
        "Ba": 0.22 * 1e-6 * CE_index_units / pyunits.kg,
        "Be": 1400 * 1e-6 * CE_index_units / pyunits.kg,
        "Bi": 8.99 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Cerium's (Ce) price is not available in USGS, but has CeO2's price.
        # Purity 99.5%.1kg price. https://www.luciteria.com/elements-for-sale/buy-cerium.
        "Ce": 100 * 1e-6 * CE_index_units / pyunits.kg,
        "Cs": 91600 * 1e-6 * CE_index_units / pyunits.kg,
        "Cr": 11.13 * 1e-6 * CE_index_units / pyunits.kg,
        "Co": 34.13 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Dysprosium's (Dy) price is not available in USGS, but has Dy2O3's price.1kg price.
        # Purity 99.99%.https://www.luciteria.com/elements-for-sale/dysprosium-metal-9999-dendritic.
        "Dy": 1900 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Erbium's (Er) price is not available in USGS, nor its oxides' price.1kg price.
        # Purity 99.99%. https://www.luciteria.com/elements-for-sale/erbium-metal-9999-dendritic.
        "Er": 39 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Europium (Eu) price is not available in USGS, but has Eu2O3's price.1kg price.
        # Purity 99.9%. https://www.luciteria.com/elements-for-sale/buy-europium.
        "Eu": 1550 * 1e-6 * CE_index_units / pyunits.kg,
        # Fluorspar, alternate Names: Fluorite, Calcium Fluoride.
        "CaF2": 0.296 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Gadolinium's (Gd) price is not available in USGS, nor its oxides' price. 1kg price.
        # Purity 99.95%. https://www.luciteria.com/elements-for-sale/gadolinium-metal-9995-dendritic
        "Gd": 850 * 1e-6 * CE_index_units / pyunits.kg,
        "Ga": 450 * 1e-6 * CE_index_units / pyunits.kg,
        "Ge": 1392 * 1e-6 * CE_index_units / pyunits.kg,
        "C": 1.08 * 1e-6 * CE_index_units / pyunits.kg,
        "Ha": 6150 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Holmium's (Ho) price is not available in USGS, nor its oxides' price.1kg price.
        # Purity 99.95%. https://www.luciteria.com/elements-for-sale/holmium-metal-9995-dendritic.
        "Ho": 1600 * 1e-6 * CE_index_units / pyunits.kg,
        "In": 244 * 1e-6 * CE_index_units / pyunits.kg,
        "Ir": 150233.37 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Lanthanum's (La) price is not available in USGS, but has La2O3's price.1kg price.
        #  https://www.luciteria.com/elements-for-sale/buy-lanthanum. Purity 99.5%.
        "La": 110 * 1e-6 * CE_index_units / pyunits.kg,
        "Li": 41.3 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Lutetium's (Lu) price is not available in USGS, nor its oxides' price.
        # Purity 99.9%. https://www.luciteria.com/elements-for-sale/buy-lutetium. 1kg price.
        "Lu": 3600 * 1e-6 * CE_index_units / pyunits.kg,
        "Mg": 11.02 * 1e-6 * CE_index_units / pyunits.kg,
        "Mn": 0.0048 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Neodymium's (Nd) price is not available in USGS, but has Nd2O3's price.1kg price.
        # Purity 99.5%. https://www.luciteria.com/elements-for-sale/buy-neodymium-metal.
        "Nd": 425 * 1e-6 * CE_index_units / pyunits.kg,
        "Ni": 21.495 * 1e-6 * CE_index_units / pyunits.kg,
        "Nb": 25 * 1e-6 * CE_index_units / pyunits.kg,
        "Pd": 43456.88 * 1e-6 * CE_index_units / pyunits.kg,
        "Pt": 31282.68 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Praseodymium's (Pr) price is not available in USGS, nor its oxides' price.
        # Purity 99.6%. https://www.luciteria.com/elements-for-sale/buy-praseodymium. 1kg price.
        "Pr": 500 * 1e-6 * CE_index_units / pyunits.kg,
        "Rh": 214142.62 * 1e-6 * CE_index_units / pyunits.kg,
        "Rb": 121000 * 1e-6 * CE_index_units / pyunits.kg,
        "Ru": 14998 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Samarium's (Sm) price is not available in USGS, nor its oxides' price.1kg price.
        # Purity 99.95%. https://www.luciteria.com/elements-for-sale/samarium-metal-999-dendritic.
        "Sm": 140 * 1e-6 * CE_index_units / pyunits.kg,
        "Sc": 153000 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Tantalum's (Ta) price is not available in USGS, but has Ta2O3's price.
        # Purity 99.95%.https://www.luciteria.com/elements-for-sale/buy-tantalum
        "Ta": 810 * 1e-6 * CE_index_units / pyunits.kg,
        "Te": 79.09 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Terbium's (Tb) price is not available in USGS, but has Tb4O7's price.1kg price.
        # Purity 99.9%. https://www.luciteria.com/elements-for-sale/terbium-metal-999-pieces.
        "Tb": 2850 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Thulium's (Tm) price is not available in USGS, nor its oxides' price.
        # Purity 99.95%. https://www.luciteria.com/elements-for-sale/buy-
        # thulium?srsltid=AfmBOoogD1Sf1Nl5FvCGmo3jYTov6EfHWAWUZ8eWkhdSlKb-jKr8csYF
        "Tm": 1500 * 1e-6 * CE_index_units / pyunits.kg,
        "Sn": 27.69 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Titanium's (Ti) price is not available in USGS, but has TiO2's price.
        # Purity 99.9%. https://www.luciteria.com/elements-for-sale/buy-titanium
        "Ti": 62 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Tungsten's (W) price is not available in USGS, but has WO3's price.
        # Purity 99.95%. https://www.luciteria.com/elements-for-sale/buy-tungsten
        "W": 155 * 1e-6 * CE_index_units / pyunits.kg,
        "V": 16.53 * 1e-6 * CE_index_units / pyunits.kg,
        # Pure Ytterbium's (Yb) price is not available in USGS, nor its oxides' price.1kg price.
        # Purity 99.99%. https://www.luciteria.com/elements-for-sale/buy-ytterbium.
        "Yb": 375 * 1e-6 * CE_index_units / pyunits.kg,
        "Y": 33 * 1e-6 * CE_index_units / pyunits.kg,
        "Zn": 3.34 * 1e-6 * CE_index_units / pyunits.kg,
        "Zr": 28 * 1e-6 * CE_index_units / pyunits.kg,
        # oxides. 11 oxides are cited from https://www.usgs.gov/centers/national
        # -minerals-information-center/minerals-yearbook-metals-and-minerals
        # year 2023, unless otherwise noticed.
        "CeO2": 1 * 1e-6 * CE_index_units / pyunits.kg,  # Purity 99.5%.
        "Dy2O3": 330 * 1e-6 * CE_index_units / pyunits.kg,  # Purity 99.5%.
        "Eu2O3": 27 * 1e-6 * CE_index_units / pyunits.kg,  # Purity 99.99%.
        "La2O3": 1 * 1e-6 * CE_index_units / pyunits.kg,  # Purity 99.5%.
        "Nd2O3": 78 * 1e-6 * CE_index_units / pyunits.kg,  # Purity 99.5%.
        "Sc2O3": 1e-6
        * pyunits.convert(
            2100 * pyunits.MUSD_2022 / pyunits.kg,
            to_units=pyunits.MUSD_2021 / pyunits.kg,
        ),  # Purity 99.99%. year 2022.
        "Ta2O5": 170 * 1e-6 * CE_index_units / pyunits.kg,  # Purity 99.9%.
        "Tb4O7": 1298 * 1e-6 * CE_index_units / pyunits.kg,  # Purity 99.9%.
        "TiO2": 1.46 * 1e-6 * CE_index_units / pyunits.kg,  # Purity 99.9%.
        "WO3": 0.26 * 1e-6 * CE_index_units / pyunits.kg,  # Purity 99.9%.
        "Y2O3": 8 * 1e-6 * CE_index_units / pyunits.kg,  # purity 99.9%.
        # Purity 99.99%. https://www.msesupplies.com/products/mse-pro-erbium-iii
        # -oxide-er-sub-2-sub-o-sub-3-sub-99-99-4n-powder?variant=31401261629498. 1kg price.
        "Er2O3": 1775.95 * 1e-6 * CE_index_units / pyunits.kg,
        # Purity 99.99%. https://www.msesupplies.com/products/mse-pro-holmium-iii-
        # oxide-ho-sub-2-sub-o-sub-3-sub-99-99-4n-powder?variant=31405247856698. 1kg price.
        "Ho2O3": 1995.95 * 1e-6 * CE_index_units / pyunits.kg,
        # Purity 99.99%. https://www.msesupplies.com/products/mse-pro-gadolinium-
        # iii-oxide-gd-sub-2-sub-o-sub-3-sub-99-999-5n-powder?variant=31402238345274. 1kg price.
        "Gd2O3": 795.95 * 1e-6 * CE_index_units / pyunits.kg,
        # Purity 99.995%. https://www.msesupplies.com/products/mse-pro-lutetium-
        # iii-oxide-lu-sub-2-sub-o-sub-3-sub-99-995-4n5-powder?variant=31401378644026. 1kg price.
        "Lu2O3": 797 * 1e-6 * CE_index_units / pyunits.kg,
        # Purity 99.9%. https://www.msesupplies.com/products/mse-pro-praseodymium-
        # iii-iv-oxide-pr-sub-6-sub-o-sub-11-sub-99-9-3n-powder?variant=31797802729530. 1kg price.
        "Pr6O11": 995.95 * 1e-6 * CE_index_units / pyunits.kg,
        # Purity 99.99%. https://www.msesupplies.com/products/mse-pro-samarium-iii-
        # oxide-sm-sub-2-sub-o-sub-3-sub-99-99-4n-powder?variant=31799234101306. 1kg price.
        "Sm2O3": 445.95 * 1e-6 * CE_index_units / pyunits.kg,
        # Purity 99.99%. https://www.msesupplies.com/products/mse-pro-thulium-oxide-
        # tm-sub-2-sub-o-sub-3-sub-99-99-4n-powder?variant=31799741677626. 1kg price.
        "Tm2O3": 1495.95 * 1e-6 * CE_index_units / pyunits.kg,
        # Purity 99.99%. https://www.msesupplies.com/products/mse-pro-ytterbium-
        # oxide-yb-sub-2-sub-o-sub-3-sub-powder-99-99-4n-high-purity. 1kg price
        "Yb2O3": 2319.5 * 1e-6 * CE_index_units / pyunits.kg,
    }
    return default_sale_prices


def convert_to_usd_2021(value_in_original_units, original_unit):
    """
    Converts a given monetary value from its original year to 2021 USD.

    Args:
        value_in_original_units (float): The numerical value of the cost.
        original_unit (pyunits object): The unit of the value, including the year.

    Returns:
        Converted value with 2021 USD as the unit.
    """
    return pyunits.convert(
        value_in_original_units * original_unit, to_units=pyunits.USD_2021
    )


def load_default_resource_prices():
    """
    Dictionary of default prices
    MUSD: the currency units are millions of USD, so its price need a 1e-6 multiplier to get USD
    """
    register_idaes_currency_units()

    # USD_2019 was the price referenced from Uky report.
    CE_index_units = convert_to_usd_2021(1, pyunits.USD_2019)

    default_resource_prices = {
        # Average industrial electricity rates (2023).
        # https://www.eia.gov/electricity/monthly/epm_table_grapher.php?t=epmt_5_3
        "power": convert_to_usd_2021(0.0804, pyunits.USD_2023) / pyunits.kWh,
        # Average industrial water rates (2021). https://www.osti.gov/servlets/purl/1975260.
        "water": 3.86e-3 * pyunits.USD_2021 / pyunits.gallon,
        # https://www.eia.gov/dnav/pet/pet_pri_gnd_dcus_nus_a.htm.
        # Diesel price annual average in U.S (2023).
        "diesel": convert_to_usd_2021(4.214, pyunits.USD_2023) / pyunits.gallon,
        "bioleaching_solution": 0.008 * 1e-6 * CE_index_units / pyunits.L,
        # Average price of year 2023. https://businessanalytiq.com/procurementanalytics/
        # index/sulfuric-acid-price-index/. Accessed 1/16/2025
        "H2SO4": convert_to_usd_2021(128.00, pyunits.USD_2023) / pyunits.tonne,
        # U.S. Annual Industrial price. https://www.eia.gov/dnav/ng/ng_pri_sum_dcu_nus_a.htm
        "natural_gas": convert_to_usd_2021(4.53e-3, pyunits.USD_2023) / pyunits.ft**3,
        "polymer": 33.61 * 1e-6 * CE_index_units / pyunits.kg,
        # (price year 2020) https://www.intratec.us/chemical-markets/caustic-soda-price.
        # Accessed 1/16/2025
        "NAOH": convert_to_usd_2021(350.00, pyunits.USD_2020) / pyunits.tonne,
        # (price year 2020) https://www.intratec.us/chemical-markets/calcium-carbonate-price.
        # Accessed 1/16/2025
        "CACO3": convert_to_usd_2021(1030.00, pyunits.USD_2020) / pyunits.tonne,
        "coal_calcite": 0.50 * 1e-6 * CE_index_units / pyunits.tonne,
        "HCL": 250.00 * 1e-6 * CE_index_units / pyunits.tonne,
        "oxalic_acid": 1.00 * 1e-6 * CE_index_units / pyunits.kg,
        "ascorbic_acid": 2.00 * 1e-6 * CE_index_units / pyunits.kg,
        # Annual average Kerosene price (2023). https://www.eia.gov/dnav/pet/hist/LeafHandler.
        # ashx?n=PET&s=EER_EPJK_PF4_RGC_DPG&f=A
        "kerosene": convert_to_usd_2021(2.699, pyunits.USD_2023) / pyunits.gallon,
        # industry grade. https://kemcore.com/products/d2ehpa-95. Accessed 1/16/2025
        "D2EHPA": convert_to_usd_2021(15.00, pyunits.USD_2023) / pyunits.kg,
        # (price year 2020) https://www.intratec.us/chemical-markets/sodium-sulfides-price.
        # Accessed 1/16/2025
        "NA2S": convert_to_usd_2021(655.00, pyunits.USD_2020) / pyunits.tonne,
        "nonhazardous_solid_waste": 1.00 * 1e-6 * CE_index_units / pyunits.ton,
        "nonhazardous_precipitate_waste": 5.00 * 1e-6 * CE_index_units / pyunits.ton,
        "dust_and_volatiles": 1.00 * 1e-6 * CE_index_units / pyunits.ton,
    }
    return default_resource_prices
