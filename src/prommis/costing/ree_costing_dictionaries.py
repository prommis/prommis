#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Python script to read costing components
This script reads the library of costing components (scaled cost, reference
parameters, costing exponents, etc.) from the json files.

"""

# TODO: Missing docstrings
# pylint: disable=missing-function-docstring

__author__ = (
    "Costing Team (B. Paul, A. Fritz, A. Ojo, A. Dasgupta, L. Deng, and M. Zamarripa)"
)
__version__ = "1.0.0"

import json
import os

from pyomo.common.fileutils import this_file_dir
from pyomo.environ import units as pyunits

import idaes.logger as idaeslog
from idaes.core import register_idaes_currency_units

directory = this_file_dir()

_log = idaeslog.getLogger(__name__)


def register_ree_currency_units():
    """
    Define conversion rates for US Dollars based on CEPCI.
    """
    register_idaes_currency_units()
    if (
        "USD_2025" in pyunits._pint_registry  # pylint: disable=protected-access
        and "USD_UKy_2019" in pyunits._pint_registry  # pylint: disable=protected-access
    ):
        # Assume that custom REE plant units have already been registered
        # Log a message and end
        _log.info(
            "Custom REE plant currency units (USD_2025, USD_UKy_2019) "
            "already appear in Pyomo unit registry. Assuming repeated "
            "call of register_ree_currency_units."
        )
    else:
        pyunits.load_definitions_from_strings(
            [
                # from UKy 2019 report
                "USD_UKy_2019 = 500/609.495 * USD_CE500",
                # from https://toweringskills.com/financial-analysis/cost-indices/ as of 9/26/2023
                # from UKy 2023 report
                "USD_2025 = 500/815.59 * USD_CE500",
            ]
        )


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
    that the reference parameter value will be given in. It includes the
    total plant cost (TPC), reference parameter value, and units for that
    value.

    This dictionary is nested with the following structure:
    tech --> ccs --> account --> property name --> property values
    """
    with open(os.path.join(directory, "REE_costing_parameters.json"), "r") as file:
        REE_costing_parameters = json.load(file)
    return REE_costing_parameters


def load_default_sale_prices():
    """
    Dictionary of default prices
    """
    register_idaes_currency_units()

    default_sale_prices = {
        # Pure elements. 33 pure elements prices are cited from year 2023, unless otherwise noticed.
        # https://www.usgs.gov/centers/national-minerals-information-center/
        # minerals-yearbook-metals-and-minerals
        # minimum purity of 99.9% for pure and oxides is assumed, unless mentioned otherwise.
        "Al": 2.78 * pyunits.USD_2023 / pyunits.kg,
        "Sb": 12.1 * pyunits.USD_2023 / pyunits.kg,
        "As": 4.52 * pyunits.USD_2023 / pyunits.kg,
        "Ba": 0.22 * pyunits.USD_2023 / pyunits.kg,
        "Be": 1400 * pyunits.USD_2023 / pyunits.kg,
        "Bi": 8.99 * pyunits.USD_2023 / pyunits.kg,
        # Pure Cerium's (Ce) price is not available in USGS, but has CeO2's price.
        # Purity 99.5%.1kg price. https://www.luciteria.com/elements-for-sale/buy-cerium.
        "Ce": 100 * pyunits.USD_2023 / pyunits.kg,
        "Cs": 91600 * pyunits.USD_2023 / pyunits.kg,
        "Cr": 11.13 * pyunits.USD_2023 / pyunits.kg,
        "Co": 34.13 * pyunits.USD_2023 / pyunits.kg,
        # Pure Dysprosium's (Dy) price is not available in USGS, but has Dy2O3's price.1kg price.
        # Purity 99.99%.https://www.luciteria.com/elements-for-sale/dysprosium-metal-9999-dendritic.
        "Dy": 1900 * pyunits.USD_2023 / pyunits.kg,
        # Pure Erbium's (Er) price is not available in USGS, nor its oxides' price.1kg price.
        # Purity 99.99%. https://www.luciteria.com/elements-for-sale/erbium-metal-9999-dendritic.
        "Er": 39 * pyunits.USD_2023 / pyunits.kg,
        # Pure Europium (Eu) price is not available in USGS, but has Eu2O3's price.1kg price.
        # Purity 99.9%. https://www.luciteria.com/elements-for-sale/buy-europium.
        "Eu": 1550 * pyunits.USD_2023 / pyunits.kg,
        # Fluorspar, alternate Names: Fluorite, Calcium Fluoride.
        "CaF2": 0.296 * pyunits.USD_2023 / pyunits.kg,
        # Pure Gadolinium's (Gd) price is not available in USGS, nor its oxides' price. 1kg price.
        # Purity 99.95%. https://www.luciteria.com/elements-for-sale/gadolinium-metal-9995-dendritic
        "Gd": 850 * pyunits.USD_2023 / pyunits.kg,
        "Ga": 450 * pyunits.USD_2023 / pyunits.kg,
        "Ge": 1392 * pyunits.USD_2023 / pyunits.kg,
        "C": 1.08 * pyunits.USD_2023 / pyunits.kg,
        "Ha": 6150 * pyunits.USD_2023 / pyunits.kg,
        # Pure Holmium's (Ho) price is not available in USGS, nor its oxides' price.1kg price.
        # Purity 99.95%. https://www.luciteria.com/elements-for-sale/holmium-metal-9995-dendritic.
        "Ho": 1600 * pyunits.USD_2023 / pyunits.kg,
        "In": 244 * pyunits.USD_2023 / pyunits.kg,
        "Ir": 150233.37 * pyunits.USD_2023 / pyunits.kg,
        # Pure Lanthanum's (La) price is not available in USGS, but has La2O3's price.1kg price.
        #  https://www.luciteria.com/elements-for-sale/buy-lanthanum. Purity 99.5%.
        "La": 110 * pyunits.USD_2023 / pyunits.kg,
        "Li": 41.3 * pyunits.USD_2023 / pyunits.kg,
        # Pure Lutetium's (Lu) price is not available in USGS, nor its oxides' price.
        # Purity 99.9%. https://www.luciteria.com/elements-for-sale/buy-lutetium. 1kg price.
        "Lu": 3600 * pyunits.USD_2023 / pyunits.kg,
        "Mg": 11.02 * pyunits.USD_2023 / pyunits.kg,
        "Mn": 0.0048 * pyunits.USD_2023 / pyunits.kg,
        # Pure Neodymium's (Nd) price is not available in USGS, but has Nd2O3's price.1kg price.
        # Purity 99.5%. https://www.luciteria.com/elements-for-sale/buy-neodymium-metal.
        "Nd": 425 * pyunits.USD_2023 / pyunits.kg,
        "Ni": 21.495 * pyunits.USD_2023 / pyunits.kg,
        "Nb": 25 * pyunits.USD_2023 / pyunits.kg,
        "Pd": 43456.88 * pyunits.USD_2023 / pyunits.kg,
        "Pt": 31282.68 * pyunits.USD_2023 / pyunits.kg,
        # Pure Praseodymium's (Pr) price is not available in USGS, nor its oxides' price.
        # Purity 99.6%. https://www.luciteria.com/elements-for-sale/buy-praseodymium. 1kg price.
        "Pr": 500 * pyunits.USD_2023 / pyunits.kg,
        "Rh": 214142.62 * pyunits.USD_2023 / pyunits.kg,
        "Rb": 121000 * pyunits.USD_2023 / pyunits.kg,
        "Ru": 14998 * pyunits.USD_2023 / pyunits.kg,
        # Pure Samarium's (Sm) price is not available in USGS, nor its oxides' price.1kg price.
        # Purity 99.95%. https://www.luciteria.com/elements-for-sale/samarium-metal-999-dendritic.
        "Sm": 140 * pyunits.USD_2023 / pyunits.kg,
        "Sc": 153000 * pyunits.USD_2023 / pyunits.kg,
        # Pure Tantalum's (Ta) price is not available in USGS, but has Ta2O3's price.
        # Purity 99.95%.https://www.luciteria.com/elements-for-sale/buy-tantalum
        "Ta": 810 * pyunits.USD_2023 / pyunits.kg,
        "Te": 79.09 * pyunits.USD_2023 / pyunits.kg,
        # Pure Terbium's (Tb) price is not available in USGS, but has Tb4O7's price.1kg price.
        # Purity 99.9%. https://www.luciteria.com/elements-for-sale/terbium-metal-999-pieces.
        "Tb": 2850 * pyunits.USD_2023 / pyunits.kg,
        # Pure Thulium's (Tm) price is not available in USGS, nor its oxides' price.
        # Purity 99.95%. https://www.luciteria.com/elements-for-sale/buy-
        # thulium?srsltid=AfmBOoogD1Sf1Nl5FvCGmo3jYTov6EfHWAWUZ8eWkhdSlKb-jKr8csYF
        "Tm": 1500 * pyunits.USD_2023 / pyunits.kg,
        "Sn": 27.69 * pyunits.USD_2023 / pyunits.kg,
        # Pure Titanium's (Ti) price is not available in USGS, but has TiO2's price.
        # Purity 99.9%. https://www.luciteria.com/elements-for-sale/buy-titanium
        "Ti": 62 * pyunits.USD_2023 / pyunits.kg,
        # Pure Tungsten's (W) price is not available in USGS, but has WO3's price.
        # Purity 99.95%. https://www.luciteria.com/elements-for-sale/buy-tungsten
        "W": 155 * pyunits.USD_2023 / pyunits.kg,
        "V": 16.53 * pyunits.USD_2023 / pyunits.kg,
        # Pure Ytterbium's (Yb) price is not available in USGS, nor its oxides' price.1kg price.
        # Purity 99.99%. https://www.luciteria.com/elements-for-sale/buy-ytterbium.
        "Yb": 375 * pyunits.USD_2023 / pyunits.kg,
        "Y": 33 * pyunits.USD_2023 / pyunits.kg,
        "Zn": 3.34 * pyunits.USD_2023 / pyunits.kg,
        "Zr": 28 * pyunits.USD_2023 / pyunits.kg,
        # oxides. 11 oxides are cited from https://www.usgs.gov/centers/national
        # -minerals-information-center/minerals-yearbook-metals-and-minerals
        # year 2023, unless otherwise noticed.
        "CeO2": 1 * pyunits.USD_2023 / pyunits.kg,  # Purity 99.5%.
        "Dy2O3": 330 * pyunits.USD_2023 / pyunits.kg,  # Purity 99.5%.
        "Eu2O3": 27 * pyunits.USD_2023 / pyunits.kg,  # Purity 99.99%.
        "La2O3": 1 * pyunits.USD_2023 / pyunits.kg,  # Purity 99.5%.
        "Nd2O3": 78 * pyunits.USD_2023 / pyunits.kg,  # Purity 99.5%.
        "Sc2O3": 2100 * pyunits.MUSD_2022 / pyunits.kg,  # Purity 99.99%. year 2022.
        "Ta2O5": 170 * pyunits.USD_2023 / pyunits.kg,  # Purity 99.9%.
        "Tb4O7": 1298 * pyunits.USD_2023 / pyunits.kg,  # Purity 99.9%.
        "TiO2": 1.46 * pyunits.USD_2023 / pyunits.kg,  # Purity 99.9%.
        "WO3": 0.26 * pyunits.USD_2023 / pyunits.kg,  # Purity 99.9%.
        "Y2O3": 8 * pyunits.USD_2023 / pyunits.kg,  # purity 99.9%.
        # Purity 99.99%. https://www.msesupplies.com/products/mse-pro-erbium-iii
        # -oxide-er-sub-2-sub-o-sub-3-sub-99-99-4n-powder?variant=31401261629498. 1kg price.
        "Er2O3": 1775.95 * pyunits.USD_2023 / pyunits.kg,
        # Purity 99.99%. https://www.msesupplies.com/products/mse-pro-holmium-iii-
        # oxide-ho-sub-2-sub-o-sub-3-sub-99-99-4n-powder?variant=31405247856698. 1kg price.
        "Ho2O3": 1995.95 * pyunits.USD_2023 / pyunits.kg,
        # Purity 99.99%. https://www.msesupplies.com/products/mse-pro-gadolinium-
        # iii-oxide-gd-sub-2-sub-o-sub-3-sub-99-999-5n-powder?variant=31402238345274. 1kg price.
        "Gd2O3": 795.95 * pyunits.USD_2023 / pyunits.kg,
        # Purity 99.995%. https://www.msesupplies.com/products/mse-pro-lutetium-
        # iii-oxide-lu-sub-2-sub-o-sub-3-sub-99-995-4n5-powder?variant=31401378644026. 1kg price.
        "Lu2O3": 797 * pyunits.USD_2023 / pyunits.kg,
        # Purity 99.9%. https://www.msesupplies.com/products/mse-pro-praseodymium-
        # iii-iv-oxide-pr-sub-6-sub-o-sub-11-sub-99-9-3n-powder?variant=31797802729530. 1kg price.
        "Pr6O11": 995.95 * pyunits.USD_2023 / pyunits.kg,
        # Purity 99.99%. https://www.msesupplies.com/products/mse-pro-samarium-iii-
        # oxide-sm-sub-2-sub-o-sub-3-sub-99-99-4n-powder?variant=31799234101306. 1kg price.
        "Sm2O3": 445.95 * pyunits.USD_2023 / pyunits.kg,
        # Purity 99.99%. https://www.msesupplies.com/products/mse-pro-thulium-oxide-
        # tm-sub-2-sub-o-sub-3-sub-99-99-4n-powder?variant=31799741677626. 1kg price.
        "Tm2O3": 1495.95 * pyunits.USD_2023 / pyunits.kg,
        # Purity 99.99%. https://www.msesupplies.com/products/mse-pro-ytterbium-
        # oxide-yb-sub-2-sub-o-sub-3-sub-powder-99-99-4n-high-purity. 1kg price
        "Yb2O3": 2319.5 * pyunits.USD_2023 / pyunits.kg,
    }
    return default_sale_prices


def load_default_resource_prices():
    """
    Dictionary of default prices
    """

    default_resource_prices = {
        # Average industrial electricity rates (2023).
        # https://www.eia.gov/electricity/monthly/epm_table_grapher.php?t=epmt_5_3
        "power": 0.0804 * pyunits.USD_2023 / pyunits.kWh,
        # Average industrial water rates (2021). https://www.osti.gov/servlets/purl/1975260.
        "water": 3.86e-3 * pyunits.USD_2021 / pyunits.gallon,
        # https://www.eia.gov/dnav/pet/pet_pri_gnd_dcus_nus_a.htm.
        # Diesel price annual average in U.S (2023).
        "diesel": 4.214 * pyunits.USD_2023 / pyunits.gallon,
        "bioleaching_solution": 0.008 * pyunits.USD_2019 / pyunits.L,
        # Average price of year 2023. https://businessanalytiq.com/procurementanalytics/
        # index/sulfuric-acid-price-index/. Accessed 1/16/2025
        "H2SO4": 128.00 * pyunits.USD_2023 / pyunits.tonne,
        # U.S. Annual Industrial price. https://www.eia.gov/dnav/ng/ng_pri_sum_dcu_nus_a.htm
        "natural_gas": 4.53e-3 * pyunits.USD_2023 / pyunits.ft**3,
        "polymer": 33.61 * pyunits.USD_2019 / pyunits.kg,
        # (price year 2020) https://www.intratec.us/chemical-markets/caustic-soda-price.
        # Accessed 1/16/2025
        "NAOH": 350.00 * pyunits.USD_2020 / pyunits.tonne,
        # (price year 2020) https://www.intratec.us/chemical-markets/calcium-carbonate-price.
        # Accessed 1/16/2025
        "CACO3": 1030.00 * pyunits.USD_2020 / pyunits.tonne,
        "coal_calcite": 0.50 * pyunits.USD_2019 / pyunits.tonne,
        "HCl": 250.00 * pyunits.USD_2019 / pyunits.tonne,
        "oxalic_acid": 1.00 * pyunits.USD_2019 / pyunits.kg,
        "ascorbic_acid": 2.00 * pyunits.USD_2019 / pyunits.kg,
        # Annual average Kerosene price (2023). https://www.eia.gov/dnav/pet/hist/LeafHandler.
        # ashx?n=PET&s=EER_EPJK_PF4_RGC_DPG&f=A
        "kerosene": 2.699 * pyunits.USD_2023 / pyunits.gallon,
        # industry grade. https://kemcore.com/products/d2ehpa-95. Accessed 1/16/2025
        "D2EHPA": 15.00 * pyunits.USD_2023 / pyunits.kg,
        # (price year 2020) https://www.intratec.us/chemical-markets/sodium-sulfides-price.
        # Accessed 1/16/2025
        "NA2S": 655.00 * pyunits.USD_2020 / pyunits.tonne,
        "nonhazardous_solid_waste": 1.00 * pyunits.USD_2019 / pyunits.ton,
        "nonhazardous_precipitate_waste": 5.00 * pyunits.USD_2019 / pyunits.ton,
        "dust_and_volatiles": 1.00 * pyunits.USD_2019 / pyunits.ton,
    }
    return default_resource_prices
