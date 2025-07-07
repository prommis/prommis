# import pyomo.environ as pyo
# from pyomo.environ import units as pyunits

# from idaes.core import register_idaes_currency_units

# register_idaes_currency_units()

# cost_old = 500 * pyunits.USD_2018

# cost_new = pyunits.convert(cost_old, to_units=pyunits.USD_2020)

# print(cost_old)
# print(pyo.value(cost_new))

### Define parameters input
# Define the final year of the plant's lifetime.
plant_start = 2023
plant_lifetime = 15
plant_end = plant_start + plant_lifetime - 1
print(plant_end)
# Define a list to hold the custom yearly costing units that need to be defined.
yearly_cost_units = []


## Define custom units for costing. Need a new costing unit for each year of the plant's lifetime.
for t in range(plant_start, plant_end + 1):
    year_str = str(t)
    unit = "MUSD_" + year_str
    print(unit + " = [currency]")
