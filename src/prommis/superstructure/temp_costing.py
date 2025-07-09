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


for opt in [(1, 1), (1, 2), (1, 3)]:
    new_opt = opt + ("test",)
    print(new_opt)


