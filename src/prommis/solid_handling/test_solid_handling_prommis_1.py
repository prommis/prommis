from solid_handling_prommis import (
    CrushAndBreakageUnit,
)  # Replace your_module_name with the actual name of your Python file
from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)

# Create a concrete model
model = ConcreteModel()
model.fs = FlowsheetBlock(default={"dynamic": False})
model.fs.properties = BTXParameterBlock(
    default={"valid_phase": ("Liq", "Vap"), "activity_coeff_model": "Ideal"}
)

# Add CrushAndBreakageUnit to the flowsheet
model.fs.unit = CrushAndBreakageUnit(default={"property_package": model.fs.properties})

# Optional: Set up model inputs, constraints, etc., here

# Solve the model
# from pyomo.environ import SolverFactory
# solver = SolverFactory('ipopt')  # Example solver
# solver.solve(model.fs.unit, tee=True)

# Print results or debug
print("BWI:", model.fs.unit.BWI.value)
