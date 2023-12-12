from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    SolverFactory,
    Suffix,
    TransformationFactory,
    value,
    Var,
)
from pyomo.core.base.suffix import SuffixFinder

from idaes.core import FlowsheetBlock
from idaes.models.unit_models.gibbs_reactor import GibbsReactor
from idaes.models.properties.activity_coeff_models.methane_combustion_ideal import (
    MethaneParameterBlock as MethaneCombustionParameterBlock,
)
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    InitializationStatus,
)
from idaes.core.util.scaling import (
    get_jacobian,
    extreme_jacobian_columns,
    extreme_jacobian_rows,
    extreme_jacobian_entries,
    jacobian_cond,
)


def autoscale_variables_by_magnitude(blk, overwrite: bool=False, zero_tolerance:float=1e-10):
    """
    Calculate scaling factors for all variables in a model based on their
    current magnitude.

    Args:
        blk - block or model to calculate scaling factors for
        overwrite - whether to overwrite existing scaling factors (default=True)
        zero_tolerance - tolerance for determining when a term is equivalent to zero
            (scaling factor=1)

    Returns:
        Suffix of all scaling factors for model

    """
    # Get scaling suffix
    try:
        sfx = blk.scaling_factor
    except AttributeError:
        # No existing suffix, create one
        sfx = blk.scaling_factor = Suffix(direction=Suffix.EXPORT)

    # Variable scaling
    for v in blk.component_data_objects(Var, descend_into=True):
        if v in sfx and not overwrite:
            # Suffix entry exists and do not overwrite
            continue
        elif v.fixed:
            # Fixed var
            continue

        if v.value is None:
            sf = 1
        else:
            val = abs(value(v))
            if val <= zero_tolerance:
                sf = 1
            else:
                sf = 1/val

        sfx[v] = sf

    return sfx


def autoscale_constraints_by_jacobian_norm(blk, overwrite=False, zero_tolerance=1e-10):
    """
    Calculate scaling factors for all constraints in a model based on the 2-norm of
    the Jacobian matrix, accounting for any variable scaling factors.

    This function is intended to be used on models which already have good variable scaling.

    Args:
        blk - block or model to calculate scaling factors for
        overwrite - whether to overwrite existing scaling factors (default=True)
        zero_tolerance - tolerance for determining when a term is equivalent to zero
            (scaling factor=1)

    Returns:
        Suffix of all scaling factors for model

    """
    # Get scaling suffix
    try:
        sfx = blk.scaling_factor
    except AttributeError:
        # No existing suffix, create one
        sfx = blk.scaling_factor = Suffix(direction=Suffix.EXPORT)

    # Constraint Scaling
    # We want to get any existing scaling
    jac, nlp = get_jacobian(blk, scaled=False)

    eq_con_list = nlp.get_pyomo_equality_constraints()
    var_list = nlp.get_pyomo_variables()

    for c_idx, c in enumerate(eq_con_list):
        if c in sfx and not overwrite:
            # Suffix entry exists and do not overwrite
            continue

        # Get norm of Jacobian row for constraint
        # Need to include variable scaling factors
        n2 = 0
        for v_idx, v_obj in enumerate(var_list):
            # Get variable scaling factor if present
            try:
                v_sf = sfx[v_obj]
            except KeyError:
                v_sf = 1

            n2 += (jac[c_idx, v_idx]/v_sf) ** 2

        # Calculate 2-norm for row
        n = n2 ** 0.5

        if n <= zero_tolerance:
            sf = 1
        else:
            sf = 1 / n

        blk.scaling_factor[c] = sf

    return sfx


# ----------------------------------------------------------------------------------------------------------------------
# Example
def gibbs_scaling_example():
    # Build model
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = MethaneCombustionParameterBlock()

    m.fs.unit = GibbsReactor(
        property_package=m.fs.properties,
        has_heat_transfer=True,
        has_pressure_change=True,
    )

    m.fs.unit.inlet.flow_mol[0].fix(230.0)
    m.fs.unit.inlet.mole_frac_comp[0, "H2"].fix(0.0435)
    m.fs.unit.inlet.mole_frac_comp[0, "N2"].fix(0.6522)
    m.fs.unit.inlet.mole_frac_comp[0, "O2"].fix(0.1739)
    m.fs.unit.inlet.mole_frac_comp[0, "CO2"].fix(1e-5)
    m.fs.unit.inlet.mole_frac_comp[0, "CH4"].fix(0.1304)
    m.fs.unit.inlet.mole_frac_comp[0, "CO"].fix(1e-5)
    m.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix(1e-5)
    m.fs.unit.inlet.mole_frac_comp[0, "NH3"].fix(1e-5)
    m.fs.unit.inlet.temperature[0].fix(1500.0)
    m.fs.unit.inlet.pressure[0].fix(101325.0)

    m.fs.unit.outlet.temperature[0].fix(2844.38)
    m.fs.unit.deltaP.fix(0)

    initializer = BlockTriangularizationInitializer(constraint_tolerance=2e-5)
    initializer.initialize(
        m.fs.unit,
        initial_guesses={
            "control_volume.properties_out[0].pressure": 101325.0,
            "control_volume.properties_out[0].flow_mol": 251.05,
            "control_volume.properties_out[0].mole_frac_comp[CH4]": 1e-5,
            "control_volume.properties_out[0].mole_frac_comp[CO]": 0.0916,
            "control_volume.properties_out[0].mole_frac_comp[CO2]": 0.0281,
            "control_volume.properties_out[0].mole_frac_comp[H2]": 0.1155,
            "control_volume.properties_out[0].mole_frac_comp[H2O]": 0.1633,
            "control_volume.properties_out[0].mole_frac_comp[N2]": 0.59478,
            "control_volume.properties_out[0].mole_frac_comp[NH3]": 1e-5,
            "control_volume.properties_out[0].mole_frac_comp[O2]": 0.0067,
        },
    )

    assert initializer.summary[m.fs.unit]["status"] == InitializationStatus.Ok

    solver = SolverFactory("ipopt")
    solver.solve(m, tee=True)

    print(f"Original Condition No.: {jacobian_cond(m, scaled=False)}")

    autoscale_variables_by_magnitude(m, overwrite=True)

    scaling = TransformationFactory('core.scale_model')
    sm_inter = scaling.create_using(m, rename=False)
    print(f"Intermediate Condition No.: {jacobian_cond(sm_inter, scaled=False)}")

    autoscale_constraints_by_jacobian_norm(m, overwrite=True)

    sm_final = scaling.create_using(m, rename=False)
    print(f"Final Condition No.: {jacobian_cond(sm_final, scaled=False)}")
