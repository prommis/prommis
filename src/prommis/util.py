#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Temporary home for assert_solution_equivalent until it gets moved to IDAES
"""
from math import ceil
import textwrap
import pytest

from pyomo.common.collections import ComponentSet
from pyomo.environ import Expression, log10, Reference, value, Var
from pyomo.dae import DerivativeVar
from pyomo.dae.flatten import flatten_dae_components, slice_component_along_sets

import idaes.logger as idaeslog
from idaes.core.scaling.util import get_scaling_factor, set_scaling_factor

_log = idaeslog.getLogger(__name__)


def assert_solution_equivalent(blk, expected_results):
    """
    Method to iterate through a structured dictionary of variables/expressions, values,
    and indices to determine whether the variables have their expected values within the
    specified tolerances. The variables that do not have their expected values are
    collected and displayed for the user, and an AssertionError is raised.

    This method is better than just writing a long series of assert statements because
    it shows *all* variables/expressions that do not have the expected values in a single
    report, rather than having to change the value in one assert statement, run the test
    again, change the next value, run the test again, etc.

    Args:
        blk: Pyomo block on which variables/expressions being tested are located
        expected_results: Dictionary of the form:
        {
            indexed_var_name: {
                index_1: (value, rel_tol, abs_tol),
                index_2: (value, rel_tol, abs_tol),
                ...
            }
            unindexed_var_name: {
                # Unindexed vars pass None as the index
                None: (value, rel_tol, abs_tol)
            }
            ...
        }
    """

    n_failures = 0
    failures = []

    for name, expected_values_dict in expected_results.items():
        recorded_var = False
        obj = blk.find_component(name)
        if obj is None:
            blk_name = blk.name
            # Pyomo ConcreteModels are named "unknown" by default
            # but seeing "unknown" show up in an error message is confusing.
            if blk_name == "unknown":
                blk_name = "model"
            failure_msg = f"  - Could not find object {name} on {blk_name}\n"
            failures.append(failure_msg)
            continue

        obj_type = None
        if isinstance(obj, Var):
            obj_type = "Variable"
        elif isinstance(obj, Expression):
            obj_type = "Expression"
        else:
            failure_msg = f"  - Error: object {name} is not a Var or Expression\n"
            failures.append(failure_msg)
            continue

        for index, (expected_value, rel, abs) in expected_values_dict.items():
            absent_index = False
            is_close = False
            if index is None:
                component_to_test = obj
            else:
                if index in obj:
                    component_to_test = obj[index]
                else:
                    absent_index = True
            if not absent_index:
                actual_value = value(component_to_test)

                # Determine if the values are approximately equal
                if actual_value == pytest.approx(expected_value, rel=rel, abs=abs):
                    is_close = True
            if (absent_index or not is_close) and not recorded_var:
                failures.append(f"  - {obj_type}: {name}")
                recorded_var = True
            if absent_index:
                failure_msg = f"    Index:    {index} is absent"
                failures.append(failure_msg)
                n_failures += 1
                continue

            # If the comparison fails, record the details
            if not is_close:
                if rel is not None:
                    n_sig_figs = ceil(-log10(rel)) + 1
                    format_spec = "." + str(n_sig_figs) + "e"
                elif abs is not None:
                    n_sig_figs = ceil(-log10(abs)) + 1
                    format_spec = "." + str(n_sig_figs) + "f"
                else:
                    format_spec = ".7e"
                failure_msg = (
                    f"    Index:    {index}\n"
                    f"    Expected: {expected_value:{format_spec}}\n"
                    f"    Actual:   {actual_value:{format_spec}}"
                )
                failures.append(failure_msg)
                n_failures += 1

        if recorded_var:
            # Extra space between variables
            failures[-1] = failures[-1] + "\n"

    # --- Final Assertion and Report Generation ---
    if len(failures) > 0:
        # Construct the final report header
        report_header = textwrap.dedent(
            f"""
        =========================== Test Value Mismatches ============================
        Found {n_failures} mismatch(es) between expected and actual model values.
        Please review the values below and update the test suite if necessary.
        ==============================================================================
        """
        )

        # Combine the header and all failure messages
        full_report = report_header + "\n\n" + "\n\n".join(failures)

        # Raise a single AssertionError with the complete report
        raise AssertionError(full_report)


def copy_first_steady_state(blk, time):
    """
    Function that propagates initial steady state guess to future time points.
    This function is used to initialize all the time discrete variables to the
    initial steady state value.
    """

    _, time_vars = flatten_dae_components(blk, time, Var, active=True)
    # Copy initial conditions forward
    for var in time_vars:
        for t in time:
            if t == time.first():
                continue
            else:
                var[t].value = var[time.first()].value


def scale_time_discretization_equations(blk, time_set, time_scaling_factor):
    """
    Scales time discretization equations generated via a Pyomo discretization
    transformation. Also scales continuity equations for collocation methods
    of discretization that require them.

    Args:
        blk: Block whose time discretization equations are being scaled
        time_set: Time set object. For an IDAES flowsheet object fs, this is fs.time.
        time_scaling_factor: Scaling factor to use for time

    Returns:
        None
    """

    tname = time_set.local_name

    # Copy and pasted from solvers.petsc.find_discretization_equations then modified
    for var in blk.component_objects(Var):
        if isinstance(var, DerivativeVar):
            cont_set_set = ComponentSet(var.get_continuousset_list())
            if time_set in cont_set_set:
                if len(cont_set_set) > 1:
                    _log.warning(
                        "IDAES presently does not support automatically scaling discretization equations for "
                        f"second order or higher derivatives like {var.name} that are differentiated at least once with "
                        "respect to time. Please scale the corresponding discretization equation yourself."
                    )
                    continue
                state_var = var.get_state_var()
                parent_block = var.parent_block()

                disc_eq = getattr(parent_block, var.local_name + "_disc_eq")
                # Look for continuity equation, which exists only for collocation with certain sets of polynomials
                try:
                    cont_eq = getattr(
                        parent_block, state_var.local_name + "_" + tname + "_cont_eq"
                    )
                except AttributeError:
                    cont_eq = None

                deriv_dict = dict(
                    (key, Reference(slc))
                    for key, slc in slice_component_along_sets(var, (time_set,))
                )
                state_dict = dict(
                    (key, Reference(slc))
                    for key, slc in slice_component_along_sets(state_var, (time_set,))
                )
                disc_dict = dict(
                    (key, Reference(slc))
                    for key, slc in slice_component_along_sets(disc_eq, (time_set,))
                )
                if cont_eq is not None:
                    cont_dict = dict(
                        (key, Reference(slc))
                        for key, slc in slice_component_along_sets(cont_eq, (time_set,))
                    )
                for key, deriv in deriv_dict.items():
                    state = state_dict[key]
                    disc = disc_dict[key]
                    if cont_eq is not None:
                        cont = cont_dict[key]
                    for t in time_set:
                        s_state = get_scaling_factor(state[t], default=1, warning=True)
                        set_scaling_factor(
                            deriv[t], s_state / time_scaling_factor, overwrite=False
                        )
                        s_deriv = get_scaling_factor(deriv[t])
                        # Check time index to decide what constraints to scale
                        if cont_eq is None:
                            if t == time_set.first() or t == time_set.last():
                                try:
                                    set_scaling_factor(
                                        disc[t], s_deriv, overwrite=False
                                    )
                                except KeyError:
                                    # Discretization and continuity equations may or may not exist at the first or last time
                                    # points depending on the method. Backwards skips first, forwards skips last, central skips
                                    # both (which means the user needs to provide additional equations)
                                    pass
                            else:
                                set_scaling_factor(disc[t], s_deriv, overwrite=False)
                        else:
                            # Lagrange-Legendre is a pain, because it has continuity equations on the edges of finite
                            # instead of discretization equations, but no intermediate continuity equations, so we have
                            # to look for both at every timepoint
                            try:
                                set_scaling_factor(disc[t], s_deriv, overwrite=False)
                            except KeyError:
                                if t != time_set.first():
                                    set_scaling_factor(
                                        # pylint: disable-next=possibly-used-before-assignment
                                        cont[t],
                                        s_state,
                                        overwrite=False,
                                    )
