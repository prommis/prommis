#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Crusher configuration validation and numerical input preparation."""

__author__ = "Daison Yancy Caballero"

import math
import operator
from collections.abc import Mapping

import numpy

import idaes.logger as idaeslog
from idaes.core.util.exceptions import ConfigurationError

from prommis.comminution.core.config_utils import _cfg_float
from prommis.comminution.core.psd_math import finest_attainable_size
from prommis.comminution.functions import distributions
from prommis.comminution.functions.breakage import breakage_matrix, ecs_from_table
from prommis.comminution.functions.selection import selection_value
from prommis.comminution.unit_models._crusher.shared import (
    SelectionFunction,
    BreakageFunction,
    RecycleMode,
    CrusherEquipment,
    CrusherStage,
    CrusherPowerLaw,
    ProductPath,
    PSDMethod,
    PreparedCrusherData,
    WORK_INDEX_BOUNDS,
    N_STRESS_EVENTS_MAX,
    _TABULAR_SUM_TOL,
    _MIN_SHAPE_SIZE,
    _SHAPE_SIZE_LIMIT_RULE,
    _shape_size_upper_bound,
    validate_recycle_denominators,
)

_log = idaeslog.getLogger("prommis.comminution.unit_models.crusher")


# Equipment presets set the stage and supply the selection and breakage
# functions the user has not chosen.
CRUSHER_EQUIPMENT_TO_STAGE = {
    CrusherEquipment.jaw: CrusherStage.primary,
    CrusherEquipment.gyratory: CrusherStage.primary,
    CrusherEquipment.cone: CrusherStage.secondary,
    CrusherEquipment.roll1: CrusherStage.secondary,
    CrusherEquipment.roll2: CrusherStage.tertiary,
    CrusherEquipment.short_head_cone: CrusherStage.tertiary,
    CrusherEquipment.hammer_mill: CrusherStage.tertiary,
}


# Equipment-specific (selection, breakage) function defaults.
# None keeps the current function and warns if the user did not set one.
CRUSHER_EQUIPMENT_FUNCTIONS = {
    CrusherEquipment.jaw: (SelectionFunction.whiten, BreakageFunction.vogel),
    CrusherEquipment.gyratory: (None, None),
    CrusherEquipment.cone: (SelectionFunction.king, BreakageFunction.vogel),
    CrusherEquipment.short_head_cone: (
        SelectionFunction.king,
        BreakageFunction.vogel,
    ),
    CrusherEquipment.roll1: (None, None),
    CrusherEquipment.roll2: (None, None),
    CrusherEquipment.hammer_mill: (
        SelectionFunction.vogel_peukert,
        BreakageFunction.vogel,
    ),
}


_SELECTION_REQUIRED = {
    SelectionFunction.whiten: ["K1", "K2"],
    SelectionFunction.king: ["CSS", "alpha1", "alpha2", "n"],
    SelectionFunction.austin: ["S1", "d1", "a"],
    SelectionFunction.vogel_peukert: ["f_mat", "xw_min", "v"],
}


_BREAKAGE_REQUIRED = {
    BreakageFunction.luckie_austin: ["phi", "gamma", "beta"],
    BreakageFunction.reid_stewart: ["phi", "gamma", "beta"],
    BreakageFunction.vogel: ["q", "dprime"],
    BreakageFunction.t10_appearance: ["appearance_table", "t10"],
}


# Options owned by each PSD method; setting one for another method is an error.
_METHOD_OPTIONS = {
    PSDMethod.tabular: ("tabular_psd", "tabular_psd_by_comp"),
    PSDMethod.distribution_function: (
        "shape_exponent",
        "target_product_p80",
        "distribution_shape",
    ),
    PSDMethod.selection_breakage: (
        "selection_params",
        "selection_params_by_comp",
        "breakage_params",
        "breakage_params_by_comp",
        "selection_function",
        "breakage_function",
    ),
}


# Machine settings must match across components in selection_params_by_comp.
_SELECTION_EQUIPMENT_KEY = {
    SelectionFunction.king: "CSS",
    SelectionFunction.vogel_peukert: "v",
}


# Each user breakage-matrix column must sum to 1 within this tolerance.
_BREAKAGE_COLUMN_SUM_TOL = 1e-9
# Avoid underflow or overflow in P80 ratios at very small shape exponents.
_MIN_SHAPE_EXPONENT = 1e-3


def validate_and_copy_configuration(config):
    """Validate options, apply equipment defaults, and return ``n_stress_events``.

    Copy ``pendulum_power_params`` because the scaler reads it after build.
    The other nested inputs are used only during build.
    """
    cfg = config
    n_events = _validate_n_stress_events(cfg.n_stress_events)
    # Capture explicit choices before equipment defaults change _userSet.
    user_set = {
        "selection_function": cfg.get("selection_function")._userSet,
        "breakage_function": cfg.get("breakage_function")._userSet,
        "distribution_shape": cfg.get("distribution_shape")._userSet,
    }
    # Use the same work-index bounds for configuration and the Var.
    _cfg_float(
        cfg.bond_work_index,
        "bond_work_index",
        lo=WORK_INDEX_BOUNDS[0],
        hi=WORK_INDEX_BOUNDS[1],
    )
    # Resolve equipment defaults before validating dependent option combinations.
    _resolve_equipment_softbind(config)
    _validate_pendulum_config(
        config, cfg.get("bond_work_index")._userSet, user_set["breakage_function"]
    )
    method = cfg.psd_method
    n_int = len(cfg.property_package.size_interval_set)
    if method == PSDMethod.selection_breakage:
        _reject_foreign_psd_options(config, method, user_set)
        _validate_selection_breakage_config(config, n_events)
        return n_events
    if cfg.get("recycle")._userSet:
        raise ConfigurationError(
            "'recycle' is consumed only by psd_method='selection_breakage'; "
            f"do not set it for psd_method='{method}'."
        )
    _reject_unused_nse(n_events, f"psd_method='{method}'")
    _reject_foreign_psd_options(config, method, user_set)
    if method == PSDMethod.tabular:
        _validate_tabular_config(config, n_int)
    if method == PSDMethod.distribution_function:
        _validate_distribution_config(config, n_int)
    return n_events


def _validate_distribution_config(config, n_int):
    """Validate the mesh size, ``shape_exponent`` and any ``target_product_p80``."""
    cfg = config
    # One interval leaves no shape after normalization and row closure.
    if n_int < 2:
        spec = "size" if cfg.target_product_p80 is not None else "power"
        raise ConfigurationError(
            f"distribution_function {spec} specification requires at "
            f"least two size intervals (got N = {n_int})."
        )
    if cfg.shape_exponent is None:
        raise ConfigurationError(
            "psd_method='distribution_function' requires 'shape_exponent'."
        )
    n_val = _cfg_float(
        cfg.shape_exponent,
        "distribution_function 'shape_exponent'",
        lo=_MIN_SHAPE_EXPONENT,
    )
    if cfg.target_product_p80 is not None:
        _check_target_product_p80(config, cfg.target_product_p80, n_val)


def _check_target_product_p80(config, target, shape_exponent):
    """Require a positive ``target_product_p80`` whose curve fits the mesh."""
    cfg = config
    target = _cfg_float(target, "target_product_p80")
    # Validate P80 at build time; initialization rechecks the mutable Param.
    pp = config.property_package
    edges = list(pp.size_edges_m)
    floor = finest_attainable_size(edges)
    if not floor <= target <= edges[-1]:
        raise ConfigurationError(
            f"target_product_p80 ({target:.6g} m) must lie within the "
            f"attainable mesh band [{floor:.6g}, {edges[-1]:.6g}] m "
            "(from 80 percent through the finest interval to the mesh "
            "top)."
        )
    shape_size = distributions.shape_size_from_p80(
        cfg.distribution_shape, target, shape_exponent
    )
    upper = _shape_size_upper_bound(cfg.distribution_shape, edges[-1], shape_exponent)
    if not _MIN_SHAPE_SIZE <= shape_size <= upper:
        rule = (
            f"; {_SHAPE_SIZE_LIMIT_RULE[cfg.distribution_shape]}"
            if shape_size > upper
            else ""
        )
        raise ConfigurationError(
            f"target_product_p80 = {target!r} m implies a "
            f"'{cfg.distribution_shape}' characteristic size {shape_size!r} m outside "
            f"[{_MIN_SHAPE_SIZE}, {upper!r}] m on this mesh{rule}."
        )


def _validate_pendulum_config(config, bond_work_index_set, breakage_function_set):
    """Require Whiten recycle, t10 breakage, and pendulum power parameters.

    Reject an explicitly set bond work index for pendulum power; other
    power laws reject ``pendulum_power_params``.
    """
    cfg = config
    if cfg.power_law != CrusherPowerLaw.pendulum:
        if cfg.pendulum_power_params is not None:
            raise ConfigurationError(
                "pendulum_power_params is consumed only by power_law='pendulum'."
            )
        return
    if (
        cfg.psd_method != PSDMethod.selection_breakage
        or cfg.recycle != RecycleMode.classification_recycle
        or cfg.breakage_function != BreakageFunction.t10_appearance
    ):
        hint = ""
        if (
            cfg.psd_method == PSDMethod.selection_breakage
            and cfg.crusher_equipment is not None
            and not breakage_function_set
            and CRUSHER_EQUIPMENT_FUNCTIONS[cfg.crusher_equipment][1] is not None
        ):
            hint = (
                f" breakage_function='{cfg.breakage_function}' came from "
                f"crusher_equipment='{cfg.crusher_equipment}'; set "
                "breakage_function='t10_appearance' explicitly."
            )
        raise ConfigurationError(
            "power_law='pendulum' requires psd_method='selection_breakage', "
            "recycle='classification_recycle' and "
            "breakage_function='t10_appearance'." + hint
        )
    if bond_work_index_set:
        raise ConfigurationError(
            "power_law='pendulum' does not use bond_work_index; do not set it."
        )
    params, keys = cfg.pendulum_power_params, ["no_load_power", "power_factor"]
    if not isinstance(params, Mapping) or set(params) != set(keys):
        raise ConfigurationError(
            "power_law='pendulum' requires pendulum_power_params to be a "
            f"mapping with exactly the keys {keys}."
        )
    for key in keys:
        _cfg_float(params[key], f"pendulum_power_params '{key}'", lo=0.0)
    # The scaler reads these after build, so the unit keeps its own copy.
    cfg.pendulum_power_params = dict(params)


def _validate_n_stress_events(val):
    """Return an integer in ``[1, N_STRESS_EVENTS_MAX]``, rejecting booleans."""
    if isinstance(val, (bool, numpy.bool_)):
        raise ConfigurationError("n_stress_events must be an integer, not a boolean.")
    try:
        ival = operator.index(val)
    except TypeError as exc:
        raise ConfigurationError(
            f"n_stress_events must be an integer in [1, {N_STRESS_EVENTS_MAX}] "
            f"(got {val!r})."
        ) from exc
    if not 1 <= ival <= N_STRESS_EVENTS_MAX:
        raise ConfigurationError(
            f"n_stress_events must be an integer in [1, {N_STRESS_EVENTS_MAX}] "
            f"(got {val!r})."
        )
    return ival


def _resolve_equipment_softbind(config):
    """Check equipment/stage consistency and fill unset stage or functions."""
    cfg = config
    eq = cfg.crusher_equipment
    stage = cfg.crusher_stage
    if eq is None and stage is None:
        return
    if eq is not None and stage is not None:
        if CRUSHER_EQUIPMENT_TO_STAGE[eq] != stage:
            raise ConfigurationError(
                f"crusher_stage='{stage}' is inconsistent with "
                f"crusher_equipment='{eq}' (which maps to stage "
                f"'{CRUSHER_EQUIPMENT_TO_STAGE[eq]}')."
            )
    elif eq is not None:
        cfg.crusher_stage = CRUSHER_EQUIPMENT_TO_STAGE[eq]
    if eq is None or cfg.psd_method != PSDMethod.selection_breakage:
        return
    for kind, default in zip(
        ("selection", "breakage"), CRUSHER_EQUIPMENT_FUNCTIONS[eq]
    ):
        option = f"{kind}_function"
        if cfg.get(option)._userSet:
            continue
        if default is None:
            _log.warning(
                "crusher_equipment='%s' has no equipment-specific %s "
                "default; using %s='%s'. Set it explicitly "
                "to override this default.",
                eq,
                kind,
                option,
                cfg[option],
            )
        else:
            cfg[option] = default


def _reject_foreign_psd_options(config, method, user_set):
    """Reject explicitly supplied options unused by the selected PSD method.

    ``user_set`` maps option names to whether the user set them; an option
    missing from it counts as set when its value is not ``None``.
    """
    cfg = config
    unused_options = sorted(
        name
        for owner, names in _METHOD_OPTIONS.items()
        if owner != method
        for name in names
        if (user_set[name] if name in user_set else cfg[name] is not None)
    )
    if unused_options:
        raise ConfigurationError(
            f"psd_method='{method}' does not use these options: {unused_options}. "
            "Remove them from the configuration."
        )


def _reject_unused_nse(n_events, context):
    """Reject ``n_stress_events != 1`` on paths that do not use it."""
    if n_events != 1:
        raise ConfigurationError(
            "n_stress_events can differ from 1 only with "
            "psd_method='selection_breakage' and recycle='single_pass'. "
            f"For {context}, omit it or set it to 1."
        )


def _validate_tabular_config(config, n_int):
    """Require valid common or per-component tabular fractions, but not both."""
    cfg = config
    _reject_conflicting_input_forms(config, "tabular_psd")
    if cfg.tabular_psd is None and cfg.tabular_psd_by_comp is None:
        raise ConfigurationError(
            f"psd_method='tabular' requires tabular_psd with {n_int} fractions, "
            "or tabular_psd_by_comp with fractions for each sized solid component."
        )
    if cfg.tabular_psd_by_comp is None:
        _check_tabular_fractions(n_int, cfg.tabular_psd, "tabular_psd")
        return
    mapping = cfg.tabular_psd_by_comp
    _validate_comp_mapping(config, mapping, "tabular_psd_by_comp")
    for species in cfg.property_package.sized_solid_list:
        _check_tabular_fractions(
            n_int, mapping[species], f"tabular_psd_by_comp[{species!r}]"
        )


def _check_tabular_fractions(n_int, raw, label):
    """Validate that ``raw`` holds ``n_int`` finite nonnegative fractions
    whose sum is within ``_TABULAR_SUM_TOL`` (1e-5) of one."""
    if isinstance(raw, Mapping):
        raise ConfigurationError(
            f"{label} must be a sequence of N = {n_int} fraction values, not a "
            "mapping."
        )
    message = f"{label} must be a sequence containing {n_int} fraction values."
    if not (hasattr(raw, "__len__") and hasattr(raw, "__getitem__")):
        raise ConfigurationError(message)
    try:
        seq = list(raw)
    except TypeError as exc:
        raise ConfigurationError(message) from exc
    if len(seq) != n_int:
        raise ConfigurationError(
            f"{label} must have length N = {n_int} (got {len(seq)})."
        )
    q = [_cfg_float(v, f"{label}[{i}]", lo=0.0) for i, v in enumerate(seq)]
    total = sum(q)
    if abs(total - 1.0) > _TABULAR_SUM_TOL:
        raise ConfigurationError(
            f"{label} sum {total:.8g} deviates from 1 by more than "
            f"{_TABULAR_SUM_TOL:.0e}."
        )


def _validate_selection_breakage_config(config, n_events):
    cfg = config
    n_int = len(config.property_package.size_interval_set)
    _reject_conflicting_input_forms(config, "selection_params")
    _reject_conflicting_input_forms(config, "breakage_params")
    if cfg.recycle == RecycleMode.classification_recycle:
        # Whiten selection supplies the shared recycle classification fractions.
        if cfg.selection_function != SelectionFunction.whiten:
            raise ConfigurationError(
                "recycle='classification_recycle' implements the Whiten method "
                "and requires selection_function='whiten' (got "
                f"'{cfg.selection_function}'). If this selection came from "
                "crusher_equipment, set recycle='single_pass' to use that "
                "equipment's non-Whiten default."
            )
        _reject_unused_nse(n_events, "recycle='classification_recycle'")
        if cfg.selection_params_by_comp is not None:
            raise ConfigurationError(
                "selection_params_by_comp is rejected on "
                "classification_recycle: the classification vector is common; "
                "supply the common selection_params."
            )
        _validate_selection_parameters(
            config, n_int, cfg.selection_params, "selection_params"
        )
        _validate_breakage_inputs(config)
        return
    if cfg.breakage_function == BreakageFunction.t10_appearance:
        raise ConfigurationError(
            "breakage_function='t10_appearance' is valid only in the Whiten "
            "(recycle='classification_recycle') method."
        )
    _validate_selection_inputs(config, n_int)
    _validate_breakage_inputs(config)


def _validate_selection_inputs(config, n_int):
    """Validate common selection data once, or every component value plus the
    equipment key that must stay equal across components."""
    cfg = config
    if cfg.selection_params_by_comp is None:
        _validate_selection_parameters(
            config, n_int, cfg.selection_params, "selection_params"
        )
        return
    mapping = cfg.selection_params_by_comp
    _validate_comp_mapping(config, mapping, "selection_params_by_comp")
    for species in cfg.property_package.sized_solid_list:
        _validate_selection_parameters(
            config,
            n_int,
            mapping[species],
            f"selection_params_by_comp[{species!r}]",
        )
    key = _SELECTION_EQUIPMENT_KEY.get(cfg.selection_function)
    if key is None:
        return
    values = {
        species: _cfg_float(
            mapping[species][key],
            f"selection_params_by_comp[{species!r}] parameter '{key}'",
        )
        for species in cfg.property_package.sized_solid_list
    }
    if len(set(values.values())) > 1:
        raise ConfigurationError(
            f"selection_params_by_comp: equipment parameter '{key}' must "
            f"have the same value for every sized solid component (got {values})."
        )


def _validate_selection_parameters(config, n_int, params, label):
    cfg = config
    sel = cfg.selection_function
    if sel == SelectionFunction.user:
        _check_user_selection_vector(n_int, params, label)
        return
    _check_parameter_keys(
        f"selection_function='{sel}'",
        params,
        label,
        _SELECTION_REQUIRED[sel],
        optional=("K3",) if sel == SelectionFunction.whiten else (),
    )
    _validate_selection_parameter_values(sel, params, label)


def _validate_selection_parameter_values(sel, params, label):
    """Validate numeric domains and parameter relationships for selection functions.

    Calibration ranges are not enforced, allowing meshes outside the
    reference data.
    """

    def _num(key, **bounds):
        return _cfg_float(
            params[key],
            f"selection_function='{sel}' parameter '{key}' in {label}",
            **bounds,
        )

    if sel == SelectionFunction.whiten:
        k1, k2 = _num("K1", positive=True), _num("K2", positive=True)
        if k2 < 1.01 * k1:
            raise ConfigurationError(f"{label}: whiten requires K2 >= 1.01*K1.")
        if "K3" in params:
            _num("K3", positive=True)
    elif sel == SelectionFunction.king:
        _num("CSS", positive=True)
        if not _num("alpha1", positive=True) < _num("alpha2"):
            raise ConfigurationError(f"{label}: king requires 0 < alpha1 < alpha2.")
        _num("n", positive=True)
    elif sel == SelectionFunction.austin:
        _num("S1", lo=0.0, hi=1.0)
        _num("d1", positive=True)
        _num("a", lo=0.0, hi=1.0)
    elif sel == SelectionFunction.vogel_peukert:
        _num("f_mat", positive=True)
        _num("xw_min", lo=0.0)
        _num("v", positive=True)


def _check_user_selection_vector(n_int, params, label):
    """Require a list or tuple with ``n_int`` selection values.

    ``prepare_crushing_data`` checks that each value is finite and in [0, 1].
    """
    if not isinstance(params, (list, tuple)):
        raise ConfigurationError(
            f"selection_function='user' requires '{label}' as a list or "
            "tuple of selection values."
        )
    if len(params) != n_int:
        raise ConfigurationError(
            f"selection_function='user' requires '{label}' to contain "
            f"{n_int} values (got {len(params)})."
        )


def _validate_breakage_inputs(config):
    """Validate common breakage data or each component's individual inputs."""
    cfg = config
    if cfg.breakage_params_by_comp is not None:
        _validate_comp_mapping(
            config, cfg.breakage_params_by_comp, "breakage_params_by_comp"
        )
    for species, params, label in _per_comp_values(
        cfg,
        "breakage_params",
        "breakage_params_by_comp",
        cfg.property_package.sized_solid_list,
    ):
        _validate_breakage_parameters(config, params, label)


def _validate_breakage_parameters(config, params, label):
    """Check parametric breakage keys and scalar values.

    ``prepare_crushing_data`` validates user matrices, appearance tables and
    Ecs tables.
    """
    cfg = config
    brk = cfg.breakage_function
    if brk == BreakageFunction.user:
        return
    pendulum_keys = ["ecs_table"] if cfg.power_law == CrusherPowerLaw.pendulum else []
    _check_parameter_keys(
        f"breakage_function='{brk}'",
        params,
        label,
        _BREAKAGE_REQUIRED[brk] + pendulum_keys,
    )
    _validate_breakage_parameter_values(brk, params, label)


def _check_parameter_keys(choice, params, label, required, optional=()):
    """Require ``params`` to be a mapping with every required key and no
    unknown key.

    ``choice`` names the selected function in messages. ``optional`` keys
    may be omitted.
    """
    if params is not None and not isinstance(params, Mapping):
        raise ConfigurationError(
            f"{choice} requires {label} as a mapping with keys {required}, "
            f"got {type(params).__name__}."
        )
    missing = [key for key in required if params is None or key not in params]
    if missing:
        raise ConfigurationError(
            f"{choice} requires {label} with keys {required}; missing {missing}."
        )
    allowed = set(required) | set(optional)
    unknown = set(params) - allowed
    if unknown:
        raise ConfigurationError(
            f"{choice} got unknown {label} keys {sorted(unknown, key=str)}; "
            f"allowed keys: {sorted(allowed)}."
        )


def _validate_breakage_parameter_values(brk, params, label):
    """Check scalar breakage parameters for finite values and valid ranges."""

    def _num(key, **bounds):
        _cfg_float(
            params[key], f"{brk} breakage parameter '{key}' in {label}", **bounds
        )

    if brk == BreakageFunction.vogel:
        _num("dprime", positive=True)
        _num("q", positive=True)
    elif brk in (BreakageFunction.luckie_austin, BreakageFunction.reid_stewart):
        _num("phi", lo=0.0, hi=1.0)
        _num("gamma", positive=True)
        _num("beta", positive=True)


def _validate_user_breakage_matrix(n_int, params, label):
    """Validate the user breakage matrix and return its rows as float tuples.

    Require finite nonnegative entries, zeros for daughter ``i > j``, and
    column sums within ``_BREAKAGE_COLUMN_SUM_TOL`` of one.
    """
    if (
        not isinstance(params, (list, tuple))
        or len(params) != n_int
        or any(
            not isinstance(row, (list, tuple)) or len(row) != n_int for row in params
        )
    ):
        raise ConfigurationError(
            f"{label} must be a {n_int}x{n_int} matrix B[i][j] using lists "
            "or tuples."
        )
    rows = tuple(
        tuple(_cfg_float(v, f"{label}[{i}][{j}]", lo=0.0) for j, v in enumerate(row))
        for i, row in enumerate(params)
    )
    for i, row in enumerate(rows):
        for j, entry in enumerate(row):
            if i > j and entry != 0.0:
                raise ConfigurationError(
                    f"{label} must be zero for daughter i > parent j "
                    f"(no coarsening; entry [{i}][{j}] = {entry!r})."
                )
    for j in range(n_int):
        defect = math.fsum([row[j] for row in rows] + [-1.0])
        if abs(defect) > _BREAKAGE_COLUMN_SUM_TOL:
            raise ConfigurationError(
                f"{label}: the sum of column {j} minus 1 is {defect!r}; its "
                f"absolute value must not exceed {_BREAKAGE_COLUMN_SUM_TOL}."
            )
    return rows


def _reject_conflicting_input_forms(config, common):
    """Reject a common input and its per-component form being set together."""
    individual = common + "_by_comp"
    cfg = config
    if getattr(cfg, common) is not None and getattr(cfg, individual) is not None:
        raise ConfigurationError(
            f"both {common} and {individual} are set; supply exactly one "
            "member of the pair."
        )


def _validate_comp_mapping(config, mapping, key):
    """Require a dict that keys every sized solid component exactly once."""
    pp = config.property_package
    if not isinstance(mapping, dict):
        raise ConfigurationError(
            f"{key} must be a dict indexed by sized solid components."
        )
    sized = list(pp.sized_solid_list)
    missing = [s for s in sized if s not in mapping]
    unknown = [k for k in mapping if k not in sized]
    if missing or unknown:
        raise ConfigurationError(
            f"{key} must map every sized solid component exactly once: missing "
            f"{missing}, unknown {unknown}."
        )


def prepare_crushing_data(config, product_path, unit_name):
    """Prepare and validate data for the selected product path.

    ``build`` has already checked the configuration. Check selection rows,
    breakage matrices and Ecs tables, normalize tabular fractions and
    breakage columns, and validate recycle denominators where needed. For
    the pendulum law, compute per-size Ecs values from its table.
    """
    cfg = config
    path = product_path
    pp = cfg.property_package
    species = list(pp.sized_solid_list)
    n_int = len(pp.size_interval_set)
    tabular_fractions = selection = breakage = classification = ecs = None
    # Reuse prepared data when multiple components share an input object.
    if path == ProductPath.tabular:
        tabular_fractions = {}
        prepared_rows = {}
        for mineral, raw, label in _per_comp_values(
            cfg, "tabular_psd", "tabular_psd_by_comp", species
        ):
            if id(raw) not in prepared_rows:
                prepared_rows[id(raw)] = _normalize_tabular_fractions(raw)
            tabular_fractions[mineral] = prepared_rows[id(raw)]
    elif cfg.psd_method == PSDMethod.selection_breakage:
        edges, size_char = pp.size_edges_m, pp.size_char_m
        selection, breakage = {}, {}
        prepared_selection, prepared_breakage, prepared_ecs = {}, {}, {}
        if cfg.power_law == CrusherPowerLaw.pendulum:
            ecs = {}
        for mineral, params, label in _per_comp_values(
            cfg, "selection_params", "selection_params_by_comp", species
        ):
            if id(params) not in prepared_selection:
                if cfg.selection_function == SelectionFunction.user:
                    row = tuple(
                        _cfg_float(v, f"{label}[{i}]", lo=0.0, hi=1.0)
                        for i, v in enumerate(params)
                    )
                else:
                    row = tuple(
                        selection_value(cfg.selection_function, d, params)
                        for d in size_char
                    )
                prepared_selection[id(params)] = row
            selection[mineral] = prepared_selection[id(params)]
        for mineral, params, label in _per_comp_values(
            cfg, "breakage_params", "breakage_params_by_comp", species
        ):
            if id(params) not in prepared_breakage:
                if cfg.breakage_function == BreakageFunction.user:
                    rows = _validate_user_breakage_matrix(n_int, params, label)
                else:
                    context = f"breakage_function='{cfg.breakage_function}' ({label}, {mineral!r})"
                    try:
                        rows = breakage_matrix(
                            cfg.breakage_function, edges, size_char, params
                        )
                    except (ArithmeticError, ValueError) as exc:
                        detail = str(exc).rstrip(".") + "."
                        raise ConfigurationError(f"{context}: {detail}") from exc
                prepared_breakage[id(params)] = _normalize_breakage_columns(rows)
                if ecs is not None:
                    try:
                        prepared_ecs[id(params)] = ecs_from_table(
                            params["ecs_table"], float(params["t10"]), size_char
                        )
                    except ValueError as exc:
                        raise ConfigurationError(f"{label}: {exc}") from exc
            breakage[mineral] = prepared_breakage[id(params)]
            if ecs is not None:
                ecs[mineral] = prepared_ecs[id(params)]
        if path == ProductPath.classification_recycle:
            # The shared Whiten row can be taken from any component.
            classification = selection[species[0]]
            minerals_by_breakage = {}
            for mineral in species:
                minerals_by_breakage.setdefault(breakage[mineral], []).append(mineral)
            for rows, minerals in minerals_by_breakage.items():
                validate_recycle_denominators(
                    rows, classification, minerals, unit_name=unit_name
                )
    return PreparedCrusherData(
        tabular_fractions, selection, breakage, classification, ecs
    )


def _normalize_tabular_fractions(raw):
    """Return validated tabular fractions divided by their sum."""
    fractions = tuple(float(v) for v in raw)
    total = math.fsum(fractions)
    return tuple(v / total for v in fractions)


def _normalize_breakage_columns(rows):
    """Return normalized breakage columns with negative roundoff set to zero."""
    rows = [[0.0 if entry < 0.0 else entry for entry in row] for row in rows]
    sums = [math.fsum(column) for column in zip(*rows)]
    return tuple(tuple(v / s for v, s in zip(row, sums)) for row in rows)


def _per_comp_values(cfg, common_key, species_key, species_list):
    """Yield (component, value, label) resolving the common/individual pair."""
    individual = getattr(cfg, species_key)
    if individual is not None:
        for s in species_list:
            yield s, individual[s], f"{species_key}[{s!r}]"
    else:
        common = getattr(cfg, common_key)
        for s in species_list:
            yield s, common, common_key
