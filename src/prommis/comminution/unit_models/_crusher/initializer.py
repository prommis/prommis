#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Initialization of the solid PSD crusher."""

__author__ = "Daison Yancy Caballero"

import math
from collections import namedtuple

from pyomo.common.config import ConfigDict, ConfigValue
from pyomo.contrib.solver.common.util import NoSolutionError
from pyomo.environ import Var, check_optimal_termination, units, value

import idaes.logger as idaeslog
from idaes.core.initialization import InitializationStatus
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import ConfigurationError, InitializationError
from idaes.core.util.model_statistics import large_residuals_set

from prommis.comminution.core.psd_math import finest_attainable_size
from prommis.comminution.core.triangular_block import back_substitution_recycle
from prommis.comminution.core.unit_utils import unit_conversion_factor
from prommis.comminution.functions import distributions
from prommis.comminution.functions.power_laws import (
    APPLICABLE_RANGE_M,
    numeric_pendulum_power,
    numeric_specific_energy,
)
from prommis.comminution.unit_models import _size_reducer_common as sr
from prommis.comminution.unit_models._crusher.shared import (
    CrusherPowerLaw,
    ProductPath,
    WORK_INDEX_BOUNDS,
    _DISTRIBUTION_PATHS,
    _MIN_SHAPE_SIZE,
    _SHAPE_SIZE_LIMIT_RULE,
    _shape_size_upper_bound,
    _close_finest,
    _numeric_crushing_coefficients,
    _p80_from_power,
)

_log = idaeslog.getLogger("prommis.comminution.unit_models.crusher")

# Calculated starting values at one time point. product_flows[s, k],
# recycle_load[s, k] and stage_flows[r, s, k] are flows in kg/s; recycle_load
# and stage_flows are None when the path has none. shape_size (m) is
# None off the distribution paths; sized_feed_flow is the sized feed in kg/s;
# feed_p80 is the feed P80 (m); product_p80 is the analytic product P80 (m) on
# the distribution paths, else None; shape_size_is_free is True only when a free
# distribution characteristic size is inverted from fixed power.
CrusherProductResult = namedtuple(
    "CrusherProductResult",
    [
        "product_flows",
        "stage_flows",
        "recycle_load",
        "shape_size",
        "sized_feed_flow",
        "feed_p80",
        "product_p80",
        "shape_size_is_free",
    ],
)


class CrusherSolidPSDInitializer(sr.SizeReducerInitializerBase):
    """Initialize the crusher and check the solution.

    The initializer checks the feed and operating inputs, calculates starting
    values for every time point, then assigns them and solves once. If an input
    check or calculation fails, it assigns no calculated starting values. IDAES
    loads ``initial_guesses`` before these checks, so those guesses remain after
    such a failure.
    """

    CONFIG = sr.SizeReducerInitializerBase.CONFIG()
    CONFIG.declare("solver_options", ConfigDict(implicit=True))

    CONFIG.solver_options.declare(
        "tol",
        ConfigValue(
            default=1e-8,
            domain=float,
            description="Convergence tolerance for the initialization solve",
        ),
    )

    CONFIG.solver_options.declare(
        "bound_push",
        ConfigValue(
            default=1e-8,
            domain=float,
            description="Minimum distance from a bound for the initial point",
        ),
    )

    def initialize(
        self,
        model,
        initial_guesses=None,
        json_file=None,
        output_level=None,
        exclude_unused_vars=False,
    ):
        """Run the crusher's IDAES initialization lifecycle and return its status.

        IDAES ``InitializerBase.initialize`` calls ``precheck`` outside its
        try/finally block, so a failed precheck would leave the inlet fixed and
        the temporary logger level set. This override restores the entry state
        on failure, clears the logger level on every exit, and resets the
        summary before each run.
        """
        self._reset_summary(model)
        try:
            return super().initialize(
                model,
                initial_guesses=initial_guesses,
                json_file=json_file,
                output_level=output_level,
                exclude_unused_vars=exclude_unused_vars,
            )
        except Exception:
            self.restore_model_state(model)
            raise
        finally:
            self._local_logger_level = None

    def precheck(self, model):
        """Validate feed and operating inputs, then check degrees of freedom."""
        try:
            for t in model.flowsheet().time:
                model.properties_in[t].validate_feed()
            model.validate_sized_feed_flow()
            self._check_work_index(model)
            self._check_power_spec(model)
        except Exception:
            self._update_summary(model, "status", InitializationStatus.PrecheckFailed)
            raise
        super().precheck(model)

    def _reset_summary(self, model):
        """Drop the per-model summary so no stale key survives; status ``none``."""
        self.summary.pop(model, None)
        self._update_summary(model, "status", InitializationStatus.none)

    @staticmethod
    def _check_work_index(model):
        """Require a work index that is available, finite and within its bounds."""
        wi = value(model.bond_work_index, exception=False)
        lo, hi = WORK_INDEX_BOUNDS
        if wi is None or not (math.isfinite(wi) and lo <= wi <= hi):
            raise ConfigurationError(
                f"{model.name}: bond_work_index value {wi!r} kWh/t must be finite "
                f"and within [{lo}, {hi}]."
            )

    @staticmethod
    def _check_power_spec(model):
        """Validate fixed power; require it on the power-spec path."""
        require_fixed = model.product_path == ProductPath.distribution_power_spec
        requirement = (
            "distribution_power_spec requires fixed, finite, nonnegative power"
            if require_fixed
            else "fixed power must be finite and nonnegative"
        )
        for t in model.flowsheet().time:
            if not model.power[t].fixed and not require_fixed:
                continue
            p_kw = value(model.power[t], exception=False)
            if not model.power[t].fixed:
                failed = "not fixed"
            elif p_kw is None or not math.isfinite(p_kw):
                failed = "not finite"
            elif p_kw < 0.0:
                failed = "negative"
            else:
                continue
            raise ConfigurationError(
                f"{model.name}: {requirement} at t = {t}; "
                f"power = {p_kw!r} kW is {failed}."
            )

    def initialization_routine(self, model):
        """Prepare and seed the crusher, solve once and check the solution.

        Return the solver results. IDAES ``postcheck`` rejects nonoptimal results.
        """
        try:
            prepared = model._crushing_data_from_params()
            assignments = self._calculate_starting_values(model, prepared)
        except Exception:
            # No starting values have been assigned, so calculation failures are
            # recorded as precheck failures.
            self._update_summary(model, "status", InitializationStatus.PrecheckFailed)
            raise
        try:
            init_log = self.get_logger(model)
            for target, starting_value in assignments:
                self._seed_if_unfixed(target, starting_value)
            for t in model.flowsheet().time:
                self._seed_passthrough(model, t)
            init_log.info_high(
                "Outlet state seeded from the calculated starting values."
            )
            options = self.config.solver_options.value()
            # IPOPT's default constr_viol_tol (1e-4) can accept a solve that
            # the residual check rejects.
            options.setdefault("constr_viol_tol", options["tol"])
            solver = get_solver("ipopt_v2", solver_options=options)
            solve_log = idaeslog.getSolveLogger(
                model.name, self.get_output_level(), tag="unit"
            )
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                results = solver.solve(model, tee=slc.tee)
            init_log.info_high(f"Initialization solve {idaeslog.condition(results)}.")
        except NoSolutionError as exc:
            # Returning None could let IDAES accept seeded values without a solution.
            self._update_summary(model, "solver_status", False)
            self._update_summary(model, "status", InitializationStatus.Failed)
            raise InitializationError(
                f"{model.name}: the initialization solve returned no solution."
            ) from exc
        except Exception:
            self._update_summary(model, "status", InitializationStatus.Error)
            raise
        if not check_optimal_termination(results):
            return results
        try:
            self._check_solution(model, options["constr_viol_tol"])
        except Exception:
            self._update_summary(model, "status", InitializationStatus.Error)
            raise
        return results

    # Initialization helpers

    def _check_solution(self, model, constr_viol_tol):
        """Check variable values and constraint residuals after an optimal solve.

        Reject NaN or infinite variable values first because
        ``large_residuals_set`` may miss NaNs. It compares each residual times
        its constraint's scaling factor (or 1 if absent) with the solve's
        ``constr_viol_tol`` and records ``None`` for constraints it cannot
        evaluate.
        """
        nonfinite = [
            v.name
            for v in model.component_data_objects(Var, descend_into=True)
            if v.value is not None and not math.isfinite(v.value)
        ]
        if nonfinite:
            raise InitializationError(
                f"{model.name}: after the solve these variables are not finite "
                f"({len(nonfinite)}): {_summarize_names(nonfinite)}."
            )
        large = large_residuals_set(
            model, tol=constr_viol_tol, return_residual_values=True
        )
        if large:
            unevaluable = [c.name for c, residual in large.items() if residual is None]
            if unevaluable:
                raise InitializationError(
                    f"{model.name}: after the solve these constraints could not be "
                    f"evaluated ({len(unevaluable)}): {_summarize_names(unevaluable)}."
                )
            worst = max(large, key=large.get)
            raise InitializationError(
                f"{model.name}: after the solve these constraints have scaled "
                f"residuals above constr_viol_tol = {constr_viol_tol:.3g} "
                f"({len(large)}, largest "
                f"{large[worst]!r} in {worst.name}): "
                f"{_summarize_names([c.name for c in large])}."
            )

    def _calculate_starting_values(self, model, crushing_data):
        """Return ``(VarData, value)`` pairs for every time point.

        Pass-through outlet states are not included; ``_seed_passthrough``
        seeds them. When power is calculated from an unfixed work index, use its current
        value as the initial estimate.
        """
        pp = model.config.property_package
        species = list(pp.sized_solid_list)
        intervals = list(pp.size_interval_set)
        sized_inlet_is_fixed = all(
            model.properties_in[t].flow_mass_sized_comp_size[s, k].fixed
            for t in model.flowsheet().time
            for s in species
            for k in intervals
        )
        # kWh/t equals Wh/kg, the basis of the numeric power-law helpers.
        work_index = value(model.bond_work_index)
        # Match power_eqn's energy floor; _build_power gives the reason.
        floored = model.product_path != ProductPath.distribution_power_spec
        pendulum = model.config.power_law == CrusherPowerLaw.pendulum
        if pendulum:
            self._warn_zero_ecs(model, crushing_data, species, intervals)
        assignments = []
        # Track specific energy when inputs determine the product, or pendulum power
        # at every time point. With fixed sized inlet flows, zero responses cannot
        # determine an unfixed work index or pendulum power factor.
        fit_responses = []
        for t in model.flowsheet().time:
            result = self._calculate_product(model, t, crushing_data)
            if result.stage_flows is not None:
                assignments.extend(
                    (model.psd_stage[t, r, s, k], flow)
                    for (r, s, k), flow in result.stage_flows.items()
                )
            if result.recycle_load is not None:
                assignments.extend(
                    (model.recycle_load[t, s, k], load)
                    for (s, k), load in result.recycle_load.items()
                )
            if result.shape_size is not None:
                assignments.append((model._shape_size_var()[t], result.shape_size))
            if hasattr(model, "feed_p80"):
                assignments.append((model.feed_p80[t], result.feed_p80))
            product = result.product_flows
            f80_live = result.feed_p80
            m_sized = result.sized_feed_flow
            product_p80 = result.product_p80
            if product_p80 is None:
                product_p80 = pp.size_at_passing_from_flows(
                    {s: [product[s, k] for k in intervals] for s in species}
                )
            out = model.properties_out[t]
            # The product rows are in kg/s; seed them in the package flow units.
            to_kg_s = unit_conversion_factor(
                out.flow_mass_sized_comp_size[species[0], intervals[0]],
                units.kg / units.s,
            )
            assignments.extend(
                (out.flow_mass_sized_comp_size[s, k], max(0.0, product[s, k]) / to_kg_s)
                for s in species
                for k in intervals
            )
            if pendulum:
                p_p = numeric_pendulum_power(
                    (crushing_data.ecs[s][k], crushing_data.classification[k], load)
                    for (s, k), load in result.recycle_load.items()
                )
                power = value(model.pendulum_power_factor) * p_p + value(
                    model.no_load_power
                )
                assignments.append((model.power[t], power))
                fit_responses.append(p_p)
            else:
                e_spec = numeric_specific_energy(
                    model.config.power_law, work_index, product_p80, f80_live
                )
                if floored:
                    e_spec = max(0.0, e_spec)
                    assignments.append((model.power[t], 3.6 * m_sized * e_spec))
                    self._log_zero_power_admission(model, t, product_p80, f80_live)
                # Do not add a range warning when fixed power sets product P80.
                if model.product_path != ProductPath.distribution_power_spec:
                    self._warn_outside_applicable_range(model, t, product_p80)
                if not result.shape_size_is_free:
                    fit_responses.append(e_spec)
        if pendulum:
            if (
                sized_inlet_is_fixed
                and not model.pendulum_power_factor.fixed
                and not any(fit_responses)
            ):
                raise ConfigurationError(
                    f"{model.name}: pendulum_power_factor is unfixed, but calculated "
                    "pendulum power is zero at every time point. Measured power cannot "
                    "determine the factor under these conditions. Fix "
                    "pendulum_power_factor."
                )
        elif (
            sized_inlet_is_fixed
            and fit_responses
            and not model.bond_work_index.fixed
            and not any(fit_responses)
        ):
            raise ConfigurationError(
                f"{model.name}: bond_work_index is unfixed and every prescribed "
                "operating condition has exactly zero power response, so no "
                "measured power can identify it. Prescribe a product strictly "
                "finer than the feed, or fix bond_work_index."
            )
        return assignments

    @staticmethod
    def _warn_zero_ecs(model, crushing_data, species, intervals):
        """Warn where component/interval pairs sent to breakage have zero Ecs."""
        zero = [
            (s, k)
            for s in species
            for k in intervals
            if crushing_data.classification[k] > 0.0 and crushing_data.ecs[s][k] == 0.0
        ]
        if zero:
            _log.warning(
                "%s: Ecs is zero for component/interval pairs sent to breakage "
                "(C > 0) %s, so they add no pendulum power. Check the ecs_table "
                "values and whether t10 and the interval sizes fall inside the table.",
                model.name,
                zero,
            )

    def _calculate_product(self, model, t, crushing_data):
        """Return calculated product data for time ``t`` without changing
        model variables."""
        pp = model.config.property_package
        species = list(pp.sized_solid_list)
        intervals = list(pp.size_interval_set)
        edges = list(pp.size_edges_m)
        feed = model.properties_in[t]
        feed_rows = {
            s: [
                value(
                    units.convert(
                        feed.flow_mass_sized_comp_size[s, k],
                        to_units=units.kg / units.s,
                    )
                )
                for k in intervals
            ]
            for s in species
        }
        feed_totals = {s: sum(feed_rows[s]) for s in species}
        m_sized = sum(sum(feed_rows[s][k] for s in species) for k in intervals)
        f80_live = pp.size_at_passing_from_flows(feed_rows)
        stages = load = shape_size = p80_power = None
        shape_size_is_free = False
        path = model.product_path
        if path in _DISTRIBUTION_PATHS:
            product, shape_size, p80_power, shape_size_is_free = (
                self._calculate_distribution(
                    model, t, edges, feed_rows, feed_totals, m_sized, f80_live
                )
            )
        elif path == ProductPath.tabular:
            product = self._calculate_tabular(
                model, crushing_data, t, feed_rows, feed_totals, intervals
            )
        elif path == ProductPath.classification_recycle:
            product, load = self._calculate_recycle(
                crushing_data, feed_totals, feed_rows, intervals
            )
        else:
            product, stages = self._calculate_stress_events(
                model, crushing_data, feed_totals, feed_rows, intervals
            )
        return CrusherProductResult(
            product_flows=product,
            stage_flows=stages,
            recycle_load=load,
            shape_size=shape_size,
            sized_feed_flow=m_sized,
            feed_p80=f80_live,
            product_p80=p80_power,
            shape_size_is_free=shape_size_is_free,
        )

    def _calculate_tabular(
        self, model, crushing_data, t, feed_rows, feed_totals, intervals
    ):
        """Calculate tabular product rows and enforce per-component no-coarsening."""
        fractions, product = {}, {}
        for s, m_s in feed_totals.items():
            table = crushing_data.tabular_fractions[s]
            row = _close_finest([m_s * table[k] for k in intervals], m_s)
            fractions[s] = _close_finest(table)
            product.update({(s, k): row[k] for k in intervals})
        model._check_no_coarsening(fractions, feed_rows, len(intervals), at_time=t)
        return product

    def _calculate_recycle(self, crushing_data, feed_totals, feed_rows, intervals):
        """Calculate Whiten closed-circuit recycle loads and product rows.

        Return ``(product, load)``, both indexed by component and interval.
        """
        product, load = {}, {}
        for s, m_s in feed_totals.items():
            b_s = crushing_data.breakage[s]
            c_vec, f_s = crushing_data.classification, feed_rows[s]
            x = back_substitution_recycle(
                intervals, lambda i, j: b_s[i][j] * c_vec[j], lambda i: f_s[i]
            )
            row = _close_finest([(1.0 - c_vec[k]) * x[k] for k in intervals], m_s)
            load.update({(s, k): x[k] for k in intervals})
            product.update({(s, k): row[k] for k in intervals})
        return product, load

    def _calculate_stress_events(
        self, model, crushing_data, feed_totals, feed_rows, intervals
    ):
        """Apply the configured stress events and return per-component product rows.

        Return ``(product, stages)``. ``stages`` maps stage, component and
        interval to the intermediate flow on the repeated-stress path and is
        ``None`` otherwise.
        """
        stages = {} if model.product_path == ProductPath.repeated_stress else None
        product = {}
        for s, m_s in feed_totals.items():
            selection = crushing_data.selection[s]
            breakage = crushing_data.breakage[s]
            current = feed_rows[s]
            for r in range(1, model.n_stress_events):
                current = self._apply_crushing_event(selection, breakage, current)
                if stages is not None:
                    for k in intervals:
                        stages[r, s, k] = current[k]
            row = _close_finest(
                self._apply_crushing_event(selection, breakage, current), m_s
            )
            product.update({(s, k): row[k] for k in intervals})
        return product, stages

    def _calculate_distribution(
        self, model, t, edges, feed_rows, feed_totals, m_sized, f80_live
    ):
        """Calculate a common product distribution on the size mesh.

        The size-spec path derives its characteristic size from
        ``target_product_p80``. The power-spec path uses a fixed characteristic
        size or estimates one from fixed power and the current work-index value.
        """
        cfg = model.config
        shape, exponent = cfg.distribution_shape, float(cfg.shape_exponent)
        intervals = list(cfg.property_package.size_interval_set)
        shape_size_var = model._shape_size_var()
        shape_size_is_free = False
        if model.product_path == ProductPath.distribution_size_spec:
            p80_power = value(units.convert(model.target_product_p80, to_units=units.m))
            shape_size = distributions.shape_size_from_p80(shape, p80_power, exponent)
        elif shape_size_var[t].fixed:
            shape_size = value(units.convert(shape_size_var[t], to_units=units.m))
            p80_power = shape_size * distributions.p80_ratio(shape, exponent)
        else:
            shape_size_is_free = True
            shape_size, p80_power = self._invert_free_shape_size(
                model, t, edges, shape, exponent, m_sized, f80_live
            )
        issue = self._distribution_checks(shape, exponent, p80_power, shape_size, edges)
        if issue is not None:
            raise ConfigurationError(f"{model.name}: at t = {t}, {issue}")
        fractions = distributions.distribution_fractions(
            shape, edges, shape_size, exponent
        )
        product = {}
        for s, m_s in feed_totals.items():
            row = _close_finest([m_s * fractions[k] for k in intervals], m_s)
            product.update({(s, k): row[k] for k in intervals})
        # If target P80 exceeds feed P80, require the aggregate product to be
        # at least as fine as the feed at every mesh edge.
        if (
            model.product_path == ProductPath.distribution_size_spec
            and p80_power > f80_live
        ):
            # Sum each bin accurately when one component dominates many smaller ones.
            aggregate_feed = [
                math.fsum(row[k] for row in feed_rows.values()) for k in intervals
            ]
            model._check_no_coarsening(
                {None: _close_finest(fractions)},
                {None: aggregate_feed},
                len(intervals),
                at_time=t,
            )
        return product, shape_size, p80_power, shape_size_is_free

    def _invert_free_shape_size(
        self, model, t, edges, shape, exponent, m_sized, f80_live
    ):
        """Estimate the characteristic size and P80 from fixed power at ``t``.

        Use the current work-index value. If it is unfixed and the implied
        distribution is invalid, use the midpoint of the valid size range as
        the starting guess. With a fixed work index, the caller reports the
        invalid distribution.
        """
        p_kw = value(units.convert(model.power[t], to_units=units.kW))
        p80 = _p80_from_power(
            model.config.power_law,
            p_kw,
            m_sized,
            # kWh/t equals Wh/kg, the basis of the numeric power-law helpers.
            value(model.bond_work_index),
            f80_live,
        )
        shape_size = distributions.shape_size_from_p80(shape, p80, exponent)
        if not model.bond_work_index.fixed and self._distribution_checks(
            shape, exponent, p80, shape_size, edges
        ):
            factor = distributions.p80_ratio(shape, exponent)
            shape_size = self._in_bounds_shape_size(
                model, edges, shape, exponent, factor
            )
            p80 = shape_size * factor
        return shape_size, p80

    @staticmethod
    def _in_bounds_shape_size(model, edges, shape, exponent, factor):
        """Return the midpoint of the valid characteristic-size range."""
        upper = _shape_size_upper_bound(shape, edges[-1], exponent)
        lower = max(_MIN_SHAPE_SIZE, finest_attainable_size(edges) / factor)
        if lower > upper:
            raise ConfigurationError(
                f"{model.name}: no valid {shape} characteristic size exists "
                f"at exponent {exponent!r} on this mesh: the calculated range "
                f"[{lower!r}, {upper!r}] m is empty. Change the shape exponent or "
                "mesh."
            )
        return 0.5 * (lower + upper)

    def _log_zero_power_admission(self, model, t, p80_power, f80_live):
        """Log when product P80 above feed P80 gives zero calculated power."""
        if model.product_path == ProductPath.tabular:
            route = "the per-component no-coarsening rule"
        elif model.product_path == ProductPath.distribution_size_spec:
            route = "the aggregate no-coarsening rule"
        else:
            return
        if p80_power > f80_live:
            self.get_logger(model).info(
                f"{model.name}: t = {t}: product P80 = {p80_power!r} m exceeds "
                f"feed P80 = {f80_live!r} m; accepted by {route}; calculated crushing "
                "power uses the zero floor."
            )

    @staticmethod
    def _warn_outside_applicable_range(model, t, p80):
        """Warn, never raise, when the product P80 leaves the law's size range."""
        law = model.config.power_law
        lower, upper = APPLICABLE_RANGE_M.get(law, (None, None))
        if (lower is not None and p80 < lower) or (upper is not None and p80 > upper):
            _log.warning(
                "%s: at t = %s, power_law='%s' is applied with product P80 = %.4g m "
                "outside its applicable range.",
                model.name,
                t,
                law,
                p80,
            )

    @staticmethod
    def _apply_crushing_event(selection, breakage, incoming):
        """Apply one selection and breakage event to an incoming size row.

        The coarser bins use the same coefficients as the symbolic crusher
        equations. Set the finest bin so the output sums to the incoming total.
        """
        n_int = len(incoming)
        a = _numeric_crushing_coefficients(selection, breakage, n_int)
        interior = [0.0] + [
            sum(a[k][i] * incoming[i] for i in range(k, n_int)) for k in range(1, n_int)
        ]
        return _close_finest(interior, sum(incoming))

    @staticmethod
    def _distribution_checks(shape, exponent, p80, shape_size, edges):
        """Return the first failed distribution check, or ``None``.

        Check the P80 mesh band and the characteristic-size bounds. Check fixed
        sizes explicitly because Pyomo only warns when a variable is fixed
        outside its bounds.
        """
        floor = finest_attainable_size(edges)
        if not (math.isfinite(p80) and floor <= p80 <= edges[-1]):
            return (
                f"product P80 {p80!r} m must be finite and within the attainable "
                f"mesh band [{floor:.6g}, {edges[-1]:.6g}] m. Review the "
                "product-size or power specification and the mesh."
            )
        upper = _shape_size_upper_bound(shape, edges[-1], exponent)
        if not _MIN_SHAPE_SIZE <= shape_size <= upper:
            rule = f"; {_SHAPE_SIZE_LIMIT_RULE[shape]}" if shape_size > upper else ""
            return (
                f"the distribution characteristic size {shape_size!r} m must lie "
                f"within [{_MIN_SHAPE_SIZE}, {upper:.6g}] m{rule}. Review the "
                "distribution specification and mesh."
            )
        return None


def _summarize_names(names, limit=10):
    """Format names for an error message, showing at most ``limit``."""
    shown = ", ".join(names[:limit])
    return shown if len(names) <= limit else f"{shown} and {len(names) - limit} more"
