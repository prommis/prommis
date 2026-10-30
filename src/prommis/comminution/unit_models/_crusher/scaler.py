#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Scaling of the solid PSD crusher."""

__author__ = "Daison Yancy Caballero"

import math

from pyomo.common.collections import ComponentSet
from pyomo.core.expr.visitor import identify_variables
from pyomo.environ import units, value

from prommis.comminution.core.psd_math import finest_attainable_size
from prommis.comminution.core.unit_utils import (
    FactorSource,
    fixed_value_or_none,
    inverse_magnitude_factor,
    inverse_sum_of_nominals,
    required_scaling_factor,
    scale_constraints_at_time,
    unit_conversion_factor,
)
from prommis.comminution.functions import distributions
from prommis.comminution.functions.power_laws import numeric_specific_energy
from prommis.comminution.unit_models import _size_reducer_common as sr
from prommis.comminution.unit_models._crusher.shared import (
    CrusherPowerLaw,
    ProductPath,
    MIN_SIZED_FEED_KG_S,
    _shape_size_upper_bound,
    _p80_from_power,
)


class CrusherSolidPSDScaler(sr.SizeReducerScalerBase):
    """Scaler for the solid PSD crusher.

    Delegates state-block scaling and scales crusher variables and constraints
    for the selected product path. With factor_source="input_based", the crusher
    scaler uses configured inputs and fixed values, without reading the current
    values of unfixed crusher variables. With factor_source="current_values",
    it uses those current values, so scaling factors depends on model
    initialization.
    """

    # Smallest power nominal (kW) used for a power scaling factor.
    POWER_NOMINAL_FLOOR_KW = 1.0
    # Smallest nominal for an unfixed work index (kWh/t), pendulum power factor
    # (dimensionless) or no-load power (kW).
    PARAMETER_NOMINAL_FLOOR = 1.0
    # Smallest size nominal (m) for the characteristic size, target P80 and
    # feed P80.
    SIZE_NOMINAL_FLOOR_M = 1.0e-6
    # Current-values working-flow floor, as a fraction of max(sized feed, 1 kg/s).
    WORKING_FLOW_RELATIVE_FLOOR = 1.0e-6
    # In current-values mode, outlet and working-flow bin nominals are at least
    # this fraction of each component's summed inlet flow magnitudes.
    CURRENT_VALUES_BIN_NOMINAL_FLOOR = 1.0e-2
    # Assumed feed-to-product P80 ratio for input-based power on the matrix paths.
    INPUT_BASED_REDUCTION_RATIO = 4.0

    def variable_scaling_routine(self, model, overwrite=False, submodel_scalers=None):
        """Delegate state scaling and set crusher variable factors.

        In input-based mode, each outlet sized flow uses its component's summed
        inlet nominal. In current-values mode, its nominal is based on the
        current flow magnitude and floored at ``CURRENT_VALUES_BIN_NOMINAL_FLOOR``
        times the sum of that component's inlet flow magnitudes. Existing outlet
        factors are kept unless ``overwrite`` is true; factors set by an outlet
        scaler in ``submodel_scalers`` are always kept.
        """
        input_based = self.config.factor_source == FactorSource.input_based
        if input_based:
            super().variable_scaling_routine(
                model, overwrite=overwrite, submodel_scalers=submodel_scalers
            )
        else:
            kept = ComponentSet()
            if not overwrite:
                kept = ComponentSet(
                    v
                    for t in model.flowsheet().time
                    for v in model.properties_out[t].flow_mass_sized_comp_size.values()
                    if self.get_scaling_factor(v) is not None
                )
            super().variable_scaling_routine(
                model, overwrite=overwrite, submodel_scalers=submodel_scalers
            )
            if not (submodel_scalers and model.properties_out in submodel_scalers):
                for t in model.flowsheet().time:
                    self._floor_outlet_current_values(model, t, kept)
        estimable = [(model.bond_work_index, model.config.bond_work_index)]
        if model.config.power_law == CrusherPowerLaw.pendulum:
            params = model.config.pendulum_power_params
            estimable += [
                (model.pendulum_power_factor, params["power_factor"]),
                (model.no_load_power, params["no_load_power"]),
            ]
        for var, configured in estimable:
            if not var.fixed:
                guess = float(configured) if input_based else value(var)
                self.set_variable_scaling_factor(
                    var,
                    inverse_magnitude_factor(guess, self.PARAMETER_NOMINAL_FLOOR),
                    overwrite=overwrite,
                )
        edges = list(model.config.property_package.size_edges_m)
        shape_size_var = model._shape_size_var()
        for t in model.flowsheet().time:
            self.set_variable_scaling_factor(
                model.power[t],
                (
                    self._input_based_power_factor(model, t, edges)
                    if input_based
                    else inverse_magnitude_factor(
                        value(model.power[t]), self.POWER_NOMINAL_FLOOR_KW
                    )
                ),
                overwrite=overwrite,
            )
            if shape_size_var is not None:
                self.set_variable_scaling_factor(
                    shape_size_var[t],
                    (
                        self._input_based_shape_size_factor(model, t, edges)
                        if input_based
                        else inverse_magnitude_factor(
                            value(shape_size_var[t]), self.SIZE_NOMINAL_FLOOR_M
                        )
                    ),
                    overwrite=overwrite,
                )
            if hasattr(model, "feed_p80"):
                self.set_variable_scaling_factor(
                    model.feed_p80[t],
                    (
                        self._input_based_feed_p80_factor(model, t, edges)
                        if input_based
                        else inverse_magnitude_factor(
                            value(model.feed_p80[t]), self.SIZE_NOMINAL_FLOOR_M
                        )
                    ),
                    overwrite=overwrite,
                )
            self._scale_working_flows(model, t, input_based, overwrite)

    def constraint_scaling_routine(self, model, overwrite=False, submodel_scalers=None):
        """Scale state, pass-through, product-path, and power constraints.

        Requires variable scaling factors to be set first.
        """
        super().constraint_scaling_routine(
            model, overwrite=overwrite, submodel_scalers=submodel_scalers
        )
        # Stress and recycle paths write product rows k >= 1 in kg/s. Their
        # finest row, and every product row on the other paths, uses the
        # property package's flow basis.
        si_rows = model.product_path in (
            ProductPath.repeated_stress,
            ProductPath.classification_recycle,
        )
        pp = model.config.property_package
        for t in model.flowsheet().time:
            feed = model.properties_in[t]
            # Use each component's feed estimate when scaling its rows.
            inv_comp = {
                s: inverse_sum_of_nominals(
                    self,
                    [
                        feed.flow_mass_sized_comp_size[s, k]
                        for k in pp.size_interval_set
                    ],
                )
                for s in pp.sized_solid_list
            }
            to_kg_s = self._solid_flow_to_kg_s(feed)
            internal_flow_components = ComponentSet(
                getattr(model, n)
                for n in ("psd_stage", "recycle_load")
                if hasattr(model, n)
            )
            for name in (
                "product_psd_eqn",
                "recycle_load_eqn",
                "psd_stage_eqn",
                "no_coarsening_comp_eqn",
            ):
                for idx, con in getattr(model, name, {}).items():
                    if idx[0] != t:
                        continue
                    in_kg_s = name in ("recycle_load_eqn", "psd_stage_eqn") or (
                        name == "product_psd_eqn" and si_rows and idx[-1] != 0
                    )
                    # Limit the factor if row flow exceeds the component feed estimate.
                    factor = min(
                        inv_comp[idx[-2]],
                        self._largest_flow_term_factor(
                            con, internal_flow_components, to_kg_s
                        ),
                    )
                    self.set_constraint_scaling_factor(
                        con,
                        factor / (to_kg_s if in_kg_s else 1.0),
                        overwrite=overwrite,
                    )
            scale_constraints_at_time(
                self, model, [("no_coarsening_eqn", 1.0)], t, overwrite=overwrite
            )
            self.set_constraint_scaling_factor(
                model.power_eqn[t],
                required_scaling_factor(self, model.power[t]),
                overwrite=overwrite,
            )
            if hasattr(model, "p80_shape_size_eqn"):
                self.set_constraint_scaling_factor(
                    model.p80_shape_size_eqn[t],
                    inverse_magnitude_factor(
                        value(model.target_product_p80), self.SIZE_NOMINAL_FLOOR_M
                    ),
                    overwrite=overwrite,
                )
            if hasattr(model, "feed_p80_eqn"):
                self.set_constraint_scaling_factor(
                    model.feed_p80_eqn[t],
                    required_scaling_factor(self, model.feed_p80[t]),
                    overwrite=overwrite,
                )

    # Scaling helpers

    @staticmethod
    def _solid_flow_to_kg_s(feed):
        """Return the factor that converts the feed's solid flow units to kg/s."""
        return unit_conversion_factor(
            next(iter(feed.flow_mass_sized_comp_size.values())), units.kg / units.s
        )

    def _largest_flow_term_factor(self, constraint, internal_flow_components, to_kg_s):
        """Return the inverse of the largest flow nominal in ``constraint``.

        Use unfixed sized flows from the state blocks and internal stage and
        recycle flows. Convert the internal flows from kg/s to the property
        package's flow units. Ignore other variables and use 1e-8 as the
        minimum nominal.
        """
        nominals = []
        for var in identify_variables(constraint.body, include_fixed=False):
            comp = var.parent_component()
            if comp in internal_flow_components:
                nominals.append(1.0 / (required_scaling_factor(self, var) * to_kg_s))
            elif comp.local_name == "flow_mass_sized_comp_size":
                nominals.append(1.0 / required_scaling_factor(self, var))
        return 1.0 / max(max(nominals, default=0.0), 1.0e-8)

    def _input_based_power_factor(self, model, t, edges):
        """Return the input-based factor for ``power[t]``.

        A fixed power uses its own magnitude. Otherwise use
        ``1 / max(P_est, POWER_NOMINAL_FLOOR_KW)`` for the power estimate in
        kW, with ``POWER_NOMINAL_FLOOR_KW`` (1 kW) as the denominator floor.
        Use ``INPUT_BASED_DEFAULT_POWER_FACTOR`` when there is no estimate or it
        is not finite.
        """
        p_kw = fixed_value_or_none(model.power[t])
        if p_kw is not None:
            return inverse_magnitude_factor(p_kw, self.POWER_NOMINAL_FLOOR_KW)
        if model.config.power_law == CrusherPowerLaw.pendulum:
            power_estimate = self._pendulum_power_estimate_kw(model, t)
        else:
            power_estimate = self._energy_law_power_estimate_kw(model, t, edges)
        if power_estimate is None or not math.isfinite(power_estimate):
            return self.INPUT_BASED_DEFAULT_POWER_FACTOR
        return 1.0 / max(power_estimate, self.POWER_NOMINAL_FLOOR_KW)

    def _pendulum_power_estimate_kw(self, model, t):
        """Estimate total crusher power in kW as if all sized feed broke once at
        the largest Ecs value; return ``None`` without fixed parameters and feed."""
        a = fixed_value_or_none(model.pendulum_power_factor)
        p_n = fixed_value_or_none(model.no_load_power)
        feed = self._fixed_sized_feed(model, t)
        if a is None or p_n is None or feed is None:
            return None
        _, sized_feed = feed
        e_max = max(value(e) for e in model.ecs.values())
        return p_n + 3.6 * a * sized_feed * e_max

    def _energy_law_power_estimate_kw(self, model, t, edges):
        """Estimate Bond, Rittinger, or Kick power in kW from an estimated P80.

        Return ``None`` on the power-spec path or if fixed sized feed is
        unavailable or below ``MIN_SIZED_FEED_KG_S`` (1e-6 kg/s).
        """
        if model.product_path == ProductPath.distribution_power_spec:
            return None
        feed = self._fixed_feed_rows_and_f80(model, t)
        if feed is None:
            return None
        rows, sized_feed, f80 = feed
        p80 = self._input_based_product_p80(model, edges, rows, sized_feed, f80)
        e_spec = max(
            0.0,
            numeric_specific_energy(
                model.config.power_law, self._input_based_work_index(model), p80, f80
            ),
        )
        return 3.6 * sized_feed * e_spec

    @classmethod
    def _input_based_product_p80(cls, model, edges, rows, sized_feed, f80):
        """Estimate product P80 in meters for input-based power scaling.

        Tabular products use fixed component totals and prepared tabular fractions;
        the size-spec path uses ``target_product_p80``. Matrix paths divide the feed F80
        by ``INPUT_BASED_REDUCTION_RATIO``, bounded by the finest attainable P80 and
        the feed F80.
        """
        if model.product_path == ProductPath.tabular:
            tables = model._crushing_data_from_params().tabular_fractions
            n_int = len(edges) - 1
            interior = [
                sum(sum(rows[s]) * tables[s][k] for s in rows) for k in range(1, n_int)
            ]
            return model.config.property_package.size_at_passing_from_flows(
                {None: [sized_feed - sum(interior)] + interior}
            )
        if model.product_path == ProductPath.distribution_size_spec:
            return value(units.convert(model.target_product_p80, to_units=units.m))
        return min(
            max(f80 / cls.INPUT_BASED_REDUCTION_RATIO, finest_attainable_size(edges)),
            f80,
        )

    def _input_based_shape_size_factor(self, model, t, edges):
        """Return the input-based factor for the distribution's characteristic size.

        Use the fixed size when available. Otherwise derive it from the
        specified P80, or on the power-spec path from fixed power, feed, and
        work index. Zero fixed power uses feed F80. If the power-spec estimate
        is unavailable or the fixed power cannot be inverted, use half the
        build upper limit. Use ``SIZE_NOMINAL_FLOOR_M`` (1e-6 m) as the
        denominator floor.
        """
        cfg = model.config
        shape = cfg.distribution_shape
        n = float(cfg.shape_exponent)
        shape_size = fixed_value_or_none(model._shape_size_var()[t])
        if shape_size is not None:
            return 1.0 / max(shape_size, self.SIZE_NOMINAL_FLOOR_M)
        if model.product_path == ProductPath.distribution_size_spec:
            p80_hat = value(units.convert(model.target_product_p80, to_units=units.m))
        else:
            p_kw = fixed_value_or_none(model.power[t])
            feed = self._fixed_feed_rows_and_f80(model, t)
            p80_hat = None
            if p_kw is not None and feed is not None:
                _, m_sized, f80_live = feed
                try:
                    p80_hat = _p80_from_power(
                        cfg.power_law,
                        p_kw,
                        m_sized,
                        self._input_based_work_index(model),
                        f80_live,
                    )
                except (ValueError, ArithmeticError):
                    p80_hat = None
            if p80_hat is None or not (math.isfinite(p80_hat) and p80_hat > 0.0):
                shape_size_ub = _shape_size_upper_bound(shape, edges[-1], n)
                return 1.0 / max(0.5 * shape_size_ub, self.SIZE_NOMINAL_FLOOR_M)
        return 1.0 / max(
            distributions.shape_size_from_p80(shape, p80_hat, n),
            self.SIZE_NOMINAL_FLOOR_M,
        )

    def _input_based_feed_p80_factor(self, model, t, edges):
        """Return the input-based scaling factor for ``feed_p80[t]``.

        Use the smoothed P80 from fixed sized-feed flows when available; otherwise
        use half the top mesh edge. Return the reciprocal after flooring the size at
        ``SIZE_NOMINAL_FLOOR_M`` (1e-6 m).
        """
        feed = self._fixed_feed_rows_and_f80(model, t)
        f80 = 0.5 * edges[-1] if feed is None else feed[2]
        return 1.0 / max(f80, self.SIZE_NOMINAL_FLOOR_M)

    @staticmethod
    def _input_based_work_index(model):
        """Return the fixed work index, or the configured value when it is unfixed."""
        wi = fixed_value_or_none(model.bond_work_index)
        if wi is None:
            wi = float(model.config.bond_work_index)
        return wi

    def _fixed_sized_feed(self, model, t):
        """Return fixed sized-feed magnitudes and their total in kg/s.

        Return ``None`` if any sized-feed value is unset or unfixed, or if the
        total is below ``MIN_SIZED_FEED_KG_S`` (1e-6 kg/s).
        """
        feed_si = self._fixed_feed_solid_si(model, t)
        if feed_si is None:
            return None
        m_sized = sum(feed_si.values())
        if m_sized < MIN_SIZED_FEED_KG_S:
            return None
        return feed_si, m_sized

    def _fixed_feed_rows_and_f80(self, model, t):
        """Return ``(rows, total, f80)`` for the fixed feed at time ``t``, or
        ``None`` when ``_fixed_sized_feed`` gives none.

        ``rows`` maps each sized solid component to its interval flows in kg/s,
        ``total`` is the sized feed in kg/s and ``f80`` the feed P80 in m.
        """
        feed = self._fixed_sized_feed(model, t)
        if feed is None:
            return None
        feed_si, sized_feed = feed
        pp = model.config.property_package
        species = list(pp.sized_solid_list)
        order = list(pp.size_interval_set)
        rows = {s: [feed_si[s, k] for k in order] for s in species}
        return rows, sized_feed, pp.size_at_passing_from_flows(rows)

    def _component_value_totals(self, model, t):
        """Return each sized component's sum of current inlet flow magnitudes."""
        totals = {}
        for (s, _), var in model.properties_in[t].flow_mass_sized_comp_size.items():
            totals[s] = totals.get(s, 0.0) + abs(value(var))
        return totals

    def _floor_outlet_current_values(self, model, t, kept):
        """Floor current-value outlet sized-flow nominals at time ``t``.

        For each Var not in ``kept``, the nominal becomes
        ``max(nominal, floor[s])``, where ``floor[s]`` is
        ``CURRENT_VALUES_BIN_NOMINAL_FLOOR`` times the sum of that
        component's current inlet flow magnitudes.
        """
        totals = self._component_value_totals(model, t)
        for (s, _), var in model.properties_out[t].flow_mass_sized_comp_size.items():
            floor = self.CURRENT_VALUES_BIN_NOMINAL_FLOOR * totals[s]
            if var not in kept and floor > 0.0:
                self.set_variable_scaling_factor(
                    var,
                    min(required_scaling_factor(self, var), 1.0 / floor),
                    overwrite=True,
                )

    def _scale_working_flows(self, model, t, input_based, overwrite):
        """Scale recycle and stress-stage flows when present at time ``t``.

        Input-based mode uses each component's summed inlet nominal, as for
        outlet bins. Current-values mode uses each working flow's magnitude,
        with a minimum set by the larger of a feed-based floor and a floor
        based on that component's current inlet flows.
        """
        feed = model.properties_in[t]
        if input_based:
            flow_si = self._solid_flow_to_kg_s(feed)
            nom_si = {
                idx: flow_si
                / required_scaling_factor(self, feed.flow_mass_sized_comp_size[idx])
                for idx in feed.flow_mass_sized_comp_size
            }
            totals = self._component_nominal_totals(nom_si)
        else:
            m_sized_si = value(
                units.convert(feed.flow_mass_sized, to_units=units.kg / units.s)
            )
            floor_s = self.WORKING_FLOW_RELATIVE_FLOOR * max(m_sized_si, 1.0)
            flow_si = self._solid_flow_to_kg_s(feed)
            floors = {
                s: max(floor_s, self.CURRENT_VALUES_BIN_NOMINAL_FLOOR * tot * flow_si)
                for s, tot in self._component_value_totals(model, t).items()
            }
        for name in ("recycle_load", "psd_stage"):
            var = getattr(model, name, None)
            if var is None:
                continue
            for idx in var:
                if idx[0] != t:
                    continue
                if input_based:
                    factor = 1.0 / totals[idx[-2]]
                else:
                    factor = inverse_magnitude_factor(value(var[idx]), floors[idx[-2]])
                self.set_variable_scaling_factor(var[idx], factor, overwrite=overwrite)
