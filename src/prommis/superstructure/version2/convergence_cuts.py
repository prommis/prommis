#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Global-Convergence Cuts for the Element-Flow Superstructure
===========================================================

Valid inequalities that tighten the non-convex MIQCP built by
:func:`new_superstructure_function.superstructure.build_model` so that Gurobi's spatial
branch-and-bound (``NonConvex=2``) proves a small optimality gap quickly. The only nonconvexity
is the composition-preserving split ``element_entering_stage = routing_fraction *
pool_element_flow`` (Equation 5); each cut below tightens its relaxation without removing any
optimal solution.

Three cuts are provided -- the ones empirically found to carry convergence:

* :func:`apply_component_invariance_cuts` -- linearizes the split of every *invariant component*
  (any element whose mass fraction is identical across all feeds, e.g. Fe at 0.70). On the
  pre-separation ("invariant") pools each such split becomes an exact linear relation, so the
  bilinear term is deactivated outright. (A two-sided generalization for elements whose fraction
  varies within a narrow band is :func:`apply_bounded_component_cuts`.)
* :func:`apply_routing_sum_cuts` -- ``sum_j routing_fraction[j,i,t] <= 1`` per source pool. A
  pool cannot route out more than it holds. Geometry-free; valid for any feed.
* :func:`apply_ratio_cuts` -- pairwise composition-ratio bounds ``a <= ratio_ab * b`` on the
  invariant pools, from the mediant inequality over the feed compositions. Needs no constant
  component.

All three are **feedstock-agnostic**: nothing is hard-coded to particular elements or feeds.
They auto-detect the relevant geometry from ``feed_composition`` and the built topology, so a
superstructure with a different element set or different feeds is handled unchanged -- a cut for
which the enabling geometry is absent simply adds nothing (it never produces an invalid cut).

These are applied by :func:`build_model` by default (toggled by its ``add_*_cuts`` arguments);
they may also be applied by hand to an already-built model.

Author: Chris Laliwala
"""

import pyomo.environ as pyo


# ---------------------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------------------
def _tracked_elements(m):
    """Return the model's tracked elements as a list."""
    return list(m.fs.tracked_elements_set)


def _detect_invariant_components(feed_composition, elements):
    """Return every element whose mass fraction is identical (and nonzero) across all feeds.

    Each such *invariant component* (e.g. Fe at 0.70 in NdFeB feeds) has that same fraction in
    every invariant node, so its split linearizes exactly. The property is not specific to the
    bulk matrix element: any element with a constant cross-feed fraction qualifies. Returns a list
    (empty if no element has a constant cross-feed fraction).
    """
    feeds = list(feed_composition.keys())
    out = []
    for e in elements:
        fracs = [feed_composition[f][e] for f in feeds]
        if max(fracs) - min(fracs) < 1e-9 and fracs[0] > 1e-12:
            out.append(e)
    return out


def _detect_bounded_components(feed_composition, elements, max_band=0.15):
    """Return ``[(element, c_lo, c_hi), ...]`` for every element with a narrow cross-feed band.

    A *banded component* has a cross-feed fraction spread ``c_hi - c_lo`` within ``max_band`` but
    not exactly zero. Exactly-constant elements are excluded -- they are handled exactly by
    :func:`apply_component_invariance_cuts`, so the two cuts compose without redundancy.

    An element absent from some feed (``c_lo == 0``) still qualifies: its fraction on every
    invariant node lies in ``[0, c_hi]``, which gives a valid *upper* bracket
    ``f <= c_hi/(1-c_hi) * sum(others)`` (the lower bracket degenerates to ``f >= 0`` and is
    dropped by :func:`apply_bounded_component_cuts`). Only an element absent from *all* feeds
    (``c_hi == 0``) is skipped -- there is nothing to bound. Returns an empty list if none qualifies.
    """
    feeds = list(feed_composition.keys())
    out = []
    for e in elements:
        fracs = [feed_composition[f][e] for f in feeds]
        lo, hi = min(fracs), max(fracs)
        if hi <= 1e-12:           # absent in ALL feeds -> nothing to bound
            continue
        band = hi - lo
        if 1e-9 < band <= max_band:   # strictly banded (constant ones go to component-invariance)
            out.append((e, lo, hi))
    return out


def _invariant_pools(m, feed_composition):
    """Return the set of pools whose composition is provably a blend of the feed compositions.

    Computed as a fixpoint over the topology: a produced pool is invariant if every stage that
    produces it is **element-uniform** (equal retention efficiency *and* equal product split
    fraction across all elements) and all of that stage's inputs are themselves invariant. Feeds
    are invariant by definition. Element-uniform stages preserve composition, so an invariant
    pool is a nonnegative combination of the feed composition vectors -- which is exactly what
    the component-invariance and ratio cuts rely on. Separation stages (element-dependent recovery
    or split) break invariance, and everything downstream of them is conservatively excluded.

    This is purely topological and component-agnostic: it makes no reference to any particular element.
    """
    fs = m.fs
    elements = _tracked_elements(m)
    feeds_of, prods_of = {}, {}
    for (j, i) in fs.stage_inlet_pairs:
        feeds_of.setdefault(j, set()).add(i)
    for (j, p) in fs.stage_outlet_pairs:
        prods_of.setdefault(j, set()).add(p)

    def uniform_eff(j):
        v = [pyo.value(fs.stage_efficiency[j, a]) for a in elements]
        return max(v) - min(v) < 1e-9

    def uniform_sf(j, p):
        v = [pyo.value(fs.split_fraction[j, p, a]) for a in elements]
        return max(v) - min(v) < 1e-9

    invariant = set(feed_composition.keys())
    all_pools = {i for (j, i) in fs.stage_inlet_pairs} | {p for (j, p) in fs.stage_outlet_pairs}
    changed = True
    while changed:
        changed = False
        for p in all_pools:
            if p in invariant:
                continue
            producers = [j for j in prods_of if p in prods_of[j]]
            if producers and all(
                uniform_eff(j) and uniform_sf(j, p) and all(inp in invariant for inp in feeds_of.get(j, ()))
                for j in producers
            ):
                invariant.add(p)
                changed = True
    return invariant


# ---------------------------------------------------------------------------------------
# Cut 1 -- component-invariance (exact linearization of every constant-fraction split)
# ---------------------------------------------------------------------------------------
def apply_component_invariance_cuts(m, feed_composition):
    """Linearize the split of every invariant component on the invariant pools.

    An *invariant component* is an element whose mass fraction ``c`` is identical across all feeds
    (auto-detected). It then has that same fraction in every invariant pool, so
    composition-preserving splitting makes its flow a fixed multiple of the others,
    ``component = r * sum(others)`` with ``r = c/(1-c)``. The bilinear ``split_con`` for each such
    component is therefore deactivated and replaced by this linear relation (at both split and pool
    level), removing its bilinear term.

    The bulk matrix element (e.g. Fe at 0.70) is simply the most common case and the one that removes
    the largest-magnitude bilinear terms, but the cut is not restricted to it: if several elements
    share a constant cross-feed fraction they are all linearized here.

    If no element has a constant cross-feed fraction, the cut adds nothing; see
    :func:`apply_bounded_component_cuts` for the banded case and :func:`apply_ratio_cuts` for the
    residual varying elements.

    Args:
        m (ConcreteModel): a built superstructure model.
        feed_composition (dict): ``{feed: {element: mass_fraction}}``.

    Returns:
        ConcreteModel: ``m``, with ``fs.component_invariance_cuts`` added.
    """
    fs = m.fs
    elements = _tracked_elements(m)
    components = _detect_invariant_components(feed_composition, elements)
    if not components:
        print("component-invariance-cuts: no element with a constant cross-feed fraction; skipping.", flush=True)
        return m

    invariant = _invariant_pools(m, feed_composition)
    fs.component_invariance_cuts = pyo.ConstraintList()
    f0 = list(feed_composition)[0]
    n = 0
    for comp in components:
        c = feed_composition[f0][comp]
        r = c / (1.0 - c)
        others = [e for e in elements if e != comp]
        for (j, i) in fs.stage_inlet_pairs:  # split level
            if i not in invariant:
                continue
            for t in fs.operating_years:
                fs.split_con[j, i, comp, t].deactivate()
                fs.component_invariance_cuts.add(
                    expr=fs.element_entering_stage[j, i, comp, t]
                    == r * sum(fs.element_entering_stage[j, i, e, t] for e in others)
                )
                n += 1
        for i in invariant:  # pool-level companion
            for t in fs.operating_years:
                if (i, comp, t) in fs.pool_element_flow:
                    fs.component_invariance_cuts.add(
                        expr=fs.pool_element_flow[i, comp, t]
                        == r * sum(fs.pool_element_flow[i, e, t] for e in others)
                    )

    print(
        f"component-invariance-cuts: linearized {components} on {len(invariant)} invariant pools "
        f"({n} splits)",
        flush=True,
    )
    return m


def apply_bounded_component_cuts(m, feed_composition, max_band=0.15):
    """Bracket every banded component's split on the invariant pools (two-sided generalization).

    When an element's mass fraction is not *identical* across feeds but stays within a narrow band
    ``[c_lo, c_hi]`` (e.g. Nd at 0.20 in one feed, 0.30 in another), the exact equality of
    :func:`apply_component_invariance_cuts` no longer holds -- but every invariant (blend-of-feeds)
    pool still has that element's fraction in ``[c_lo, c_hi]`` (a convex combination of the feed
    fractions). That gives the bracket

        r_lo * sum(others) <= component <= r_hi * sum(others),
        r_lo = c_lo/(1-c_lo),  r_hi = c_hi/(1-c_hi),

    at both pool and split level. Unlike the exact cut this does **not** deactivate the bilinear
    split (an inequality cannot eliminate the variable) -- it brackets it, tightening the McCormick
    envelope. Every element with a band within ``max_band`` is bracketed.

    When an element is absent from some feed (``c_lo == 0``) the lower bracket degenerates to
    ``component >= 0`` (already implied by nonnegativity) and is omitted; only the valid **upper**
    bracket ``component <= r_hi * sum(others)`` is added. So such an element gets a one-sided cut.

    It is **complementary** to :func:`apply_component_invariance_cuts`: exactly-constant elements
    are handled there (and excluded here), so the two can both be enabled without redundancy.
    ``max_band`` is the widest cross-feed fraction spread accepted as a banded component.

    Args:
        m (ConcreteModel): a built superstructure model.
        feed_composition (dict): ``{feed: {element: mass_fraction}}``.
        max_band (float): widest cross-feed fraction spread accepted as a banded component.

    Returns:
        ConcreteModel: ``m``, with ``fs.bounded_component_cuts`` added (when a banded component exists).
    """
    fs = m.fs
    elements = _tracked_elements(m)
    components = _detect_bounded_components(feed_composition, elements, max_band)
    if not components:
        print(f"bounded-component-cuts: no banded element within band {max_band}; skipping.", flush=True)
        return m

    invariant = _invariant_pools(m, feed_composition)
    fs.bounded_component_cuts = pyo.ConstraintList()
    n = 0
    for comp, c_lo, c_hi in components:
        others = [e for e in elements if e != comp]
        r_lo = c_lo / (1.0 - c_lo)
        r_hi = c_hi / (1.0 - c_hi)
        add_lower = r_lo > 1e-12   # absent in some feed (c_lo == 0) -> lower bound is the trivial f >= 0
        for (j, i) in fs.stage_inlet_pairs:  # split level
            if i not in invariant:
                continue
            for t in fs.operating_years:
                others_sum = sum(fs.element_entering_stage[j, i, e, t] for e in others)
                if add_lower:
                    fs.bounded_component_cuts.add(expr=fs.element_entering_stage[j, i, comp, t] >= r_lo * others_sum)
                    n += 1
                fs.bounded_component_cuts.add(expr=fs.element_entering_stage[j, i, comp, t] <= r_hi * others_sum)
                n += 1
        for i in invariant:  # pool level
            for t in fs.operating_years:
                if (i, comp, t) in fs.pool_element_flow:
                    others_sum = sum(
                        fs.pool_element_flow[i, e, t] for e in others if (i, e, t) in fs.pool_element_flow
                    )
                    if add_lower:
                        fs.bounded_component_cuts.add(expr=fs.pool_element_flow[i, comp, t] >= r_lo * others_sum)
                        n += 1
                    fs.bounded_component_cuts.add(expr=fs.pool_element_flow[i, comp, t] <= r_hi * others_sum)
                    n += 1

    bands = ", ".join(
        f"{c}[{lo:.4g},{hi:.4g}]" + ("" if lo > 1e-12 else " (upper-only)")
        for c, lo, hi in components
    )
    print(
        f"bounded-component-cuts: {n} cuts on {{{bands}}} ({len(invariant)} invariant pools)",
        flush=True,
    )
    return m


# ---------------------------------------------------------------------------------------
# Cut 2 -- routing-sum conservation (geometry-free)
# ---------------------------------------------------------------------------------------
def apply_routing_sum_cuts(m):
    """Add ``sum_j routing_fraction[j,i,t] <= 1`` for every source pool ``i`` and year ``t``.

    A pool cannot route out more than it holds. This is implied by the mass balance only after
    dividing out the (variable) pool mass, so the McCormick relaxation does not see it; stating
    it explicitly tightens the relaxation in routing-fraction space and blocks "over-routing"
    relaxation points. It is the linear parent of the RLT split-conservation cut.

    Requires no feed geometry, so it is valid for **any** feedstock and never skips.

    Args:
        m (ConcreteModel): a built superstructure model.

    Returns:
        ConcreteModel: ``m``, with ``fs.routing_sum_cuts`` added.
    """
    fs = m.fs
    consumers = {}
    for (j, i) in fs.stage_inlet_pairs:
        consumers.setdefault(i, set()).add(j)
    fs.routing_sum_cuts = pyo.ConstraintList()
    n = 0
    for i, js in consumers.items():
        for t in fs.operating_years:
            fs.routing_sum_cuts.add(expr=sum(fs.routing_fraction[j, i, t] for j in js) <= 1)
            n += 1

    print(f"routing-sum-cuts: {n} sum-routing<=1", flush=True)
    return m


# ---------------------------------------------------------------------------------------
# Cut 3 -- pairwise composition ratios (component-independent)
# ---------------------------------------------------------------------------------------
def apply_ratio_cuts(m, feed_composition):
    """Add pairwise composition-ratio bounds on the invariant pools.

    Every invariant pool is a nonnegative combination of the feed composition vectors, so for any
    ordered element pair ``(a, b)`` the mediant inequality gives ``a <= ratio_ab * b`` with
    ``ratio_ab = max_f comp[f][a]/comp[f][b]`` -- valid for every blend. A pair is skipped when it
    cannot be bounded (some feed has ``comp[b]=0`` while ``comp[a]>0``, making ``a/b`` unbounded).
    Added at both pool and split level.

    Needs **no constant component**: it works whether or not any element has a constant cross-feed
    fraction (e.g. it still fires when Fe is 0.70 in one feed and 0.71 in another). Any *exact*
    invariant components are excluded from the pairs, since :func:`apply_component_invariance_cuts`
    already pins them by equality -- so in the constant-component case this reduces to the tight
    ratios between the remaining (varying) elements.

    Args:
        m (ConcreteModel): a built superstructure model.
        feed_composition (dict): ``{feed: {element: mass_fraction}}``.

    Returns:
        ConcreteModel: ``m``, with ``fs.ratio_cuts`` added.
    """
    fs = m.fs
    elements = _tracked_elements(m)
    feeds = list(feed_composition)
    # Exact invariant components are handled by apply_component_invariance_cuts; drop them from the pairs.
    components = set(_detect_invariant_components(feed_composition, elements))
    pair_elements = [e for e in elements if e not in components]

    bounds = {}
    for a in pair_elements:
        for b in pair_elements:
            if a == b:
                continue
            # a is unbounded by b if some feed has b == 0 < a -> no finite ratio.
            if any(feed_composition[f][b] <= 1e-12 and feed_composition[f][a] > 1e-12 for f in feeds):
                continue
            cand = [feed_composition[f][a] / feed_composition[f][b] for f in feeds if feed_composition[f][b] > 1e-12]
            if not cand:
                continue
            bounds[(a, b)] = max(cand)

    invariant = _invariant_pools(m, feed_composition)
    fs.ratio_cuts = pyo.ConstraintList()
    n = 0
    for (a, b), ratio in bounds.items():
        for (j, i) in fs.stage_inlet_pairs:  # split level
            if i in invariant:
                for t in fs.operating_years:
                    fs.ratio_cuts.add(
                        expr=fs.element_entering_stage[j, i, a, t]
                        <= ratio * fs.element_entering_stage[j, i, b, t]
                    )
                    n += 1
        for i in invariant:  # pool level
            for t in fs.operating_years:
                if (i, a, t) in fs.pool_element_flow and (i, b, t) in fs.pool_element_flow:
                    fs.ratio_cuts.add(
                        expr=fs.pool_element_flow[i, a, t] <= ratio * fs.pool_element_flow[i, b, t]
                    )
                    n += 1

    tag = f"excluding exact components {sorted(components)}" if components else "no constant component; all element pairs"
    print(
        f"ratio-cuts: {n} pairwise a<=r*b cuts ({tag}); {len(bounds)} finite pairs, "
        f"{len(invariant)} invariant pools",
        flush=True,
    )
    return m
