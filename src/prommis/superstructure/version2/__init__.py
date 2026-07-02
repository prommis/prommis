"""
new_superstructure_function.superstructure
==========================================

The **generic, feedstock-agnostic** element-flow (bilinear stream-mixing) superstructure
function. This subpackage is self-contained — it has no dependency on any case-study data — so
it can be reused for any feedstock or lifted into PrOMMiS as ``prommis.superstructure``. The
included global-convergence cuts (:mod:`.convergence_cuts`) are likewise generic.

Public API:
    build_model               -- construct the generic superstructure model
    ObjectiveFunctionChoice   -- objective-function enum
    SuperstructureScaler      -- optional IDAES scaler
    define_custom_units       -- register the custom Pyomo units
    apply_component_invariance_cuts / apply_bounded_component_cuts / apply_routing_sum_cuts /
    apply_ratio_cuts          -- generic global-convergence cuts (applied by build_model by
                                 default; also callable on an already-built model)
"""

from .objective_function_enums import ObjectiveFunctionChoice
from .superstructure_function import build_model, define_custom_units, SuperstructureScaler
from .convergence_cuts import (
    apply_component_invariance_cuts,
    apply_bounded_component_cuts,
    apply_routing_sum_cuts,
    apply_ratio_cuts,
)

__all__ = [
    "build_model",
    "define_custom_units",
    "SuperstructureScaler",
    "ObjectiveFunctionChoice",
    "apply_component_invariance_cuts",
    "apply_bounded_component_cuts",
    "apply_routing_sum_cuts",
    "apply_ratio_cuts",
]
