#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Crusher
=======

Author: Lingyan Deng

The Crusher module includes power consumption for solid crushing. It is a function
of particle size distribution, mass flow rate, and bond work index.

This model contains stage-specific applicability guidance so the unit does not accept
an arbitrary feed/product particle-size combination without context.

The stage categories are interpreted using product-P80 intervals:

- primary: product P80 >= 10 cm
- secondary: 2 cm <= product P80 <= 10 cm - eps
- tertiary: 0.5 cm <= product P80 <= 2 cm - eps

Here eps is a small positive tolerance used so the algebraic constraints
remain non-overlapping in Pyomo.
    
If the requested product P80 is below 0.5 cm, the target is treated as finer than
this modeled crushing-stage range.

The stage categories is referenced from: Upadhyay, R. K. "Mining, Mineral Beneficiation,
and Environment." Geology and Mineral Resources. Singapore: Springer Nature Singapore, 2025. 799-858.

The model can be configured either by crusher stage or by equipment name. When an
unambiguous equipment name is supplied, the corresponding stage is inferred.

A utility method, ``check_applicability()``, is provided to show warning messages
when the current product size falls outside the configured stage range, or when the
requested target is finer than the modeled crushing-stage range.

Degrees of Freedom
------------------

A Crusher module has two degrees of freedom, which are the output of
"particle_size_median" and "particle_size_width".

Model Structure
---------------

The Crusher model includes one inlet Port (inlet) and one outlet Port (outlet).
The properties of the Crusher Unit model is mainly the particle size distribution.

Additional Constraints
----------------------

Crusher adds additional constraints to calculate the work required to crush the
particles, enforce crushing direction (feed P80 >= product P80), and optionally
enforce stage-specific product-P80 applicability limits.

.. math:: work_{t} = 10 * m_{t, in} * BWI * \left(\frac{1}{\sqrt{P_{t, prod, 80}}} - \frac{1}{\sqrt{P_{t, feed, 80}}}\right)

where :math:`work_{t}` is the work required to crush the particles at :math:`t`
time, 10 is an empirical value and should not be changed, :math:`m_{t, in}` is
the inlet mass flow rate at :math:`t` time, :math:`BWI` is the bond work index of
particles, :math:`P_{t, prod, 80}` is production particle size with 80% passing
the mesh at :math:`t` time, :math:`P_{t, feed, 80}` is feed particle size with
80% passing the mesh at :math:`t` time.

Expressions
-----------

Crusher includes two expressions to calculate the size of particles that has 80% 
passing the mesh for both feed and product particles.

.. math:: P_{t, feed, 80} = \frac{S_{t, in, median}}{unit} * \left(-\log(1 - 0.8)\right)^{\frac{SW_{t, in}}{2}}

where :math:`P_{t, feed, 80}` is the feed particle size that has 80% passing the mesh
at :math:`t` time, :math:`\frac{S_{t, in, median}}{unit}` is the median particle size
of input at :math:`t` time and unitless. The default particle size is in micrometer.
The :math:`SW_{t, in}` is the particle size width of input at :math:`t` time.

.. math:: P_{t, prod, 80} = \frac{S_{t, out, median}}{unit} * \left(-\log(1 - 0.8)\right)^{\frac{SW_{t, out}}{2}}

where :math:`P_{t, prod, 80}` is the product particle size that has 80% passing the mesh
at :math:`t` time, :math:`\frac{S_{t, out, median}}{unit}` is the median particle size of
output at :math:`t` time and unitless. The default particle size is in micrometer.
The :math:`SW_{t, in}` is the particle size width of output at :math:`t` time.

Variables
---------

Crusher add the following additional variables beyond those created in property packages.

================ ====== =============================================================================
Variable         Name   Notes
================ ====== =============================================================================
:math:`work_{t}`  work
================ ====== =============================================================================

"""

from functools import partial
import math

from pyomo.common.config import ConfigDict, ConfigValue, In
from pyomo.environ import Constraint, Var, log, value
from pyomo.environ import units as pyunits

import idaes.logger as idaeslog
from idaes.core import UnitModelBlockData, declare_process_block_class, useDefault
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.tables import create_stream_table_dataframe

_log = idaeslog.getLogger(__name__)

STAGE_P80_EPS = 1e-6 * pyunits.cm
CRUSHER_STAGE_ORDER = ("primary", "secondary", "tertiary")

CRUSHER_STAGE_DATA = {
    "primary": {
        "description": "Primary crushing",
        "equipment_names": ("jaw", "gyratory"),
        "applicable_product_p80_lower": 10 * pyunits.cm,
        "applicable_product_p80_upper": None,
        "default_product_p80": 12 * pyunits.cm,
    },
    "secondary": {
        "description": "Secondary crushing",
        "equipment_names": ("cone", "roll1"),
        "applicable_product_p80_lower": 2 * pyunits.cm,
        "applicable_product_p80_upper": 10 * pyunits.cm,
        "default_product_p80": 5 * pyunits.cm,
    },
    "tertiary": {
        "description": "Tertiary crushing",
        "equipment_names": ("roll2", "short_head_cone", "hammer_mill"),
        "applicable_product_p80_lower": 0.5 * pyunits.cm,
        "applicable_product_p80_upper": 2 * pyunits.cm,
        "default_product_p80": 1 * pyunits.cm,
    },
}

CRUSHER_EQUIPMENT_TO_STAGE = {
    "jaw": "primary",
    "gyratory": "primary",
    "cone": "secondary",
    "roll1": "secondary",
    "roll2": "tertiary",
    "short_head_cone": "tertiary",
    "hammer_mill": "tertiary",
}

def _size_to_native_units(size_value, native_units):
    """
    Convert a size quantity to the numeric magnitude associated with native_units.

    Example: if size_value is 1 cm and native_units is um, return 10000.
    """
    return value(pyunits.convert(size_value, to_units=native_units) / native_units)

def _stage_eps_to_native_units(native_units):
    return _size_to_native_units(STAGE_P80_EPS, native_units)

def _equipment_label(equipment_names):
    return ", ".join(name.replace("_", " ") for name in equipment_names)

def _range_bound_to_native_units(size_value, native_units):
    if size_value is None:
        return None
    return _size_to_native_units(size_value, native_units)

def _get_stage_bounds_native(stage_name, native_units):
    stage_data = CRUSHER_STAGE_DATA[stage_name]
    lower = _range_bound_to_native_units(
        stage_data["applicable_product_p80_lower"], native_units
    )
    upper = _range_bound_to_native_units(
        stage_data["applicable_product_p80_upper"], native_units
    )
    return lower, upper

def _product_p80_in_stage_range(stage_name, product_p80_native, native_units):
    """
    Check whether product P80 lies in the stated stage-classification interval.

    The interval rules follow the user-specified classification:
    - primary: product P80 >= 10 cm
    - secondary: 2 cm <= product P80 <= 10 cm - eps
    - tertiary: 0.5 cm <= product P80 <= 2 cm - eps
    
    Here eps is a small positive tolerance used so the algebraic constraints
    remain non-overlapping in Pyomo.
    """
    lower, upper = _get_stage_bounds_native(stage_name, native_units)
    eps = _stage_eps_to_native_units(native_units)

    if stage_name == "primary":
        return product_p80_native >= lower

    if stage_name in ("secondary", "tertiary"):
        return lower <= product_p80_native <= upper - eps

    raise ValueError(f"Unsupported crusher stage '{stage_name}'.")

def _recommend_stage_from_product_p80_native(product_p80_native, native_units):
    """
    Recommend a crusher stage for a product P80 using explicit stage intervals.

    Returns
    -------
    tuple[str, bool]
        (recommended_stage, in_modeled_range), where in_modeled_range is False if
        the requested product is finer than 0.5 cm.
    """
    if _product_p80_in_stage_range("primary", product_p80_native, native_units):
        return "primary", True

    if _product_p80_in_stage_range("secondary", product_p80_native, native_units):
        return "secondary", True

    if _product_p80_in_stage_range("tertiary", product_p80_native, native_units):
        return "tertiary", True

    return "tertiary", False

@declare_process_block_class("Crusher")
class CrusherData(UnitModelBlockData):
    CONFIG = ConfigDict()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Crush unit is steady-state only""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Crush unit has no holdup.""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigDict(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigDict with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "crusher_stage",
        ConfigValue(
            default=None,
            domain=In([None] + list(CRUSHER_STAGE_DATA.keys())),
            description="Crusher stage",
            doc="""Select the applicable crushing stage.

Valid values:
- None: infer from equipment, or default to secondary if no equipment is given
- primary: jaw, gyratory crusher
- secondary: cone, roll1 crusher
- tertiary: roll2, short-head cone crusher, hammer mill
""",
        ),
    )
    CONFIG.declare(
        "crusher_equipment",
        ConfigValue(
            default=None,
            domain=In([None] + list(CRUSHER_EQUIPMENT_TO_STAGE.keys())),
            description="Optional crusher equipment name used to infer stage",
            doc="""Optional equipment selector.

Valid values:
- None
- jaw
- gyratory
- cone
- roll1
- roll2
- short_head_cone
- hammer_mill
""",
        ),
    )
    CONFIG.declare(
        "enforce_stage_size_limits",
        ConfigValue(
            default=True,
            domain=In([True, False]),
            description="Whether to enforce stage-specific product size constraints",
            doc="""If True, enforce inclusive interval constraints for the selected
stage. Note that the stated primary-stage classification is product P80 > 10 cm;
for the algebraic constraint, the lower bound is enforced as product P80 >= 10 cm
because strict inequalities cannot be imposed directly in Pyomo constraints.""",
        ),
    )

    def _resolve_stage_selection(self):
        selected_stage = self.config.crusher_stage
        selected_equipment = self.config.crusher_equipment

        # Return values are:
        # (resolved_stage, resolved_equipment, selection_basis)
        if selected_stage is None and selected_equipment is None:
            return "secondary", None, "default"

        if selected_stage is not None and selected_equipment is None:
            return selected_stage, None, "stage"

        selected_equipment = selected_equipment.strip().lower()
        inferred_stage = CRUSHER_EQUIPMENT_TO_STAGE[selected_equipment]

        if selected_stage is None:
            return inferred_stage, selected_equipment, "equipment"

        if selected_stage != inferred_stage:
            raise ValueError(
                "Inconsistent crusher specification: crusher_stage="
                f"'{selected_stage}' conflicts with crusher_equipment="
                f"'{selected_equipment}', which `{selected_equipment}`"
                f" belongs to stage '{inferred_stage}'."
            )

        return selected_stage, selected_equipment, "stage_and_equipment"


    def build(self):
        """
        Begin building model (pre-DAE transformation).

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super().build()

        # Build state blocks.
        self.properties_in = self.config.property_package.build_state_block(
            self.flowsheet().time,
            defined_state=True,
            **self.config.property_package_args,
        )
        self.properties_out = self.config.property_package.build_state_block(
            self.flowsheet().time,
            defined_state=True,
            **self.config.property_package_args,
        )

        tref = self.flowsheet().time.first()
        statevars = self.properties_in[tref].define_state_vars()

        # A constraint of In=Out for all state variables which are not related to particle size.
        for k, v in statevars.items():
            if k not in ["particle_size_median", "particle_size_width"]:
                idx = v.index_set()
                c = Constraint(
                    self.flowsheet().time,
                    idx,
                    doc=f"{k} constraint",
                    rule=partial(_state_rule, state=k),
                )
                self.add_component(k + "_constraint", c)

        # Add Ports
        self.add_port("inlet", self.properties_in)
        self.add_port("outlet", self.properties_out)

        self.work = Var(
            self.flowsheet().time,
            units=pyunits.W,
            initialize=3915.17,
            doc="Work required to increase crush the solid",
        )

        """ Breakage Distribution calculation as a constraint. 
        This is the equation for accumulative fraction of solid breakage 
        probability distribution smaller than size x=feed_p80
        """

        sunit = self.properties_in[tref].particle_size_median.get_units()
        stage_name, equipment_name, selection_basis = self._resolve_stage_selection()
        stage_spec = CRUSHER_STAGE_DATA[stage_name]

        self._crusher_stage = stage_name
        self._crusher_equipment = equipment_name
        self._stage_selection_basis = selection_basis
        self._stage_description = stage_spec["description"]
        self._stage_equipment_names = stage_spec["equipment_names"]
        self._stage_equipment = _equipment_label(stage_spec["equipment_names"])
        self._applicable_product_p80_lower = _range_bound_to_native_units(
            stage_spec["applicable_product_p80_lower"], sunit
        )
        self._applicable_product_p80_upper = _range_bound_to_native_units(
            stage_spec["applicable_product_p80_upper"], sunit
        )
        self._stage_product_p80_eps = _stage_eps_to_native_units(sunit)

        default_product_p80 = _size_to_native_units(
            stage_spec["default_product_p80"], sunit
        )

        for t in self.flowsheet().time:
            width_guess = value(self.properties_out[t].particle_size_width)
            p80_factor = (-math.log(1 - 0.8)) ** (width_guess / 2)
            if (not self.properties_out[t].particle_size_median.fixed) and p80_factor > 0:
                self.properties_out[t].particle_size_median.set_value(
                    default_product_p80 / p80_factor
                )
        
        @self.Expression(self.flowsheet().time, doc="Feed P80 size")
        def feed_p80(self, t):
            return (
                self.properties_in[t].particle_size_median
                / sunit
                * (-log(1 - 0.8)) ** (self.properties_in[t].particle_size_width / 2)
            )

        @self.Expression(self.flowsheet().time, doc="Product P80 size")
        def prod_p80(self, t):
            return (
                self.properties_out[t].particle_size_median
                / sunit
                * (-log(1 - 0.8)) ** (self.properties_out[t].particle_size_width / 2)
            )

        @self.Constraint(self.flowsheet().time, doc="Crusher work constraint")
        def crush_work_eq(self, t):
            return self.work[t] == (
                10  # 10 is an empirical correlation, this should not be changed.
                * self.properties_in[t].flow_mass
                * self.config.property_package.bond_work_index
                * (1 / (self.prod_p80[t]) ** 0.5 - 1 / (self.feed_p80[t]) ** 0.5)
            )
        
        @self.Constraint(
            self.flowsheet().time,
            doc="Crusher must reduce or maintain particle size (feed P80 >= product P80)",
        )
        def crushing_direction_constraint(self, t):
            return self.feed_p80[t] >= self.prod_p80[t]

        @self.Constraint(
            self.flowsheet().time,
            doc="Selected stage lower product-P80 bound",
        )
        def product_p80_lower_bound(self, t):
            if not self.config.enforce_stage_size_limits:
                return Constraint.Skip
            if self._applicable_product_p80_lower is None:
                return Constraint.Skip
            return self.prod_p80[t] >= self._applicable_product_p80_lower

        @self.Constraint(
            self.flowsheet().time,
            doc="Selected stage upper product-P80 bound",
        )
        def product_p80_upper_bound(self, t):
            if not self.config.enforce_stage_size_limits:
                return Constraint.Skip
            if self._applicable_product_p80_upper is None:
                return Constraint.Skip
            return self.prod_p80[t] <= self._applicable_product_p80_upper - self._stage_product_p80_eps

    def recommend_stage_for_product_p80(self, t=None):
        """
        Recommend a crusher stage for the current product P80.

        Parameters
        ----------
        t : hashable or None
            Time index to inspect. If None, use the first time point.

        Returns
        -------
        tuple[str, bool]
            (recommended_stage, in_modeled_range)
        """
        if t is None:
            t = self.flowsheet().time.first()

        tref = self.flowsheet().time.first()
        sunit = self.properties_in[tref].particle_size_median.get_units()
        product_p80_native = value(self.prod_p80[t])
        return _recommend_stage_from_product_p80_native(product_p80_native, sunit)

    def check_applicability(self, time_points=None, emit_warning=True):
        """
        Check whether the current product P80 is consistent with the configured
        stage classification interval.

        Parameters
        ----------
        time_points : iterable or None
            Time points to check. If None, check all model time points.
        emit_warning : bool
            If True, log a warning for each issue.

        Returns
        -------
        list[str]
            A list of warning messages. Empty if no issues are found.
        """
        if time_points is None:
            time_points = list(self.flowsheet().time)

        tref = self.flowsheet().time.first()
        sunit = self.properties_in[tref].particle_size_median.get_units()

        lower_cm = None
        if self._applicable_product_p80_lower is not None:
            lower_cm = value(
                pyunits.convert(
                    self._applicable_product_p80_lower * sunit, to_units=pyunits.cm
                )
                / pyunits.cm
            )

        upper_cm = None
        if self._applicable_product_p80_upper is not None:
            upper_cm = value(
                pyunits.convert(
                    self._applicable_product_p80_upper * sunit, to_units=pyunits.cm
                )
                / pyunits.cm
            )

        warning_messages = []
        for t in time_points:
            try:
                product_p80_native = value(self.prod_p80[t])
            except Exception:  # pragma: no cover - defensive fallback
                continue

            if product_p80_native is None:
                continue

            product_p80_cm = value(
                pyunits.convert(product_p80_native * sunit, to_units=pyunits.cm)
                / pyunits.cm
            )

            recommended_stage, in_modeled_range = _recommend_stage_from_product_p80_native(
                product_p80_native, sunit
            )
            recommended_equipment = _equipment_label(
                CRUSHER_STAGE_DATA[recommended_stage]["equipment_names"]
            )

            equipment_fragment = ""
            if self._crusher_equipment is not None:
                equipment_fragment = f", equipment='{self._crusher_equipment}'"

            in_selected_stage_range = _product_p80_in_stage_range(
                self._crusher_stage, product_p80_native, sunit
            )

            if not in_modeled_range:
                message = (
                    f"Crusher stage '{self._crusher_stage}'{equipment_fragment} may be outside "
                    f"its modeled applicability at time {t}: product P80 = "
                    f"{product_p80_cm:.6g} cm is finer than the minimum modeled crushing "
                    f"range (0.5 cm). The finest configured stage is '{recommended_stage}' "
                    f"({recommended_equipment})."
                )
                warning_messages.append(message)
                if emit_warning:
                    _log.warning(message)
                continue

            if not in_selected_stage_range:
                if upper_cm is None:
                    selected_range_text = f">= {lower_cm:.6g} cm"
                else:
                    selected_range_text = f"[{lower_cm:.6g}, {upper_cm:.6g}) cm"

                message = (
                    f"Crusher stage '{self._crusher_stage}'{equipment_fragment} may not be the "
                    f"most appropriate category at time {t}: product P80 = {product_p80_cm:.6g} "
                    f"cm falls in the '{recommended_stage}' crushing range "
                    f"({recommended_equipment}), while the current stage"
                    f" range is {selected_range_text}."
                )
                warning_messages.append(message)
                if emit_warning:
                    _log.warning(message)

        return warning_messages

    def _get_stream_table_contents(self, time_point=0):
        # Dictionary to hold data for all streams
        io_dict = {"Inlet": self.properties_in, "Outlet": self.properties_out}
        return create_stream_table_dataframe(io_dict, time_point=time_point)

    def _get_performance_contents(self, time_point=0):
        # Report
        var_dict = {
            "Work Required (W)": value(self.work[time_point]),
            "Feed P80": value(self.feed_p80[time_point]),
            "Product P80": value(self.prod_p80[time_point]),
        }
        return {"vars": var_dict}

def _state_rule(b, time, index, state):
    sin = b.properties_in[time].define_state_vars()[state]
    sout = b.properties_out[time].define_state_vars()[state]
    return sin[index] == sout[index]
