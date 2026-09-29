#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
whims_separator.py

Wet High-Intensity Magnetic Separator (WHIMS / HGMS) unit model built on top
of :mod:`mineral_size_psd`, the joint (size, mineral) PSD property package.

Total recovery is the sum of a magnetically-captured fraction R_M and a
physically-entrapped fraction R_P (Dobby & Finch's two-mechanism framing;
R_M = R_T - R_P, i.e. forward: R_T = R_M + R_P):

    R_T[k, j] = clip[ R_M[k, j] + Rp[k],  0, 1 ]

R_M follows King [1], chapter 8, eq. 8.53 -- clipped linear-in-log10 of a
magnetic force group M relative to a cut point M50:

    R_M = clip[ 0.5 + B * log10(M / M50),  0, 1 ]                     (8.53)

**Group M's exact functional form is STILL NOT independently verified.** The
"Model theory" section of the Met Dynamics wiki page (which would give the
precise combination of H, G [field gradient, T], Hs, rho, chi, beta, d, u,
Lm into M) remains behind a login I don't have legitimate access to
What's implemented below for M PHYSICALLY MOTIVATED reconstruction using only the
confirmed public inputs:

    M = H_eff * G * (rho*chi)^beta * d_char^2.5 / (u^1.8 * Lm^0.8)

  where H_eff = H*Hs/(H+Hs) is a smooth saturation-limited field (-> H for
  H << Hs, -> Hs for H >> Hs -- my own construction, not sourced), G enters
  multiplicatively as the field-gradient term (physically consistent with
  HGMS capture-force theory, where force ~ H * grad(H)), and beta (now a
  genuine calibratable Param, not fixed) replaces the fixed 1.2 exponent
  King's simplified eq. 8.52 used on (rho*chi). The d_char^2.5, u^1.8,
  Lm^0.8 exponents are carried over from King's eq. 8.52 since that's the
  only sourced numeric value I have for them; the full model may use
  different exponents once G and Hs are included explicitly -- treat these
  as placeholders to recalibrate, not literature-verified constants.

**UNIT-CONVENTION CAVEAT.** The public Excel/SysCAD parameter tables give
Met Dynamics' actual working units: H, Hs, G in Tesla; u in m/s; Lm in
kg/kg; M50, B, beta dimensionless; particle size in **mm** (not m); feed
mass flow in **t/h**; density in **t/m3**. This property package
(`mineral_size_psd`) and this unit model both use SI throughout (m, kg/s,
kg/m3), so **Met Dynamics report feed rates and densities are in t/h and t/m3 and must be converted to kg/s and kg/m3** before
comparing against or fixing this model's state. Any B/M50/G values
calibrated from a Met Dynamics/Excel run are tied to that mm/t·h-1/t·m-3
convention and are not directly portable into this SI implementation
without either explicit unit conversion at the calibration step or a
straight refit in SI.

The numerical example data used by the corresponding test are transcribed
from the Met Dynamics WHIMS example [3]. Its Model theory section is
login-gated; the exact vendor definition of M was not available.

REFERENCES:
  [1] King, R. P. "8 - Magnetic Separation." In Modeling and Simulation of
      Mineral Processing Systems. Butterworth-Heinemann, 2001.
      https://doi.org/10.1016/B978-0-08-051184-9.50012-2
  [2] Dobby, G. and J. A. Finch. "Capture of Mineral Particles in a High
      Gradient Magnetic Field." Powder Technology 17, no. 1 (1977): 73-82.
      https://doi.org/10.1016/0032-5910(77)85044-4
  [3] Met Dynamics. "WHIMS (Dobby and Finch)."
      https://wiki.metdynamics.com.au/view/WHIMS_(Dobby_and_Finch)

Author: Carolina Tristan
"""

from pyomo.common.config import ConfigBlock, ConfigValue
from pyomo.environ import Param, Var, log, log10, units as pyunits

from idaes.core import UnitModelBlockData, declare_process_block_class, useDefault
from idaes.core.initialization import InitializerBase
from idaes.core.scaling import CustomScalerBase
from idaes.core.util.config import is_physical_parameter_block
import idaes.logger as idaeslog

from prommis.psd.core.smooth_functions import EPS_SMOOTH, smooth_clamp
from prommis.psd.core.unit_utils import declare_factor_source, required_scaling_factor

__author__ = "Carolina Tristan"
__all__ = ["WHIMSSeparator"]

_log = idaeslog.getLogger(__name__)

# Exponents carried over from King [1], eq. 8.52 -- see the module
# docstring's caveat: sourced for a SIMPLIFIED version of the correlation
# that omits G and Hs, so treat these as placeholders pending recalibration
# against own test data, not verified constants of the full model.
_EXP_DP = 2.5
_EXP_U = 1.8
_EXP_LM = 0.8

# Numerical floor inside log10() to keep R_M's equation finite as M -> 0
# (e.g. a fully non-magnetic gangue mineral, rho*chi -> 0) during
# initialization from a flat starting point.
_EPS_LOG = 1e-12

class WHIMSSeparatorScaler(CustomScalerBase):
    """Scaler for the WHIMS separator.

    Delegates variable scaling for the three attached state blocks to each
    block's own default scaler (``MineralSizePSDScaler``) via
    ``call_submodel_scaler_method``, then derives factors for the mass
    balance Constraints from those declared flow factors via
    ``required_scaling_factor``.
    """

    CONFIG = declare_factor_source(CustomScalerBase.CONFIG())

    RECOVERY_SCALING_FACTOR = 1.0

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        for state_name in ("feed_state", "mags_state", "nonmags_state"):
            self.call_submodel_scaler_method(
                model,
                submodel=getattr(model, state_name),
                method="variable_scaling_routine",
                submodel_scalers=submodel_scalers,
                overwrite=overwrite,
            )
        for t in model.flowsheet().time:
            for k in model.config.property_package.size_interval_set:
                for j in model.config.property_package.solid_component_set:
                    self.set_variable_scaling_factor(
                        model.recovery[t, k, j],
                        self.RECOVERY_SCALING_FACTOR,
                        overwrite=overwrite,
                    )
                    self.set_variable_scaling_factor(
                        model.recovery_M[t, k, j],
                        self.RECOVERY_SCALING_FACTOR,
                        overwrite=overwrite,
                    )

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        tset = model.flowsheet().time
        pp = model.config.property_package
        for t in tset:
            for k in pp.size_interval_set:
                for j in pp.solid_component_set:
                    feed_var = model.feed_state[t].flow_mass_size_comp[k, j]
                    sf = required_scaling_factor(self, feed_var)
                    self.set_constraint_scaling_factor(
                        model.mags_solid_balance[t, k, j], sf, overwrite=overwrite
                    )
                    self.set_constraint_scaling_factor(
                        model.nonmags_solid_balance[t, k, j], sf, overwrite=overwrite
                    )
                    self.set_constraint_scaling_factor(
                        model.recovery_M_eqn[t, k, j],
                        1.0 / self.RECOVERY_SCALING_FACTOR,
                        overwrite=overwrite,
                    )
                    self.set_constraint_scaling_factor(
                        model.recovery_eqn[t, k, j],
                        1.0 / self.RECOVERY_SCALING_FACTOR,
                        overwrite=overwrite,
                    )


class WHIMSSeparatorInitializer(InitializerBase):
    """Initializer for the WHIMS separator.

    Both recovery equations depend only on fixed parameters (chi, rho, H,
    Hs, G, beta, u, Lm, M50, B, Rp) -- not on any flow -- so they're
    evaluated directly rather than solved iteratively:

      1. Initialize ``feed_state`` via its own default initializer.
      2. Evaluate ``M``, then ``recovery_M`` (eq. 8.53), then total
         ``recovery = clip(recovery_M + Rp, 0, 1)`` directly from current
         parameter values for every (size, mineral) cell.
      3. Propagate ``mags_state`` / ``nonmags_state`` from
         ``recovery * feed_state`` directly, then hand off to the caller's
         solve of the full square system.
    """

    CONFIG = InitializerBase.CONFIG()

    def initialization_routine(self, model):
        from pyomo.environ import value

        pp = model.config.property_package
        tset = model.flowsheet().time

        model.feed_state.default_initializer().initialize(model.feed_state)

        for t in tset:
            for k in pp.size_interval_set:
                for j in pp.solid_component_set:
                    # M = value(model.M[t, k, j])
                    # B = value(model.B)
                    # M50 = value(model.M50)
                    Rp = value(model.Rp[k])

                    r_m_raw = value(
                        0.5 + model.B * log10((model.M[t, k, j] + _EPS_LOG) / model.M50)
                    )
                    r_m = min(1.0, max(0.0, r_m_raw))
                    model.recovery_M[t, k, j].set_value(r_m)

                    r_total = min(1.0, max(0.0, r_m + Rp))
                    model.recovery[t, k, j].set_value(r_total)

                    feed_val = value(model.feed_state[t].flow_mass_size_comp[k, j])
                    model.mags_state[t].flow_mass_size_comp[k, j].set_value(
                        r_total * feed_val
                    )
                    model.nonmags_state[t].flow_mass_size_comp[k, j].set_value(
                        (1 - r_total) * feed_val
                    )
        return None


@declare_process_block_class("WHIMSSeparator")
class WHIMSSeparatorData(UnitModelBlockData):
    """
    Wet High-Intensity Magnetic Separator (WHIMS) unit model.

    Splits a mineral-by-size solid feed stream (a
    :class:`mineral_size_psd.MineralSizePSDStateBlock`) into a Magnetics
    stream and a Non-Magnetics stream, using the Dobby & Finch recovery
    correlation (magnetic capture, eq. 8.53) plus a physical-entrapment
    term Rp, combined additively (R_T = R_M + Rp), computed independently
    for every (size interval, mineral) cell of the joint
    ``flow_mass_size_comp`` state.

    Ports:
        feed          (inlet)
        magnetics     (outlet)
        non_magnetics (outlet)
    """

    default_scaler = WHIMSSeparatorScaler
    default_initializer = WHIMSSeparatorInitializer

    # ------------------------------------------------------------------
    # CONFIG
    # ------------------------------------------------------------------
    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Mineral-by-size PSD property package "
            "(mineral_size_psd.MineralSizePSDParameterBlock)",
            doc="Must expose `size_interval_set`, `solid_component_set`, and "
            "`d_char[k]`, and build state blocks with a joint "
            "`flow_mass_size_comp[k, j]` state Var.",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
        ),
    )

    # ------------------------------------------------------------------
    # BUILD
    # ------------------------------------------------------------------
    def build(self):
        super().build()

        pp = self.config.property_package
        size_set = pp.size_interval_set
        mineral_set = pp.solid_component_set
        tset = self.flowsheet().time

        # --------------------------------------------------------------
        # State blocks / ports.
        # --------------------------------------------------------------
        self.feed_state = pp.build_state_block(
            tset, defined_state=True, **self.config.property_package_args
        )
        self.mags_state = pp.build_state_block(
            tset, defined_state=False, **self.config.property_package_args
        )
        self.nonmags_state = pp.build_state_block(
            tset, defined_state=False, **self.config.property_package_args
        )

        self.add_port(name="feed", block=self.feed_state, doc="Feed solids")
        self.add_port(name="magnetics", block=self.mags_state, doc="Magnetics outlet")
        self.add_port(
            name="non_magnetics",
            block=self.nonmags_state,
            doc="Non-magnetics outlet",
        )

        # --------------------------------------------------------------
        # Reference scales for non-dimensionalizing M (see UNIT-CONVENTION
        # CAVEAT). Defaulted to 1 in SI units.
        # --------------------------------------------------------------
        self.H_ref = Param(
            initialize=1.0,
            mutable=True,
            units=pyunits.T,
            doc="Reference field strength for non-dimensionalizing H.",
        )
        self.d_ref = Param(
            initialize=1.0,
            mutable=True,
            units=pyunits.m,
            doc="Reference length for non-dimensionalizing d_p.",
        )
        self.u_ref = Param(
            initialize=1.0,
            mutable=True,
            units=pyunits.m / pyunits.s,
            doc="Reference velocity for non-dimensionalizing u.",
        )

        # --------------------------------------------------------------
        # WHIMS operating parameters -- reinstated H, Hs, beta, G, and
        # per-size Rp, all confirmed as real Met Dynamics WHIMS-page /
        # Susceptibility-page inputs.
        # --------------------------------------------------------------
        self.H = Param(
            initialize=1.0,
            mutable=True,
            units=pyunits.T,
            doc="Applied magnetic field strength, H",
        )
        self.Hs = Param(
            initialize=2.0,
            mutable=True,
            units=pyunits.T,
            doc="Saturation magnetic field strength, Hs",
        )
        self.G = Param(
            initialize=1.0,
            mutable=True,
            units=pyunits.T,
            doc="Magnetic field gradient, G (T)",
        )
        self.u = Param(
            initialize=0.05,
            mutable=True,
            units=pyunits.m / pyunits.s,
            doc="Average slurry velocity through the matrix, u",
        )
        self.Lm = Param(
            initialize=0.05,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Matrix loading Lm = ratio of mass of magnetic solids to "
            "mass of matrix in the separation chamber (kg/kg).",
        )
        self.M50 = Param(
            initialize=1.269e-2,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Magnetic cut point. Literature example value -- see "
            "UNIT-CONVENTION CAVEAT before trusting this number here.",
        )
        self.B = Param(
            initialize=0.348,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Recovery-curve slope (coefficient of the log term). "
            "Literature example value -- see UNIT-CONVENTION CAVEAT.",
        )
        self.beta = Param(
            initialize=1.2,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Susceptibility exponent, confirmed real input "
            "(SysCAD 'SusceptibilityExponent / Beta'). Default 1.2 "
            "carried over from King's simplified fixed-exponent example "
            "-- calibrate against your own test data.",
        )

        # --------------------------------------------------------------
        # Per (size, mineral) susceptibility and density, and per-size
        # physical-entrapment recovery Rp -- all confirmed inputs on the
        # public SysCAD Susceptibility page.
        # --------------------------------------------------------------
        self.chi = Param(
            size_set,
            mineral_set,
            initialize=1e-7,
            mutable=True,
            units=pyunits.m**3 / pyunits.kg,
            doc="Mass magnetic susceptibility, chi[k, j]",
        )
        self.rho = Param(
            size_set,
            mineral_set,
            initialize=2700.0,
            mutable=True,
            units=pyunits.kg / pyunits.m**3,
            doc="Mineral density, rho[k, j] (kg/m3) -- Met Dynamics reports "
            "this in t/m3; convert before setting.",
        )
        self.Rp = Param(
            size_set,
            initialize=0.0,
            mutable=True,
            # bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Fraction recovered to magnetics by physical entrapment, "
            "Rp[k] -- a real per-size-interval input (SysCAD Susceptibility "
            "page), not a fabricated term.",
        )

        # --------------------------------------------------------------
        # Recovery variables: magnetic-capture-only R_M, and total R_T
        # (named `recovery` for continuity with earlier versions).
        # --------------------------------------------------------------
        self.recovery_M = Var(
            tset,
            size_set,
            mineral_set,
            initialize=0.5,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Magnetically-captured recovery fraction, R_M (eq. 8.53), "
            "per (size interval, mineral).",
        )
        self.recovery = Var(
            tset,
            size_set,
            mineral_set,
            initialize=0.5,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Total recovery to magnetics, R_T = clip(R_M + Rp, 0, 1), "
            "per (size interval, mineral).",
        )

        def _field_function(b):
            # Smooth saturation-limited field: -> H for H << Hs, -> Hs for
            # H >> Hs. Own construction (see module docstring caveat).
            return b.H * b.Hs / (b.H + b.Hs)

        @self.Expression(
            tset,
            size_set,
            mineral_set,
            doc="Magnetic force group M -- best-effort reconstruction, "
            "see module docstring caveat on exact functional form.",
        )
        def M(b, t, k, j):
            H_star = _field_function(b) / b.H_ref
            d_star = pp.d_char[k] / b.d_ref
            u_star = b.u / b.u_ref
            rhochi = b.rho[k, j] * b.chi[k, j]  # dimensionless: kg/m3 * m3/kg
            return (
                H_star
                * b.G
                * rhochi**b.beta
                * d_star**_EXP_DP
                / (u_star**_EXP_U * b.Lm**_EXP_LM)
            )

        # --------------------------------------------------------------
        # R_M, eq. 8.53 -- clipped linear-in-log10.
        # --------------------------------------------------------------
        @self.Constraint(
            tset,
            size_set,
            mineral_set,
            doc="Magnetic-capture recovery, eq. 8.53",
        )
        def recovery_M_eqn(b, t, k, j):
            r_raw = 0.5 + b.B * log10((b.M[t, k, j] + _EPS_LOG) / b.M50)
            return b.recovery_M[t, k, j] == smooth_clamp(
                r_raw, 0.0, 1.0, eps=EPS_SMOOTH
            )

        # --------------------------------------------------------------
        # R_T = clip(R_M + Rp, 0, 1) -- confirmed additive combination.
        # --------------------------------------------------------------
        @self.Constraint(
            tset,
            size_set,
            mineral_set,
            doc="Total recovery R_T = clip(R_M + Rp, 0, 1)",
        )
        def recovery_eqn(b, t, k, j):
            return b.recovery[t, k, j] == smooth_clamp(
                b.recovery_M[t, k, j] + b.Rp[k], 0.0, 1.0, eps=EPS_SMOOTH
            )

        # --------------------------------------------------------------
        # Mass balances on the single joint state variable.
        # --------------------------------------------------------------
        @self.Constraint(
            tset,
            size_set,
            mineral_set,
            doc="Magnetics solid mass flow per (size, mineral)",
        )
        def mags_solid_balance(b, t, k, j):
            return b.mags_state[t].flow_mass_size_comp[k, j] == (
                b.recovery[t, k, j] * b.feed_state[t].flow_mass_size_comp[k, j]
            )

        @self.Constraint(
            tset,
            size_set,
            mineral_set,
            doc="Non-magnetics solid mass flow per (size, mineral)",
        )
        def nonmags_solid_balance(b, t, k, j):
            return (
                b.nonmags_state[t].flow_mass_size_comp[k, j]
                == (1 - b.recovery[t, k, j]) * b.feed_state[t].flow_mass_size_comp[k, j]
            )

        # --------------------------------------------------------------
        # Convenience: overall (mass-weighted) recovery.
        # --------------------------------------------------------------
        @self.Expression(tset, doc="Overall (mass-weighted) total recovery")
        def overall_recovery(b, t):
            total_in = b.feed_state[t].flow_mass
            total_mags = sum(
                b.recovery[t, k, j] * b.feed_state[t].flow_mass_size_comp[k, j]
                for k in size_set
                for j in mineral_set
            )
            return total_mags / (total_in + pp.flow_eps)
