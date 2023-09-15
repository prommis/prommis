import os
import json
import enum

import pyomo.environ as pyo
from pyomo.common.config import ConfigValue, In
from pyomo.common.fileutils import this_file_dir

from idaes.core import (
    StateBlock,
    StateBlockData,
    PhysicalParameterBlock,
    LiquidPhase,
    VaporPhase,
    Phase,
    Component,
    ElectrolytePropertySet,
    MaterialBalanceType,
    MaterialFlowBasis,
)
from idaes.core import declare_process_block_class
from idaes.core.util.math import smooth_min, smooth_abs


_default_units = {
    "time": pyo.units.s,
    "length": pyo.units.m,
    "mass": pyo.units.kg,
    "amount": pyo.units.mol,
    "temperature": pyo.units.K,
}


def get_data_path():
    return os.path.join(this_file_dir(), "precip_prop.json")


def read_aqueous_equilibrium_parameters(data_file, key_components, skip_set=None):
    """Read component data from json a database file.  Limit the components
    included to the ones selected based on the key components and solids. In
    addition to key components, H2O, H^+, and OH^-, can be used in reactions.
    The database also includes parameters for equilibrium constants and activity
    coefficients.
    """
    if skip_set is None:
        skip_set = set()

    with open(data_file, "r") as fp:
        dat = json.load(fp)
    key_dict = {}
    aq_dict = {}
    solid_dict = {}

    # In addition to key components, H2O, H^+, and OH^-, can be used in reactions

    # Selected key components plus water
    keys_plus = key_components | {"H2O", "H^+", "OH^-"}

    # All keys defined in the database
    all_plus = set(dat["key_components"].keys()) | {"H2O"}

    # selected key component data dictionary, if selecting a key component
    # that is not defined in the database, a key error will be raised.
    # force H^+ in key components,
    # OH^- is not key but convenient to use in reactions, so will remove later
    key_components = set(key_components) | {"H^+", "OH^-"}
    key_dict = {k: dat["key_components"][k] for k in key_components}

    # Verify DB no misspellings in stoich and using key components only
    for k, d in dat["aq_components"].items():
        for k2 in d["stoich"]:
            if k2 not in all_plus:
                raise KeyError(f"{k2} is not a key component")
    for k, d in dat["solid_components"].items():
        for k2 in d["stoich"]:
            if k2 not in all_plus:
                raise KeyError(f"{k2} is not a key component")

    # Select aqueous components that can be formed from selected keys
    for k, d in dat["aq_components"].items():
        if k in skip_set:
            keep = False
        else:
            keep = True
            if k not in ["OH^-"]:  # "OH^-" automatically
                for k2 in d["stoich"]:
                    if k2 not in keys_plus:
                        keep = False
                        break
        if keep:
            aq_dict[k] = d

    # Select solid components that can be formed from selected keys
    for k, d in dat["solid_components"].items():
        if k in skip_set:
            keep = False
        else:
            keep = True
            for k2 in d["stoich"]:
                if k2 not in keys_plus:
                    keep = False
                    break
        if keep:
            solid_dict[k] = d

    # Get all molecular weights
    mw_dict = {"H2O": dat["H2O"]["mw"] / 1000.0}
    for k, d in key_dict.items():
        mw_dict[k] = d["mw"] / 1000.0
    for k, d in aq_dict.items():
        mw_dict[k] = 0
        for ks, scoeff in d["stoich"].items():
            mw_dict[k] -= mw_dict[ks] * scoeff
    for k, d in solid_dict.items():
        mw_dict[k] = 0
        for ks, scoeff in d["stoich"].items():
            mw_dict[k] -= mw_dict[ks] * scoeff

    # OH isn't really a key comp, but want to use it in reactions anyway
    # so take it out of the key list.
    del key_dict["OH^-"]

    return key_dict, aq_dict, solid_dict, mw_dict


def build_parameter_block(blk):
    """Although the solids property package and the aqueous property package
    do not track the same components, we'll use the same parameter block for
    both to simplify things.

    Args:
        blk (Block): add parameters to this block

    Returns:
        None

    """
    data_file = get_data_path()
    (
        blk.key_dict,
        blk.aq_dict,
        blk.solid_dict,
        mw_dict,
    ) = read_aqueous_equilibrium_parameters(
        data_file, blk.config.key_components, blk.config.skip_components
    )
    blk.mw_comp = pyo.Param(
        set(mw_dict.keys()),
        initialize=mw_dict,
        doc="Component molecular weight",
        units=pyo.units.kg / pyo.units.mol,
    )

    blk.key_set = pyo.Set(initialize=list(blk.key_dict.keys()))
    blk.aq_set = pyo.Set(initialize=list(blk.aq_dict.keys()))
    blk.solid_set = pyo.Set(initialize=list(blk.solid_dict.keys()))

    blk.charge_dict = {}
    for k, d in blk.key_dict.items():
        blk.charge_dict[k] = d["charge"]
    for k, d in blk.aq_dict.items():
        blk.charge_dict[k] = d["charge"]

    blk.charge = pyo.Param(
        blk.aq_set | blk.key_set,
        initialize=blk.charge_dict,
        doc="Component Charge",
    )

    blk.a_act = pyo.Param(
        initialize=0.5108,
        mutable=True,
        doc="A parameter in the Davies activity coefficient equation",
        units=pyo.units.mol / pyo.units.kg,
    )
    blk.b_act = pyo.Param(
        initialize=0.3,
        mutable=True,
        doc="B parameter in the Davies activity coefficient equation",
    )


def _config_blk_build(blk):
    blk.declare(
        "key_components",
        ConfigValue(
            default=None,
            domain=set,
            description="Set of key components",
            doc="(set) Set of key components",
        ),
    )

    blk.declare(
        "skip_components",
        ConfigValue(
            default=None,
            domain=set,
            description="Set of reactions to ignore",
            doc="(set) Set of reactions to ignore.  This can be used to"
            " simplify equations by ignoring components that are not"
            " expected to form.",
        ),
    )


class ActivityCoefficientType(enum.Enum):
    """This Enum is used to set the activity coefficient calculation type.

    * CONSTANT: User settable parameter, 1 by default
    * DAVIES: Calculate using Davies equation best for low ionic strength
    * WATEQ: Calculate using WATEQ equation
    * SIT: Calculate using SIT equation
    """

    CONSTANT = 1
    DAVIES = 2
    WATEQ = 3
    SIT = 4


@declare_process_block_class("AqueousStateParameterBlock")
class AqueousStateParameterBlockData(PhysicalParameterBlock):
    CONFIG = PhysicalParameterBlock.CONFIG()
    _config_blk_build(CONFIG)

    def build(self):
        super().build()
        build_parameter_block(self)
        self._state_block_class = AqueousStateBlock
        self.AqueousPhase = Phase()
        self.phase_list = ["AqueousPhase"]
        self.component_list = list(self.key_set | self.aq_set)
        print(self.get_metadata())

    @classmethod
    def define_metadata(cls, obj):
        obj.define_property_set(ElectrolytePropertySet)
        obj.add_properties(
            {
                "temperature": {"method": None, "units": "K"},
                "pressure": {"method": None, "units": "Pa"},
                "flow_mass": {"method": None, "units": "kg/s"},
                "log10_molality_comp": {
                    "method": None,
                    "units": "dimensionless",
                },
                "molality_comp": {"method": None, "units": "mol/kg"},
                "log10_act_coeff": {
                    "method": "_log10_act_coeff_method",
                    "units": "dimensionless",
                },
                "pH": {"method": "_pH_method", "units": "dimensionless"},
            }
        )
        obj.add_default_units(_default_units)


@declare_process_block_class("PrecipitateStateParameterBlock")
class PrecipitateStateParameterBlockData(PhysicalParameterBlock):
    CONFIG = PhysicalParameterBlock.CONFIG()
    _config_blk_build(CONFIG)

    def build(self):
        super().build()
        build_parameter_block(self)
        self._state_block_class = PrecipitateStateBlock
        self.SolidPhase = Phase()
        self.phase_list = ["SolidPhase"]
        self.component_list = list(self.solid_set)

    @classmethod
    def define_metadata(cls, obj):
        obj.define_property_set(ElectrolytePropertySet)
        obj.add_properties(
            {
                "temperature": {"method": None, "units": "K"},
                "flow_mol_comp": {"method": None, "units": "mol/s"},
            }
        )
        obj.add_default_units(_default_units)


class _AqueousStateBlock(StateBlock):
    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        # Fix state variables
        fix_state_vars(self)


class _PrecipitateStateBlock(StateBlock):
    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        # Fix state variables
        fix_state_vars(self)


@declare_process_block_class("AqueousStateBlock", block_class=_AqueousStateBlock)
class AqueousStateBlockData(StateBlockData):
    def build(self):
        super().build()
        params = self.params

        # State Variables
        self.temperature = pyo.Var(
            initialize=300,
            doc="Temperature",
            units=pyo.units.kelvin,
        )
        self.pressure = pyo.Var(
            initialize=101325,
            doc="Pressure",
            units=pyo.units.Pa,
        )
        self.flow_mass = pyo.Var(
            initialize=1, doc="Total mass flow rate", units=pyo.units.kg / pyo.units.s
        )

        self.log10_molality_comp = pyo.Var(
            params.key_set | params.aq_set,
            initialize=-15,
            doc="Log10 molality of components in solution and precipitates",
            units=pyo.units.dimensionless,
            bounds=(-50, 50),
        )

        @self.Expression(params.key_set | params.aq_set)
        def molality_comp(blk, comp):
            return 10 ** blk.log10_molality_comp[comp] * pyo.units.mol / pyo.units.kg 

        self.charge_concentration = pyo.Expression(
            expr=sum(
                self.molality_comp[i] * params.charge[i]
                for i in params.key_set | params.aq_set
                if params.charge[i] != 0
            ),
            doc="The charge molality of the solution.",
        )

        self.eps_unit = pyo.Param(
            initialize=1e-4,
            mutable=True,
            doc="A parameter in the Davies activity coefficient equation",
            units=pyo.units.mol / pyo.units.kg,
        )


        # Usually the ionic charge should be balanced, but in some cases it
        # may be convenient not to enforce a charge balance.  In the case
        # where charge is not balanced, we just assume that ions exist to
        # balance the charge, but they do not participate in reactions. The
        # unspecified ions are assumed to have a charge of -1 or 1.
        self.ionic_strength = pyo.Expression(
            expr=0.5
            * (
                smooth_abs(self.charge_concentration, eps=self.eps_unit)
                + sum(
                    self.molality_comp[i] * params.charge[i] ** 2
                    for i in params.key_set | params.aq_set
                    if params.charge[i] != 0
                )
            )
        )

    def _log10_act_coeff_method(self):
        """Add variable and equation for activity coefficients."""
        self.log10_act_coeff = pyo.Var(
            self.params.aq_set | self.params.key_set,
            initialize=0,
            doc="log10 activity coefficient",
            units=pyo.units.dimensionless,
        )

        @self.Constraint(
            self.params.aq_set | self.params.key_set,
            doc="Davies equation for activity coefficient.",
        )
        def log10_act_coeff_eqn(blk, comp):
            B = blk.params.b_act
            A = blk.params.a_act
            I = blk.ionic_strength
            sqrtI = pyo.sqrt(I)
            z = blk.params.charge[comp]
            return blk.log10_act_coeff[comp] == (B * I - A * z**2 * sqrtI / (1* pyo.units.mol**0.5 / pyo.units.kg**0.5 + sqrtI))/pyo.units.mol * pyo.units.kg

        """
        @self.Expression(
            self.params.aq_set | self.params.key_set,
            doc="Davies equation for activity coefficient.",
        )
        def log10_act_coeff(blk, comp):
            B = blk.params.b_act
            A = blk.params.a_act
            I = blk.ionic_strength
            sqrtI = pyo.sqrt(I)
            z = blk.params.charge[comp]
            return B * I - A * z**2 * sqrtI / (1 + sqrtI)
        """

    def _pH_method(self):
        self.pH = pyo.Var(
            initialize=7,
            doc="-log10 H^+ activity",
            units=pyo.units.dimensionless,
        )
        self.pH_eqn = pyo.Constraint(
            expr=-self.pH
            == self.log10_molality_comp["H^+"] + self.log10_act_coeff["H^+"]
        )

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mass

    def get_material_flow_terms(self, p, j):
        """Create material flow terms for control volume."""
        return self.log10_molality_comp[j]

    def default_material_balance_type(self):
        return MaterialBalanceType.componentTotal

    def define_state_vars(self):
        return {
            "flow_mass": self.flow_mass,
            "temperature": self.temperature,
            "pressure": self.pressure,
            "log10_molality_comp": self.log10_molality_comp,
        }

    def _flow_mol_comp(self):
        pass

    def _flow_mass_comp(self):
        pass


@declare_process_block_class(
    "PrecipitateStateBlock", block_class=_PrecipitateStateBlock
)
class PrecipitateStateBlockData(StateBlockData):
    def build(self):
        super().build()
        params = self.params

        # State Variables
        self.temperature = pyo.Var(
            initialize=300,
            doc="Temperature",
            units=pyo.units.kelvin,
        )
        self.flow_mol_comp = pyo.Var(
            params.solid_set,
            initialize=1,
            doc="Component molar flow rates",
            units=pyo.units.mol / pyo.units.s,
        )

    def define_state_vars(self):
        return {
            "flow_mass": self.flow_mol_comp,
            "temperature": self.temperature,
        }
