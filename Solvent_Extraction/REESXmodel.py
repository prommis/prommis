from pyomo.environ import Constraint, RangeSet, Reals, Var, Param
from pyomo.common.config import ConfigDict, ConfigValue, Bool, In
from pyomo.dae import DerivativeVar

from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
    FlowDirection,
    MaterialFlowBasis,
)
from idaes.core.util.config import (
    is_physical_parameter_block
)
from idaes.core.util.exceptions import (
    ConfigurationError,
    BurntToast,
)

Stream_Config = ConfigDict()

Stream_Config.declare(
    "property_package",
    ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for given stream",
        doc="""Property parameter object used to define property calculations for given stream,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}"""
    )
)

Stream_Config.declare(
    "property_package_args",
    ConfigDict(
        implicit=True,
        description="Dict of arguments to use for constructing property package",
        doc="""A ConfigDict with arguments to be passed to property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
    )
)

Stream_Config.declare(
    "flow_direction",
    ConfigValue(
        default=FlowDirection.forward,
        domain=In(FlowDirection),
        doc="Direction of flow for stream",
        description="FlowDirection Enum indicating direction of "
        "flow for given stream. Default=FlowDirection.forward.",
    )
)

Stream_Config.declare(
    "has_feed",
    ConfigValue(
        default=True,
        domain=Bool,
        doc="Bool indicating whether stream has a feed.",
        description="Bool indicating whether stream has a feed Port and inlet "
        "state, or if all flow is provided via mass transfer. Default=True.",
    )
)

Stream_Config.declare(
    "side_streams",
    ConfigValue(
        default=None,
        domain=list,
        doc="List of finite elements at which a side stream should be included.",
    )
)


@declare_process_block_class("REESX")
class REESXData(UnitModelBlockData):

    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "aqueous_streams",
        ConfigDict(
            implicit=True,
            implicit_domain=Stream_Config,
            description="Dict of aqueous streams and associated property packages",
            doc="ConfigDict with keys indicating names for each stream in system and "
            "values indicating property package and associated arguments.",
        )
    )

    CONFIG.declare(
        "organic_streams",
        ConfigDict(
            implicit=True,
            implicit_domain=Stream_Config,
            description="Dict of organic streams and associated property packages",
            doc="ConfigDict with keys indicating names for each stream in system and "
            "values indicating property package and associated arguments.",
        )
    )


    CONFIG.declare(
        "number_of_finite_elements",
        ConfigValue(domain=int, description="Number of finite elements to use"),
    )


    def build(self):

        super().build()

        self._verify_inputs()

        for stream, sconfig in self.config.aqueous_streams.items():
            aqflow_basis, aquom = self._build_state_blocks(stream)
        for stream, sconfig in self.config.organic_streams.items():
            ogflow_basis,oguom = self._build_state_blocks(stream)
        
        self._build_material_balance_constraints(aqflow_basis, aquom, ogflow_basis, oguom)
        self._build_ports()
    
   
    def _verify_inputs(self):
        # Check that at least two streams were declared
        if len(self.config.aqueous_streams) < 1:
            raise ConfigurationError(
                f"REESX models must define at least one stream; received "
                f"{list(self.config.aqueous_streams.keys())}"
            )
        
        if len(self.config.organic_streams) < 1:
            raise ConfigurationError(
                f"REESX models must define at least one stream; received "
                f"{list(self.config.organic_streams.keys())}"
            )
        
        self.elements = RangeSet(
            1,
            self.config.number_of_finite_elements,
            doc="Set of finite elements in cascade (1 to number of elements)",
        )
    
    # Creation of state blocks for a component

    def _build_state_blocks(self, stream):

        if stream in self.config.aqueous_streams.keys():
         streams = self.config.aqueous_streams
        elif stream in self.config.organic_streams.keys():
         streams = self.config.organic_streams
        else:
         raise BurntToast("If/else overrun when constructing balances")
        
        ppack = streams[stream].property_package

        arg_dict2 = dict(**streams[stream].property_package_args)
        arg_dict2["defined_state"] = False

        state = ppack.build_state_block(
                    self.flowsheet().time,
                    self.elements,
                    doc=f"States for stream {stream} in organic phase ",
                    **arg_dict2,
                )
            
        self.add_component(stream, state)

        if streams[stream].has_feed:
                arg_dict0 = dict(**streams[stream].property_package_args)
                arg_dict0["defined_state"] = True

                inlet_state = ppack.build_state_block(
                    self.flowsheet().time,
                    doc=f"Inlet states for stream {stream}.",
                    **arg_dict0,
                )
                self.add_component(stream + "_inlet_state", inlet_state)
        
        tref = self.flowsheet().time.first()
        sref = self.elements.first()
        
        flow_basis = state[tref, sref].get_material_flow_basis()
        uom = state[tref, sref].params.get_metadata().derived_units

        if flow_basis is None:
            # Set unit level flow basis and units from first stream

            flow_basis = state[tref, sref].get_material_flow_basis()
            uom = state[tref, sref].params.get_metadata().derived_units
        else:
            # Check that flow bases are consistent
            if not state[tref, sref].get_material_flow_basis() == flow_basis:
                raise ConfigurationError(
                    f"Property packages use different flow bases: ExtractionCascade "
                    f"requires all property packages to use the same basis. "
                    f"{stream} uses {state[tref, sref].get_material_flow_basis()}, "
                        f"whilst first stream uses {flow_basis}."
                    )

        return flow_basis, uom

    # Creation of material balance

    def _build_material_balance_constraints(self, aqflow_basis, aquom, ogflow_basis, oguom):
            # Get units for transfer terms
            if aqflow_basis is MaterialFlowBasis.molar:
                amb_units = aquom.FLOW_MOLE
            elif aqflow_basis is MaterialFlowBasis.mass:
                amb_units = aquom.FLOW_MASS
            else:
                # Flow type other, so cannot determine units
                amb_units = None
            

            for stream, sconfig in self.config.aqueous_streams.items():
                state_block = getattr(self, stream)
                ppack = sconfig.property_package
                component_list = state_block.component_list
                #component_list = self.parent_block().prop.dissolved_elements

                if stream in self.config.aqueous_streams.keys():
                    in_state = getattr(self, stream + "_inlet_state")
                
                #Dynamic part needs to be repaired

                #dC_a = Var(self.flowsheet().time,
                    #self.elements,
                    #ppack.dissolved_elements,
                    #domain=Reals,
                        #initialize=0.0)
                
                #self.add_component(
                        #stream + "dC_a",
                        #dC_a
                    #)
                
                #def dynamic_term_a(b,t,s,j):
                    #for j in ppack.dissolved_elements:
                       # if self.config.dynamic == True:
                           # return dC_a[t,s,j] == state_block[t, s].conc_mass_comp[j]
                        #else:
                            #return dC_a[t,s,j] == 0
                    #return Constraint.Skip
                
                #dynamic_term_a_constraint = Constraint(self.flowsheet().time, self.elements,
                                                            #ppack.dissolved_elements, rule=dynamic_term_a)
                
                #self.add_component(
                    #stream + "_dynamic_term_a_constraint",
                    #dynamic_term_a_constraint,
                #)

                #dC_dt_a = DerivativeVar(dC_a, wrt=self.flowsheet().time)

                #self.add_component(
                        #stream + "_dC_dt_a",
                        #dC_dt_a
                    #)

                # state_block.display()
                distribution_extent = Var(
                    self.flowsheet().time,
                    self.elements,
                    ppack.dissolved_elements,
                    domain=Reals,
                        initialize=0.0,
                        doc=f"Extent of transfer in stream {stream}",
                        units=amb_units,
                )
                self.add_component(
                        stream + "_distribution_extent",
                        distribution_extent
                    )
                
                partition_coefficient = Param(ppack.dissolved_elements, initialize = {
                    "Al":3.6/100,
                    "Ca":3.7/100,
                    "Fe":2.1/100,
                    "Si":0/100,
                    "Sc":100/100,
                    "Y":100/100,
                    "La":75.2/100,
                    "Ce":95.7/100,
                    "Pr":96.5/100,
                    "Nd":99.2/100,
                    "Pm":100/100,
                    "Sm":100/100,
                    "Eu":99.9/100,
                    "Gd":98.6/100,
                    "Tb":99.3/100,
                    "Dy":99.9/100,
                    "Ho":99.5/100,
                    "Er":99.5/100,
                    "Tm":98.6/100,
                    "Yb":80.7/100,
                    "Lu":99.5/100,
                    "Th":5/100,
                    "U":99.5/100
                        }, mutable=True,
                        doc="The fraction of component that goes from aqueous to organic phase")
                
                self.add_component(
                        stream + "_partition_coefficient",
                        partition_coefficient
                    )

                
                def distribution_extent_rule(b, t, s, j):
                    if j in ppack.dissolved_elements:
                      if s == self.elements.first():
                        return distribution_extent[t, s, j] == in_state[t].get_material_flow_terms(j)*partition_coefficient[j]
                      else:
                        return distribution_extent[t, s, j] == state_block[t, s-1].get_material_flow_terms(j)*partition_coefficient[j]
                    return Constraint.Skip
                
                distribution_extent_constraint = Constraint(self.flowsheet().time, self.elements,
                                                            ppack.dissolved_elements, rule=distribution_extent_rule)

                self.add_component(
                    stream + "_distribution_extent_constraint",
                    distribution_extent_constraint,
                )

                def material_balance_aq_rule(b, t, s, j):
                    for aqstream, pconfig in b.config.aqueous_streams.items():
                            
                            in_state_a, out_state_a, side_state_a = _get_state_blocks(b, t, s, aqstream)


                            if in_state_a is not None:
                                rhsa =  in_state_a.get_material_flow_terms(j) - out_state_a.get_material_flow_terms(j)
                                
                            else:
                                rhsa = -out_state_a.get_material_flow_terms(j)
                        
                            # Aq streams always have a distribution extent
                            if j != 'H2SO4':
                                rhsa += -distribution_extent[t, s, j]
                            
                            #if j != 'H2SO4':
                                #return rhsa == state_block[t, s].aqueous_vol*dC_dt_a[t, s, j]
                            #else:
                            return rhsa == 0

                      
                
                mbal = Constraint(self.flowsheet().time, self.elements, component_list, rule=material_balance_aq_rule)
                self.add_component(stream + "_material_balance", mbal)

            for stream, sconfig in self.config.organic_streams.items():
                state_block = getattr(self, stream)
                ppack = sconfig.property_package
                component_list = state_block.component_list

                # Dynamic part needs to be repaired

                #dC_o = Var( self.flowsheet().time,
                    #self.elements,
                    #ppack.dissolved_elements,
                    #domain=Reals,
                        #initialize=0.0)
                
                #self.add_component(
                        #stream + "dC_o",
                        #dC_o
                    #)
                
                #def dynamic_term_o(b,t,s,j):
                    #for j in ppack.dissolved_elements:
                        #if self.config.dynamic == True:
                            #return dC_o[t,s,j] == state_block[t, s].conc_mass_comp[j]
                        #else:
                            #return dC_o[t,s,j] == 0
                    #return Constraint.Skip
                
                #dynamic_term_o_constraint = Constraint(self.flowsheet().time, self.elements,
                                                            #ppack.dissolved_elements, rule=dynamic_term_o)
                
                #self.add_component(
                    #stream + "_dynamic_term_o_constraint",
                    #dynamic_term_o_constraint,
                #)

                #dC_dt_o = DerivativeVar(dC_o, wrt=self.flowsheet().time)

                #self.add_component(
                        #stream + "_dC_dt_o",
                        #dC_dt_o
                    #)


                def material_balance_og_rule(b, t, s, j):
                    for ogstream, pconfig in b.config.organic_streams.items():
                            in_state_o, out_state_o, side_state_o = _get_state_blocks(b, t, s, ogstream)

                            if in_state_o is not None:
                                rhso =  in_state_o.get_material_flow_terms(j) - out_state_o.get_material_flow_terms(j)
                                
                            else:
                                rhso = -out_state_o.get_material_flow_terms(j)
                            
                            if j not in ['H2SO4','DEHPA']:
                                rhso += distribution_extent[t, s, j]
                            
                            # Dynamic part to be repaired

                            #if j not in ['H2SO4','DEHPA']:
                                #return rhso == state_block[t, s].organic_vol*dC_dt_o[t, s, j] 
                            #else:
                            return rhso == 0

                                    
                mbal_og = Constraint(self.flowsheet().time, self.elements, component_list, rule=material_balance_og_rule)
                
                self.add_component(stream + "_material_balance", mbal_og)

    # Creation of ports          
                    
    def _build_ports(self):
            # Add Ports
            for aqstream, pconfig in self.config.aqueous_streams.items():
                sblock = getattr(self, aqstream)
                flow_dir = pconfig.flow_direction

                if pconfig.has_feed:
                    inlet_state = getattr(self, aqstream + "_inlet_state")
                    in_port, _ = inlet_state.build_port(
                        f"{aqstream} Inlet", slice_index=(slice(None))
                    )
                    self.add_component(aqstream + "_inlet", in_port)

                if flow_dir == FlowDirection.forward:
                    outlet = self.elements.last()
                elif flow_dir == FlowDirection.backward:
                    outlet = self.elements.first()
                else:
                    raise BurntToast("If/else overrun when constructing Ports")

                out_port, _ = sblock.build_port(
                    f"{aqstream} Outlet", slice_index=(slice(None), outlet)
                )
                self.add_component(aqstream + "_outlet", out_port)
            
            for ogstream, pconfig in self.config.organic_streams.items():
                sblock = getattr(self, ogstream)
                flow_dir = pconfig.flow_direction

                if pconfig.has_feed:
                    inlet_state = getattr(self, ogstream + "_inlet_state")
                    in_port, _ = inlet_state.build_port(
                        f"{ogstream} Inlet", slice_index=(slice(None))
                    )
                    self.add_component(ogstream + "_inlet", in_port)

                if flow_dir == FlowDirection.forward:
                    outlet = self.elements.last()
                elif flow_dir == FlowDirection.backward:
                    outlet = self.elements.first()
                else:
                    raise BurntToast("If/else overrun when constructing Ports")

                out_port, _ = sblock.build_port(
                    f"{ogstream} Outlet", slice_index=(slice(None), outlet)
                )
                self.add_component(ogstream + "_outlet", out_port)

# Creation of inlet and outlet state blocks for the streams

def _get_state_blocks(b, t, s, stream):
    """
    Utility method for collecting states representing flows into and out of
    a stage for a given stream.
    """
    if stream in b.config.aqueous_streams.keys():
        streams = b.config.aqueous_streams
    elif stream in b.config.organic_streams.keys():
        streams = b.config.organic_streams
    else:
        raise BurntToast("If/else overrun when constructing balances")

    state_block = getattr(b, stream)
    if streams[stream].flow_direction == FlowDirection.forward:
        if s == b.elements.first():
            if not streams[stream].has_feed:
                in_state = None
            else:
                in_state = getattr(b, stream + "_inlet_state")[t]
        else:
            in_state = state_block[t, b.elements.prev(s)]
    elif streams[stream].flow_direction == FlowDirection.backward:
        if s == b.elements.last():
            if not streams[stream].has_feed:
                in_state = None
            else:
                in_state = getattr(b, stream + "_inlet_state")[t]
        else:
            in_state = state_block[t, b.elements.next(s)]
    else:
        raise BurntToast("If/else overrun when constructing balances")

    out_state = state_block[t, s]

    # Look for side state
    side_state = None
    if hasattr(b, stream + "_side_stream_state"):
        try:
            side_state = getattr(b, stream + "_side_stream_state")[t, s]
        except KeyError:
            pass

    return in_state, out_state, side_state
