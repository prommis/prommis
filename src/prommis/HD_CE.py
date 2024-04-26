from pyomo.environ import ConcreteModel, Constraint, Param, Var, units,PositiveReals, log
from pyomo.environ import units as pyunits


from idaes.core import (
    FlowsheetBlock,
    UnitModelCostingBlock,
    UnitModelBlock,
)
from idaes.models.costing.SSLW import (
    SSLWCosting,
    SSLWCostingData,
    VesselMaterial,
)
from idaes.core.solvers import get_solver

from math import pi 



class HDFurnaceCostEstimator:

    """Cost Estimation of an Hydrogen Decrepitation Furnace"""

    def __init__(self,
                ramp_up_time = 20,                                                  #Ramp Up time (in minutes)   
                sample_shape = "cylindrical",
                operating_temperature = 443.15,                                       #in Kelvin
                decrepitation_duration = 3,                                     #Amount of time the furnace runs at operating temperature(in hr)
                preparation_time = 1,                                           # in hr
                cool_down_time = 1,                                             # in hr
                sample_cp = 0.44,                                               #in kJ/(kg*K)
                diameter = 0.152,
                length=0.456,
                breadth=0.1524, 
                width=0.2191,
                sample_mass =  62.025,                                          #kg 
                length_insulation1 = 7.62,                                      #m
                width_insulation1 = 0.6096,                                     #m
                thickness_insulation1 = 0.0254,                                 #m
                price_insulation1 = 65.16,                                      #USD
                weight_insulation1 = 15.12,                                     #kg
                length_insulation2 = 1.2192,                                    #m
                width_insulation2 = 0.6096,                                     #m
                thickness_insulation2 = 0.0508,                                 #m
                price_insulation2 = 35.76,                                      #USD
                weight_insulation2 = 1.814,                                     #kg                                                                  
                metal_material1 = VesselMaterial.StainlessSteel304,
                metal_material2 = VesselMaterial.CarbonSteel,
                hours_per_shift=8,
                shifts_per_day=3,
                specific_heat_capacity_insulation1 = 1.08,                        #KJ/(kg*K)                                       
                specific_heat_capacity_Metal1 = 0.468,                            #KJ/(kg*K)
                specific_heat_capacity_insulation2 =0.9,                          #KJ/(kg*K)
                specific_heat_capacity_Metal2 = 0.502416,                         #KJ/(kg*K)
                operating_days_per_year=336,
                efficiency = 0.95,
                electricity_rate = 0.0683,
                labor_rate = 75,                                                  #USD
                transformer_cost = 2647.50, 
                temperature_controller_price = 129.00
                ):
        
        """length_insulation1&2, width_insulation1&2, thickness_insulation1&2 ,                                                
        price_insulation1&2, weight_insulation1&2 are vendor's specification                                                                                                                                                         

        Sample shapes: cylindrical - 1, rectangular cuboid - 2, square faced cuboid - 3, 
        spherical - 4
        """

        self.sample_shape = sample_shape                                    
        self.operating_temperature = operating_temperature
        self.decrepitation_duration = decrepitation_duration
        self.sample_cp = sample_cp
        self.diameter = diameter
        self.length = length
        self.breadth = breadth
        self.width = width
        self.sample_mass = sample_mass
        self.length_insulation1 = length_insulation1
        self.width_insulation1 = width_insulation1
        self.thickness_insulation1 = thickness_insulation1
        self.price_insulation1 = price_insulation1
        self.weight_insulation1 = weight_insulation1
        self.length_insulation2 = length_insulation2
        self.width_insulation2 = width_insulation2
        self.thickness_insulation2 = thickness_insulation2
        self.price_insulation2 = price_insulation2
        self.weight_insulation2 = weight_insulation2
        self.metal_material1 = metal_material1
        self.metal_material2 = metal_material2
        self.ramp_up_time = ramp_up_time
        self.hours_per_shift = hours_per_shift
        self.shifts_per_day = shifts_per_day
        self.specific_heat_capacity_insulation1 = specific_heat_capacity_insulation1
        self.specific_heat_capacity_Metal1 = specific_heat_capacity_Metal1
        self.specific_heat_capacity_insulation2 = specific_heat_capacity_insulation2
        self.specific_heat_capacity_Metal2 = specific_heat_capacity_Metal2
        self.operating_days_per_year = operating_days_per_year
        self.efficiency = efficiency
        self.electricity_rate = electricity_rate
        self.labor_rate = labor_rate
        self.transformer_cost = transformer_cost
        self.temperature_controller_price = temperature_controller_price
        self.preparation_time = preparation_time
        self.cool_down_time = cool_down_time


       
    #Internal volume of furnace
        
    def furnace_internal_volume(self):                                            
        #cylindrical - 1, rectangular cuboid - 2, square faced cuboid - 3, spherical - 4
        
        if self.sample_shape == 1:
            if self.length == None:
                l = 3 * self.diameter
            else:
                l = self.length
            volume = ((pi * (self.diameter)**2)/4)*l

        elif self.sample_shape == 2:
            volume = self.length*self.breadth*self.width
        elif self.sample_shape == 3:
            volume = self.length**3
        elif self.sample_shape == 4:
            volume = (pi*(self.diameter**3))/6
        else: #didn't match any of the expected strings    
            raise TypeError(f"Shape type {self.sample_shape} is not a valid type.")
        
        internal_volume = 2 * volume
       
        return  internal_volume
    
    
    
    #Heat loss and optimum thickness of each materials

    def heat_loss(self):
        m = ConcreteModel()

        m.temp1 = Var(within=PositiveReals, units = units.K, bounds = (950.13, 1173.15) )
        m.T1 = Var(within=PositiveReals, units = units.m,bounds = (1e-6,None))
        m.T2 = Var(within=PositiveReals, units = units.m,bounds = (1e-6,None))
        m.T3 = Var(within=PositiveReals, units = units.m,bounds = (1e-6,None))
        m.T4 = Var(within=PositiveReals, units = units.m,bounds = (1e-6,None))
        m.q = Var(within=PositiveReals, units = units.W)



        m.r = Param(initialize=(self.furnace_internal_volume()/(6*pi))**(1/3), units=units.m, mutable=True)
        m.K1 = Param(initialize=0.33,units = units.W/units.m/units.K, mutable=True)
        m.K2 = Param(initialize=13.53, units = units.W/units.m/units.K, mutable=True)
        m.K3 = Param(initialize=0.069,units = units.W/units.m/units.K, mutable=True)
        m.K4 = Param(initialize=45,units = units.W/units.m/units.K, mutable=True)
        m.l = Param(initialize=3*(2 * m.r), units=units.m, mutable=True)                     


        m.c1 = Constraint(expr=m.q == 0.70 * 2 * pi * (m.r-(1e-3)) *m.l*(5.67*1e-8)*((1173.15**4) - (m.temp1**4)))                          
        m.c2 = Constraint(expr=m.q == (2 * pi * m.K1 *m.l* (m.temp1 - 950.13)) / (log((m.r + m.T1) / m.r)))
        m.c3 = Constraint(expr=m.q == ((2 * pi * m.K2*m.l) * (950.13-950.00)) / (log((m.r + m.T1 + m.T2) / (m.r + m.T1))))
        m.c4 = Constraint(expr=m.q == ((2 * pi * m.K3*m.l) * (950.00-333.18)) / (log((m.r + m.T1 + m.T2 + m.T3) / (m.r + m.T1 + m.T2))))
        m.c5 = Constraint(expr=m.q == ((2 * pi * m.K4*m.l) * (333.18-333.15)) / (log((m.r + m.T1 + m.T2 + m.T3 + m.T4) / (m.r + m.T1 + m.T2 + m.T3))))
        m.c6 = Constraint(expr=m.q == 5 * 2 * pi *m.l* (333.15 - 298.15) * (m.r + m.T1 + m.T2 + m.T3 + m.T4))

        Solver = get_solver('ipopt')
        result = Solver.solve(m, tee= False)

        v = [m.temp1.value, m.T1.value, m.T2.value, m.T3.value, m.T4.value, m.q.value]
        return v
    
    #Cost of the high temperature insulation
    def C_Ins1(self):                                               
        T1 = self.heat_loss()[1]
        r = (self.furnace_internal_volume()/(6 * pi))**(1/3)
        d_ext = 2 * (r+T1)
        l = 3 * (2 * r)
        T1 = self.heat_loss()[1]
        V = pi * l * ((r+T1)**2) - (r**2)
        N = V/(self.length_insulation1*self.width_insulation1*self.thickness_insulation1)                                                                 
        N = max(1,N)
        wt = N*self.weight_insulation1
        ext_surface_area = (pi*d_ext*l)+((pi*(d_ext**2))/2)
        material_cost = N*self.price_insulation1
        labor_hours_rqd = max(1,1.956*ext_surface_area)
        labor_cost = (4+labor_hours_rqd)*self.labor_rate
        total_cost = material_cost + labor_cost
        p = [wt, total_cost, material_cost, ext_surface_area]
        return p

    #Cost of metal material 1 

    def C_Metal1(self):
        T1, T2 = self.heat_loss()[1:3]
        r = (self.furnace_internal_volume()/(6*pi))**(1/3)
        d = 2*(r+T1)
        d_ext = 2*(r+T1+T2)
        l = 3*(r*2)
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.costing = SSLWCosting()
        m.fs.unit = UnitModelBlock()


        m.fs.unit.costing = UnitModelCostingBlock(
                flowsheet_costing_block=m.fs.costing,
                costing_method=SSLWCostingData.cost_vessel,
                costing_method_arguments={
                    "vertical": False,
                    "material_type": self.metal_material1,
                    "include_platforms_ladders": False,
                    "shell_thickness": T2*pyunits.m,
                    "vessel_diameter": d * pyunits.m,
                    "vessel_length": l * pyunits.m,                        
                    "weight_limit": 1
                },
            )
        solver = get_solver()
        results = solver.solve(m, tee=True)
        M1_material_cost = m.fs.unit.costing.capital_cost.value
        M1_weight = m.fs.unit.costing.weight.value
        ext_surface_area = (pi*d_ext*l)+((pi*(d_ext**2))/2)
        int_surface_area = (pi*d*l)+((pi*(d**2))/2)
        labor_hours_assembling_and_welding = (6.111*1e-3*M1_weight) + 174.89
        labor_hours_sandblasting = max(1,(ext_surface_area + int_surface_area)*0.474)
        labor_cost = (labor_hours_assembling_and_welding + labor_hours_sandblasting)*self.labor_rate
        total_M1_cost =  M1_material_cost + labor_cost
        M = [M1_weight, total_M1_cost, M1_material_cost, ext_surface_area]
        
        return M
    
    #Cost of the low temperature insulation

    def C_Ins2(self): 
        T1, T2, T3 = self.heat_loss()[1:4]
        r = (self.furnace_internal_volume()/(6*pi))**(1/3)
        d_ext = 2*(r+T1+T2+T3)
        l = 3*(2*r)
        V = pi*l*((r+T3+T2+T1)**2)-((r+T2+T1)**2)
        N = V/(self.length_insulation2*self.width_insulation2*self.thickness_insulation2) 
        N = max(1, N)
        wt = N*self.weight_insulation2
        material_cost = N*self.price_insulation2
        ext_surface_area = (pi*d_ext*l)+((pi*(d_ext**2))/2)
        labor_hours_rqd = 1.956*ext_surface_area
        labor_hours_rqd = max(1, 1.956 * ext_surface_area)
        labor_cost = (4+labor_hours_rqd)*self.labor_rate
        total_cost = material_cost + labor_cost
        p = [wt, total_cost, material_cost, ext_surface_area]
        return p

    def C_Metal2(self):
        T1, T2, T3, T4 = self.heat_loss()[1:5]
        r = (self.furnace_internal_volume()/(6*pi))**(1/3)
        d = 2*(r+T1+T2+T3)
        d_ext = 2*(r+T1+T2+T3+T4)
        l = 3*(2*r)                                                                 
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.costing = SSLWCosting()
        m.fs.unit = UnitModelBlock()

        m.fs.unit.costing = UnitModelCostingBlock(
                flowsheet_costing_block=m.fs.costing,
                costing_method=SSLWCostingData.cost_vessel,
                costing_method_arguments={
                    "vertical": False,
                    "material_type": self.metal_material2,
                    "include_platforms_ladders": False,
                    "shell_thickness": T4*pyunits.m,
                    "vessel_diameter": d * pyunits.m,
                    "vessel_length": l * pyunits.m,                            
                    "weight_limit": 1
                },
            )
        solver = get_solver()
        results = solver.solve(m, tee=True)
        M2_material_cost = m.fs.unit.costing.capital_cost.value
        M2_weight = m.fs.unit.costing.weight.value
        M2_volume = (pi*(d_ext**2)*l)/4
        ext_surface_area = (pi*d_ext*l)+((pi*(d_ext**2))/2)
        int_surface_area = (pi*d*l)+((pi*(d**2))/2)
        labor_hours_assembling_and_welding = (6.111*1e-3*M2_weight) + 174.89
        labor_hours_sandblasting_and_painting = 30+(2.77*ext_surface_area)+(int_surface_area*0.474)
        labor_cost = (labor_hours_assembling_and_welding + labor_hours_sandblasting_and_painting)*self.labor_rate
        total_M2_cost =  M2_material_cost + labor_cost
        M = [M2_weight, total_M2_cost, M2_material_cost, ext_surface_area, M2_volume]
        return M
    
    def energy_req(self):                                                                                     
        #	Energy required to raise the interior temperature of the furnace to operating_temperature
        Q1 = ((self.C_Ins1()[0])*self.specific_heat_capacity_insulation1*((490.13 - 298.15)))+((((self.C_Metal1()[0])*0.453592))*self.specific_heat_capacity_Metal1*(490.00-298.15))+((self.C_Ins2()[0])*self.specific_heat_capacity_insulation2*(333.18-298.15))+(((self.C_Metal2()[0])*0.453592)*self.specific_heat_capacity_Metal2*(333.15-298.15))
        Q1_spec = Q1/self.sample_mass
        Q1_R = (Q1/(self.ramp_up_time*60))                                                   #kW
        
        Q2 = self.heat_loss()[5]                                                        #W                                                            #Heat loss from the external surface area of a furnace
        p2 = Q2/self.efficiency                                                         #W                                                                            #power
        E = p2*self.decrepitation_duration                                              #W-hr                                                               #W-hr
        Q2_spec = (E*3.6)/self.sample_mass                                              #kJ/kg of Sample                                                                                                         #kJ/kg

        #Energy (specific) required to raise the temperature of the sample from room temperature to operating temperature
        Q3 = self.sample_mass*self.sample_cp*(self.operating_temperature - 298.15)
        Q3_spec = self.sample_cp*(self.operating_temperature - 298.15)
        Q3_R =  (Q3/(self.ramp_up_time*60))                                                         #KJ/KG of sample

        Actual_coil_rating = (Q1_R+Q3_R)/self.efficiency

        M = [Q1_spec, Q2_spec, Q3_spec, Actual_coil_rating]
        return  M
    
    def HDF_power(self):
        
        Q2 = self.heat_loss()[5]                                                                                #Heat loss from the external surface area of a furnace
        power = Q2/self.efficiency
        return power
    
    def HDF_OPEX(self):
        M = self.energy_req()

        Q_total = M[0] + M[1] + M[2]                                                                                                                        #KJ/KG of sample
        total_runs = (self.hours_per_shift*self.shifts_per_day*self.operating_days_per_year)/self.decrepitation_duration                                 #Total number of runs in a year 
        total_energy_required = (Q_total * total_runs)/3600
        utility = (total_energy_required * self.electricity_rate)                                                                                        #USD Per kg of sample
        processing_time = (self.preparation_time + (self.ramp_up_time/60) + self.decrepitation_duration+self.cool_down_time)/self.sample_mass                 # hr/kg
        opex = (1/processing_time) * (utility) * (8000)                                            #USD per year

        return opex                                                              
        
    def heating_coil(self):
        M = self.energy_req()[3]
        if 0<M <=5:
            price = 80
        elif M>5:                                                                
            price = (14.81*M) + 5.93
        return price
    
    def HDF_CAPEX(self):

        """
        transformer_cost specification - Transformer, dry type, DOE 2016, 30kVA, 3 phase, 480V delta primary, 208Y/120V secondary, 18M, 150C rise
        temperature_pressure_controller_price - N480D PID Temperature controller
        """

        vessel_weight = self.C_Ins1()[0]+self.C_Metal1()[0]+self.C_Ins2()[0]+self.C_Metal2()[0]
        vessel_volume = self.C_Metal2()[4]*35.315
        vessel_density = vessel_weight/vessel_volume
        transportation_cost = (10+(0.57*vessel_density))*self.labor_rate 
        inspection_cleaning_and_testing = (50+(1.892*self.C_Metal2()[3]))*self.labor_rate

        CAPEX = self.C_Ins1()[1] + self.C_Metal1()[1] + self.C_Ins2()[1]+ self.C_Metal2()[1] + self.transformer_cost + self.temperature_controller_price+self.heating_coil() + transportation_cost + inspection_cleaning_and_testing                                             
        return CAPEX
                                                                                               #USD per kg of sample

HD1 = HDFurnaceCostEstimator(sample_shape=1)
print(HD1.HDF_OPEX()) 
  

