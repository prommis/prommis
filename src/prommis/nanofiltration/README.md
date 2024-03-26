# membrane examples
The following files serve as a demonstration of membrane unit models (sourced from WaterTAP) implemented for lithium brine purification.

## nanofiltration
### nf_brine.py
This is an example flowsheet for the separation of lithium and magnesium ions using nanofiltration (nf). The unit model is the DSPM-DE and the property model is the MCAS, both of which can be found in WaterTAP. Currently, the inlet feed is modeling the Salar de Atacama brine. 

The objective of the pyomo model is to limit the amount of lithum lost in the retentate stream, with 2 DOF (pressure, area). Operating at this extreme end of the tradeoff bewteen the selectivity and permeability, combined with the low solvent flux for this example, the magnesium rejection is low (~14%). An additional constraint is enacted on the model to limit the feed pressure to 70 bar, per typical experimental capabilities.

good references: https://doi.org/10.1021/acs.est.2c08584 ; https://doi.org/10.1016/j.memsci.2020.118809

### test_nf_brine.py
This file tests the above nf flowsheet.

### nf_brine_plot.py
This file uses the nf_brine.py flowsheet to perform a sensitivity analysis on the volume recovery of the nanofiltration system. The file outputs four subplots:

1. Ion Rejection versus Volume Recovery: As the volume recovery (i.e., the relative amount of total solution recovered in the permeate) increases, the lithium rejeciton monotonically decreases and the magnesium rejection shows a local maxima. This is the result of the tradeoff between the permeability and selectivity of the membrane.

2. Mg:Li Ratio versus Volume Recovery: The ratio of Mg:Li in the feed is constant, which is represented by the red line. All results tested show that the permeate stream has smaller Mg:Li ratios than the feed, which is desirable as the end goal is purifying lithium. There is a local minima at a volume recovery of ~55%.

3. Membrane Area (m2) versus Volume Recovery: The membrane area increases as the system capacity (volume recovery) increases.

4. Feed Pressure (bar) versus Volume Recovery: The optimal feed pressure increases as the system capacity increases and the membrane area reaches its upper bound.
