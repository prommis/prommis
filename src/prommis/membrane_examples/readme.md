##### membrane examples #####
The following files serve as a demonstration of membrane unit models (sourced from WaterTAP) implemented for lithium brine purification.

Based on the PR checks, there may be some dependencies to sort out (e.g., watertap).

### nanofiltration
# nf_brine.py
This is an example flowsheet for the separation of lithium and magnesium ions using nanofiltration. The unit model is DSPM-DE and the property model is MCAS, both of which can be found in WaterTAP. Currently, the inlet feed is modeling the Salar de Atacama brine. The magnesium rejection is low, which is likely due to the low water flux.
good references: https://doi.org/10.1021/acs.est.2c08584 ; https://doi.org/10.1016/j.memsci.2020.118809

The optimization objective is maximizing the reduction of Mg:Li in the permeate compared to the feed, with 2 DOF (pressure, area). The model struggles with additional bounds.

# test_nf_brine.py
This file tests the above nf flowsheet.


### reverse osmosis
# LiCl_prop_pack.py
This is the property package for the solution diffusion model needed for reverse osmosis. This package was based on the WaterTAP package written for NaCl for use in the flowsheet below.

# ro_brine.py
This is an example flowsheet for the concentration of lithium using reverse osmosis. The RO unit model (which can be found in WaterTAP) can only handle total dissolved solids (TDS), rather than individual ions. As such, this example focuses solely on lithium, using chlorine in equimolar quantities as the counter ion. The inlet flows are based on the below reference, which achieves a lithium ion recovery of 2.3%.
reference: https://doi.org/10.1016/j.desal.2020.114620

The optimization objective is minimizing the specific energy consumption with 1 DOF (area). The model struggles with additional bounds and with dilute concentrations.

There are still some unit inconsistencies (see output from model diagnostics) regarding the osmotic pressure which need to be addressed in a manner that does not cause failed initialization of the unit.

# test_ro_brine.py
This file tests the above ro flowsheet.