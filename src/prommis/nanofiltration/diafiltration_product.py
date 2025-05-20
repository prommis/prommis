#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""
Custom product block used with the two-salt diafiltration unit model.

Author: Molly Dougher
"""

from pyomo.environ import Reference
from pyomo.network import Port

from idaes.core import declare_process_block_class
from idaes.models.unit_models.product import ProductData


@declare_process_block_class("DiafiltrationProduct")
class DiafiltrationProductData(ProductData):
    """
    Modification of the Product Unit Model for use in the Two-Salt
    Diafiltration Unit Model.
    """

    CONFIG = ProductData.CONFIG()

    def build(self):
        super().build()
        self.add_inlet_port()

    def add_inlet_port(self):
        """
        Updates the inlet Port to the product block to incorporate the time
        index within the Reference.
        """
        # remove defined self.inlet Port and References from super().build()
        self.del_component(self.inlet)
        self.del_component(self._flow_vol_inlet_ref)
        self.del_component(self._conc_mass_lithium_inlet_ref)
        self.del_component(self._conc_mass_cobalt_inlet_ref)
        self.del_component(self._conc_mass_chlorine_inlet_ref)

        # create self.inlet Port with adjusted References
        self.inlet = Port(doc="Inlet Port")
        self._flow_volume_ref = Reference(self.flow_vol[0])
        self.inlet.add(self._flow_volume_ref, "flow_vol")
        self._conc_mass_lithium_ref = Reference(self.conc_mass_lithium[0])
        self.inlet.add(self._conc_mass_lithium_ref, "conc_mass_lithium")
        self._conc_mass_cobalt_ref = Reference(self.conc_mass_cobalt[0])
        self.inlet.add(self._conc_mass_cobalt_ref, "conc_mass_cobalt")
        self._conc_mass_chlorine_ref = Reference(self.conc_mass_chlorine[0])
        self.inlet.add(self._conc_mass_chlorine_ref, "conc_mass_chlorine")
