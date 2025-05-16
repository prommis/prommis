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

from pyomo.common.config import ConfigValue
from pyomo.dae import ContinuousSet
from pyomo.environ import TransformationFactory

from idaes.core import declare_process_block_class
from idaes.models.unit_models.product import ProductData


@declare_process_block_class("DiafiltrationProduct")
class DiafiltrationProductData(ProductData):
    """
    Modification of the Product Unit Model for use in the Two-Salt
    Diafiltration Unit Model.

    Discretizes the property package over the width of the membrane.
    """

    CONFIG = ProductData.CONFIG()

    CONFIG.declare(
        "NFEx",
        ConfigValue(
            doc="Number of discretization points in the x-direction",
        ),
    )

    def build(self):
        super().build()
        self.add_length()
        self.add_state_blocks()
        self.discretize_length()
        self.add_inlet_port()

    def add_length(self):
        """
        Creates an equivalent length scale as the memrbane width in the
        diafiltration membrane unit model with the product block.
        """
        self.x_bar = ContinuousSet(bounds=(0, 1))

    def add_state_blocks(self):
        """
        Adds to state blocks for the property package that will be indexed over
        time and the length scale (membrane width).
        """
        # remove defined self.properties from super().build()
        self.del_component(self.properties)

        # define self.properties such that it is discretized over t, x
        self.properties = self.config.property_package.build_state_block(
            self.flowsheet().time,
            self.x_bar,
            doc="Material properties",
            initialize=dict(**self.config.property_package_args),
        )

    def discretize_length(self):
        """
        Discretizes the product block over the same number of finite elements
        as the diafiltration membrane such that the ports have properly
        matched variables.
        """
        discretizer = TransformationFactory("dae.finite_difference")
        discretizer.apply_to(
            self, wrt=self.x_bar, nfe=self.config.NFEx, scheme="BACKWARD"
        )

    def add_inlet_port(self):
        """
        Adds the inlet Port to the product block that will include the
        References to the discretizes state block properties,
        """
        # remove defined self.inlet Port and References from super().build()
        self.del_component(self.inlet)
        self.del_component(self._flow_vol_inlet_ref)
        self.del_component(self._conc_mass_comp_inlet_ref)

        # create self.inlet Port such that it is discretized over t, x
        self.add_port(name="inlet", block=self.properties, doc="Inlet Port")
