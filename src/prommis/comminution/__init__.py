#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Particle-size reduction and classification models for PrOMMiS.

The package provides particle size distribution (PSD) calculations, a
solid PSD property package, comminution and classification functions,
and IDAES unit models for particle-size reduction and classification operations.
The package-level imports are experimental and may change before release.
"""

from prommis.comminution.properties.solid_psd_properties import (
    SolidPSDParameterBlock,
)
from prommis.comminution.unit_models.crusher import CrusherSolidPSD

__all__ = ["CrusherSolidPSD", "SolidPSDParameterBlock"]
