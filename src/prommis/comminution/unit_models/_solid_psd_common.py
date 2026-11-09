#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Helpers for solid PSD unit models.

They declare property-package options, check the selected package, report
applicability warnings, and register or rebuild precomputed Params.
"""

__author__ = "Daison Yancy Caballero"


from idaes.core.util.exceptions import ConfigurationError

from prommis.comminution.properties.solid_psd_properties import SolidPSDParameterData


def require_solid_psd_package(config, unit_name):
    """Return the configured solid PSD package or raise ``ConfigurationError``."""
    pp = config.property_package
    if not isinstance(pp, SolidPSDParameterData):
        raise ConfigurationError(
            f"{unit_name} requires a SolidPSDParameterBlock property package."
        )
    return pp


def report_applicability_warnings(logger, unit_name, messages, emit_warning=True):
    """Return applicability messages and optionally log them to ``logger``."""
    messages = list(messages)
    if emit_warning:
        for msg in messages:
            logger.warning("%s applicability: %s", unit_name, msg)
    return messages
