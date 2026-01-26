#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Two-Salt Diafiltration Unit Model
=================================

Author: Molly Dougher

This membrane unit model is for the multi-component diafiltration of a two-salt system with a common anion. The membrane is designed for use in a diafiltration cascade, i.e., the model represents one spiral-wound membrane module.

Configuration Arguments
-----------------------

The Two-Salt Diafiltration unit model requires a property package that provides the valency (:math:`z_i`), single solute diffusion coefficient (:math:`D_i`) in :math:`\mathrm{mm}^2 \, \mathrm{h}^{-1}`, reflection coefficient (:math:`\sigma_i`), partition coefficient (:math:`H_i`) at the solution-membrane interface(s), and number of dissolved species (:math:`n_i`) for each ion :math:`i` in solution. When used in a flowsheet, the user can provide separate property packages for the feed and product streams.

There are two required arguments, ``NFE_module_length`` and ``NFE_membrane_thickness``, to specify the desired number of finite elements across the width (module length) and thickness of the membrane, respectively.

Degrees of Freedom
------------------

The Two-Salt Diafiltration unit model has three degrees of freedom (variable names and default values are provided in parentheses):

#. the length of the membrane module (``total_module_length``; :math:`4 \, \mathrm{m}`)
#. the length of the membrane (``total_membrane_length``; :math:`41 \, \mathrm{m}`)
#. the pressure applied to the membrane system (``applied_pressure``; :math:`10 \, \mathrm{bar}`)

To run a simulation (with zero degrees of freedom) in a flowsheet, the following variables must be fixed to obtain zero degrees of freedom (variable names and default values are provided in parentheses):

#. the volumetric flow rate of the feed (``feed_flow_volume``; :math:`12.5 \, \mathrm{m}^3 \, \mathrm{h}^{-1}`)
#. the lithium concentration in the feed (``feed_conc_mol_comp[t,"Li"]``; :math:`245 \, \mathrm{mol} \, \mathrm{m}^{-3}`)
#. the cobalt concentration in the feed (``feed_conc_mol_comp[t,"Co"]``; :math:`288 \, \mathrm{mol} \, \mathrm{m}^{-3}`)
#. the volumetric flow rate of the diafiltrate (``diafiltrate_flow_volume``; :math:`3.75 \, \mathrm{m}^3 \, \mathrm{h}^{-1}`)
#. the lithium concentration in the diafiltrate (``diafiltrate_conc_mol_comp[t,"Li"]``; :math:`14 \, \mathrm{mol} \, \mathrm{m}^{-3}`)
#. the cobalt concentration in the diafiltrate (``diafiltrate_conc_mol_comp[t,"Co"]``; :math:`3 \, \mathrm{mol} \, \mathrm{m}^{-3}`)

Model Structure
---------------

There are three phases in the Two-Salt Diafiltration model: the retentate, the membrane, and the permeate. The retentate and the permeate are only discretized with respect to :math:`x` (parallel to the membrane surface), while the membrane is discretized with respect to both :math:`x` and :math:`z` (perpendicular to the membrane surface). The resulting system of partial differential algebraic equations is solved by discretizing with the backward finite difference method.

Assumptions
-----------

The membrane module dimensions, maximum applied pressure, and inlet flow rates assume that one tube (one instance of this model) consists of 4 NF270-440 membranes in series.

The partitioning relationships, which describe how the solutes transition (partition) across the solution-membrane interfaces, are derived assuming Donnan equilibrium. The partitioning coefficients incorporate both steric and Donnan effects.

The default value for the membrane's surface charge (:math:`-140 \, \mathrm{mM}`), was calculated using zeta potential measurements for NF270 membranes. (See `this reference <https://doi.org/10.1021/acs.iecr.4c04763>`_). Currently, the default property package only supports negatively charged membranes.

The membrane is assumed to be :math:`100 \, \mathrm{nm}` thick.

The default value for the membrane permeability (:math:`0.01 \, \mathrm{m} \, \mathrm{h}^{-1} \, \mathrm{bar}^{-1}`) is based off of parameter estimation results from `this reference <https://doi.org/10.1021/acs.iecr.4c04763>`_ for NF270 membranes.

The formation of a boundary layer at the membrane surface due to concentration polarization is neglected for mathematical simplicity.

The dominating transport mechanism within the bulk/retentate solution is convection in the :math:`x`-direction (parallel to the membrane surface). The dominating transport mechanism within the permeate solution is convection in the :math:`z`-direction (perpendicular to the membrane surface).

The transport mechanisms modeled within the membrane are convection, diffusion, and electromigration. Diffusion within the membrane that is normal to the pore walls is ignored, meaning the concentration gradient of ion :math:`i` within the membrane only has a :math:`z`-component (perpendicular to the membrane surface).

Sets
----

The Two-Salt Diafiltration model defines the following discrete sets for ions and cations in the system, respectively:

.. math:: \mathcal{I}=\{\mathrm{Li,Co,Cl}\}
.. math:: \mathcal{K}=\{\mathrm{Li,Co}\}

There are 2 continuous sets for each length dimension: ``dimensionless_module_length`` (in the :math:`x`-direction parallel to the membrane surface) and ``dimensionless_membrane_thickness`` (in the :math:`z`-direction perpendicular to the membrane surface). :math:`x` and :math:`z` are non-dimensionalized (denoted as :math:`\bar{x}` and :math:`\bar{z}`, respectively) using the module length or (:math:`w`) and membrane thickness (:math:`l`), respectively, to improve numerics.

.. math:: \bar{x} \in \mathbb{R} \| 0 \leq \bar{x} \leq 1
.. math:: \bar{z} \in \mathbb{R} \| 0 \leq \bar{z} \leq 1

Some variables have a time domain to be compatible with the property package, even though this is not a dynamic model. Thus, the following set is defined for time.

.. math:: t \in [0]

Default Model Parameters
------------------------

The Two-Salt Diafiltration model has the following parameters.

================ =============================================== ============================ ============= ==========================================================
Parameter        Description                                     Name                         Default Value Units
================ =============================================== ============================ ============= ==========================================================
:math:`\epsilon` numerical tolerance for zero values             ``numerical_zero_tolerance`` 1e-10
:math:`l`        thickness of the membrane                       ``total_membrane_thickness`` 1e-07         :math:`\mathrm{m}`
:math:`L_p`      hydraulic permeability of the membrane          ``membrane_permeability``    0.01          :math:`\mathrm{m} \, \mathrm{h}^{-1} \, \mathrm{bar}^{-1}`
:math:`T`        temperature of the system                       ``temperature``              298           :math:`\mathrm{K}`
:math:`\chi`     concentration of surface charge on the membrane ``membrane_fixed_charge``    -140          :math:`\mathrm{mol} \, \mathrm{m}^{-3}`
================ =============================================== ============================ ============= ==========================================================

Variables
---------

The Two-Salt Diafiltration model adds the following variables.

=========================== ============================================================== ================================================= =========================================================================== =================================================================================================
Variable                    Description                                                    Name                                              Units                                                                       Indexed over
=========================== ============================================================== ================================================= =========================================================================== =================================================================================================
:math:`c_{i,d}`             ion concentration in the diafiltrate                           ``diafiltrate_conc_mol_comp``                     :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                                     :math:`t` and :math:`i \in \mathcal{I}`
:math:`c_{i,f}`             ion concentration in the feed                                  ``feed_conc_mol_comp``                            :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                                     :math:`t` and :math:`i \in \mathcal{I}`
:math:`c_{i,m}`             ion concentration in the membrane                              ``membrane_conc_mol_comp``                        :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                                     :math:`t`, math:`\bar{x}`, :math:`\bar{z}`, amd :math:`i \in \mathcal{I}`
:math:`c_{i,p}`             ion concentration in the permeate                              ``permeate_conc_mol_comp``                        :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                                     :math:`t`, :math:`\bar{x}`, amd :math:`i \in \mathcal{I}`
:math:`c_{i,r}`             ion concentration in the retentate                             ``retentate_conc_mol_comp``                       :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                                     :math:`t`, :math:`\bar{x}`, amd :math:`i \in \mathcal{I}`
:math:`\tilde{D}`           diffusion & convection coefficient denominator in the membrane ``membrane_D_tilde``                              :math:`\mathrm{mm}^2 \, \mathrm{h}^{-1} \, \mathrm{mol} \, \mathrm{m}^{-3}` :math:`t`, :math:`\bar{x}`, and :math:`\bar{z}`
:math:`D_{i,j}^{bilinear}`  bilinear cross-diffusion coefficient in the membrane           ``membrane_cross_diffusion_coefficient_bilinear`` :math:`\mathrm{mm}^4 \, \mathrm{h}^{-2} \, \mathrm{mol} \, \mathrm{m}^{-3}` :math:`t`, :math:`\bar{x}`, :math:`\bar{z}`, :math:`i \in \mathcal{I}`, :math:`j \in \mathcal{I}`
:math:`\alpha_i^{bilinear}` bilinear convection coefficient in the membrane                ``membrane_convection_coefficient_bilinear``      :math:`\mathrm{mm}^2 \, \mathrm{h}^{-1} \, \mathrm{mol} \, \mathrm{m}^{-3}` :math:`t`, :math:`\bar{x}`, :math:`\bar{z}`, and :math:`i \in \mathcal{I}`
:math:`D_{i,j}`             cross-diffusion coefficient in the membrane                    ``membrane_cross_diffusion_coefficient``          :math:`\mathrm{mm}^2 \, \mathrm{h}^{-1}`                                    :math:`t`, :math:`\bar{x}`, :math:`\bar{z}`, :math:`i \in \mathcal{I}`, :math:`j \in \mathcal{I}`
:math:`\alpha_i`            convection coefficient in the membrane                         ``membrnane_convection_coefficient``              :math:`\mathrm{dimensionless}`                                              :math:`t`, :math:`\bar{x}` and :math:`\bar{z}`
:math:`j_i`                 molar flux of ions across the membrane                         ``molar_ion_flux``                                :math:`\mathrm{mol} \, \mathrm{m}^{-2} \, \mathrm{h}^{-1}`                  :math:`t`, :math:`\bar{x}`, amd :math:`i \in \mathcal{I}`
:math:`J_w`                 water flux across the membrane                                 ``volume_flux_water``                             :math:`\mathrm{m}^3 \, \mathrm{m}^{-2} \, \mathrm{h}^{-1}`                  :math:`t` and :math:`\bar{x}`
:math:`L`                   length of the membrane                                         ``total_membrane_length``                         :math:`\mathrm{m}`
:math:`\Delta \pi`          osmotic pressure of feed-side fluid                            ``osmotic_pressure``                              :math:`\mathrm{bar}`                                                        :math:`t` and :math:`\bar{x}`
:math:`\Delta P`            applied pressure to the membrane                               ``applied_pressure``                              :math:`\mathrm{bar}`
:math:`q_d`                 volumetric flow rate of the diafiltrate                        ``diafiltrate_flow_volume``                       :math:`\mathrm{m}^3 \, \mathrm{h}^{-1}`                                     :math:`t`
:math:`q_f`                 volumetric flow rate of the feed                               ``feed_flow_volume``                              :math:`\mathrm{m}^3 \, \mathrm{h}^{-1}`                                     :math:`t`
:math:`q_p`                 volumetric flow rate of the permeate                           ``permeate_flow_volume``                          :math:`\mathrm{m}^3 \, \mathrm{h}^{-1}`                                     :math:`t` and :math:`\bar{x}`
:math:`q_r`                 volumetric flow rate of the retentate                          ``retentate_flow_volume``                         :math:`\mathrm{m}^3 \, \mathrm{h}^{-1}`                                     :math:`t` and :math:`\bar{x}`
:math:`w`                   length of the membrane module                                  ``total_module_length``                           :math:`\mathrm{m}`
=========================== ============================================================== ================================================= =========================================================================== =================================================================================================

Derivative Variables
--------------------

The Two-Salt Diafiltration model adds the following derivative variables.

=================================================== =========================================== ================================= ======================================= ==========================================================================
Variable                                            Description                                 Name                              Units                                   Indexed over
=================================================== =========================================== ================================= ======================================= ==========================================================================
:math:`\frac{\mathrm{d}c_{k,r}}{\mathrm{d}\bar{x}}` ion concentration gradient in the retentate ``d_retentate_conc_mass_comp_dx`` :math:`\mathrm{kg} \, \mathrm{m}^{-3}`  :math:`t`, :math:`\bar{x}`, and :math:`k \in \mathcal{K}`
:math:`\frac{\mathrm{d}q_r}{\mathrm{d}\bar{x}}`     retentate flow rate gradient                ``d_retentate_flow_volume_dx``    :math:`\mathrm{m}^3 \, \mathrm{h}^{-1}` :math:`t` and :math:`\bar{x}`
:math:`\frac{\partial c_{k,m}}{\partial \bar{z}}`   ion concentration gradient in the membrane  ``d_membrane_conc_mass_comp_dz``  :math:`\mathrm{kg} \, \mathrm{m}^{-3}`  :math:`t`, :math:`\bar{x}`, :math:`\bar{z}`, and :math:`k \in \mathcal{K}`
=================================================== =========================================== ================================= ======================================= ==========================================================================

Constraints
-----------

Differential mole balances:

.. math:: \frac{\mathrm{d}q_r(\bar{x})}{\mathrm{d}\bar{x}} = - J_w(\bar{x}) wL  \qquad \forall \, \bar{x} \neq 0
.. math:: q_r(\bar{x}) \frac{\mathrm{d}c_{\mathrm{Li},r}(\bar{x})}{\mathrm{d}\bar{x}} = wL (J_w(\bar{x}) c_{\mathrm{Li},r}(\bar{x}) - j_{\mathrm{Li}}(\bar{x}))  \qquad \forall \, \bar{x} \neq 0
.. math:: q_r(\bar{x}) \frac{\mathrm{d}c_{\mathrm{Co},r}(\bar{x})}{\mathrm{d}\bar{x}} = wL (J_w(\bar{x}) c_{\mathrm{Co},r}(\bar{x}) - j_{\mathrm{Co}}(\bar{x}))  \qquad \forall \, \bar{x} \neq 0

Bulk flux balances:

.. math:: q_p(\bar{x}) = \bar{x} wL J_w(\bar{x}) \qquad \forall \, \bar{x} \neq 0
.. math:: j_{\mathrm{Li}}(\bar{x}) = c_{\mathrm{Li},p}(\bar{x}) J_w(\bar{x}) \qquad \forall \, \bar{x} \neq 0
.. math:: j_{\mathrm{Co}}(\bar{x}) = c_{\mathrm{Co},p}(\bar{x}) J_w(\bar{x}) \qquad \forall \, \bar{x} \neq 0

Overall water flux through the membrane:

.. math:: J_w (\bar{x}) = L_p (\Delta P - \Delta \pi (\bar{x})) \qquad \forall \, \bar{x} \neq 0
.. math:: \Delta \pi (\bar{x}) = \mathrm{R} \mathrm{T} \sum_{i \in \mathcal{I}} n_i \sigma_i (c_{i,r}(\bar{x})-c_{i,p}(\bar{x})) \qquad \forall \, \bar{x} \neq 0

Solute flux through the membrane (extended Nernst-Planck equation):

.. math:: j_{\mathrm{Li}}(\bar{x}) = \alpha_{Li}(\bar{x},\bar{z}) c_{\mathrm{Li},m}(\bar{x},\bar{z}) J_w(\bar{x}) + \frac{D_{\mathrm{Li,Li}}(\bar{x},\bar{z})}{l} \frac{\partial c_{\mathrm{Li},m}(\bar{x},\bar{z})}{\partial \bar{z}} + \frac{D_{\mathrm{Li,Co}}(\bar{x},\bar{z})}{l} \frac{\partial c_{\mathrm{Co},m}(\bar{x},\bar{z})}{\partial \bar{z}} \qquad \forall \, \bar{z} \neq 0
.. math:: j_{\mathrm{Co}}(\bar{x}) = \alpha_{Co}(\bar{x},\bar{z}) c_{\mathrm{Co},m}(\bar{x},\bar{z}) J_w(\bar{x}) + \frac{D_{\mathrm{Co,Li}}(\bar{x},\bar{z})}{l} \frac{\partial c_{\mathrm{Li},m}(\bar{x},\bar{z})}{\partial \bar{z}} + \frac{D_{\mathrm{Co,Co}}(\bar{x},\bar{z})}{l} \frac{\partial c_{\mathrm{Co},m}(\bar{x},\bar{z})}{\partial \bar{z}} \qquad \forall \, \bar{z} \neq 0

Note that the single solute diffusion coefficients are provided in :math:`\mathrm{mm}^2\ \, \mathrm{h}^{-1}` to improve numerical stability, but the diffusion coefficients in the Nernst-Planck equations must be converted to :math:`\mathrm{m}^2\ \, \mathrm{h}^{-1}`. The convection and cross-diffusion coefficients are defined as:

.. math:: \tilde{D}(\bar{x},\bar{z}) = z_{\mathrm{Li}}(z_{\mathrm{Li}} D_{\mathrm{Li}} - z_{\mathrm{Cl}} D_{\mathrm{Cl}})c_{\mathrm{Li},m}(\bar{x},\bar{z}) + z_{\mathrm{Co}}(z_{\mathrm{Co}} D_{\mathrm{Co}} - z_{\mathrm{Cl}} D_{\mathrm{Cl}})c_{\mathrm{Co},m}(\bar{x},\bar{z}) - z_{\mathrm{Cl}} D_{\mathrm{Cl}} \chi \qquad \forall \, \bar{x} \neq 0
.. math:: \alpha_{\mathrm{Li}}(\bar{x},\bar{z}) = 1 + \frac{z_{\mathrm{Li}} D_{\mathrm{Li}} \chi}{\tilde{D}(\bar{x},\bar{z})} \qquad \forall \, \bar{x} \neq 0
.. math:: \alpha_{\mathrm{Co}}(\bar{x},\bar{z}) = 1 + \frac{z_{\mathrm{Co}} D_{\mathrm{Co}} \chi}{\tilde{D}(\bar{x},\bar{z})} \qquad \forall \, \bar{x} \neq 0
.. math:: D_{\mathrm{Li,Li}}(\bar{x},\bar{z}) = \frac{z_{\mathrm{Li}} D_{\mathrm{Li}} D_{\mathrm{Cl}}(z_{\mathrm{Cl}} - z_{\mathrm{Li}})c_{\mathrm{Li},m}(\bar{x},\bar{z}) + z_{\mathrm{Co}} D_{\mathrm{Li}}(z_{\mathrm{Cl}} D_{\mathrm{Cl}} - z_{\mathrm{Co}} D_{\mathrm{Co}})c_{\mathrm{Co},m}(\bar{x},\bar{z}) + z_{\mathrm{Cl}} D_{\mathrm{Li}} D_{\mathrm{Cl}} \chi}{\tilde{D}(\bar{x},\bar{z})} \qquad \forall \, \bar{x} \neq 0
.. math:: D_{\mathrm{Li,Co}}(\bar{x},\bar{z}) = \frac{z_{\mathrm{Li}} z_{\mathrm{Co}} D_{\mathrm{Li}}(D_{\mathrm{Co}} - D_{\mathrm{Cl}})c_{\mathrm{Li},m}(\bar{x},\bar{z})}{\tilde{D}(\bar{x},\bar{z})} \qquad \forall \, \bar{x} \neq 0
.. math:: D_{\mathrm{Co,Li}}(\bar{x},\bar{z}) = \frac{z_{\mathrm{Li}} z_{\mathrm{Co}} D_{\mathrm{Co}}(D_{\mathrm{Li}} - D_{\mathrm{Cl}})c_{\mathrm{Co},m}(\bar{x},\bar{z})}{\tilde{D}(\bar{x},\bar{z})} \qquad \forall \, \bar{x} \neq 0
.. math:: D_{\mathrm{Co,Co}}(\bar{x},\bar{z}) = \frac{z_{\mathrm{Li}} D_{\mathrm{Co}} (z_{\mathrm{Cl}} D_{\mathrm{Cl}} - z_{\mathrm{Li}} D_{\mathrm{Li}})c_{\mathrm{Li},m}(\bar{x},\bar{z}) + z_{\mathrm{Co}} D_{\mathrm{Co}} D_{\mathrm{Cl}} (z_{\mathrm{Cl}} - z_{\mathrm{Co}})c_{\mathrm{Co},m}(\bar{x},\bar{z}) + z_{\mathrm{Cl}} D_{\mathrm{Co}} D_{\mathrm{Cl}} \chi}{\tilde{D}(\bar{x},\bar{z})} \qquad \forall \, \bar{x} \neq 0

The diffusion and convection coefficients are reformulated to bilinear constraints:

.. math:: \alpha_{\mathrm{Li}}^{bilinear}(\bar{x},\bar{z}) = \alpha_{\mathrm{Li}}(\bar{x},\bar{z}) \tilde{D}(\bar{x},\bar{z}) = \tilde{D}(\bar{x},\bar{z}) + z_{\mathrm{Li}} D_{\mathrm{Li}} \chi \qquad \forall \, \bar{x} \neq 0
.. math:: \alpha_{\mathrm{Co}}^{bilinear}(\bar{x},\bar{z}) = \alpha_{\mathrm{Co}}(\bar{x},\bar{z}) \tilde{D}(\bar{x},\bar{z}) = \tilde{D}(\bar{x},\bar{z}) + z_{\mathrm{Co}} D_{\mathrm{Co}} \chi \qquad \forall \, \bar{x} \neq 0
.. math:: D_{\mathrm{Li,Li}}^{bilinear}(\bar{x},\bar{z}) = D_{\mathrm{Li,Li}}(\bar{x},\bar{z}) \tilde{D}(\bar{x},\bar{z}) = z_{\mathrm{Li}} D_{\mathrm{Li}} D_{\mathrm{Cl}}(z_{\mathrm{Cl}} - z_{\mathrm{Li}})c_{\mathrm{Li},m}(\bar{x},\bar{z}) + z_{\mathrm{Co}} D_{\mathrm{Li}}(z_{\mathrm{Cl}} D_{\mathrm{Cl}} - z_{\mathrm{Co}} D_{\mathrm{Co}})c_{\mathrm{Co},m}(\bar{x},\bar{z}) + z_{\mathrm{Cl}} D_{\mathrm{Li}} D_{\mathrm{Cl}} \chi \qquad \forall \, \bar{x} \neq 0
.. math:: D_{\mathrm{Li,Co}}^{bilinear}(\bar{x},\bar{z}) = D_{\mathrm{Li,Co}}(\bar{x},\bar{z}) \tilde{D}(\bar{x},\bar{z}) = z_{\mathrm{Li}} z_{\mathrm{Co}} D_{\mathrm{Li}}(D_{\mathrm{Co}} - D_{\mathrm{Cl}})c_{\mathrm{Li},m}(\bar{x},\bar{z}) \qquad \forall \, \bar{x} \neq 0
.. math:: D_{\mathrm{Co,Li}}^{bilinear}(\bar{x},\bar{z}) = D_{\mathrm{Co,Li}}(\bar{x},\bar{z}) \tilde{D}(\bar{x},\bar{z}) = z_{\mathrm{Li}} z_{\mathrm{Co}} D_{\mathrm{Co}}(D_{\mathrm{Li}} - D_{\mathrm{Cl}})c_{\mathrm{Co},m}(\bar{x},\bar{z}) \qquad \forall \, \bar{x} \neq 0
.. math:: D_{\mathrm{Co,Co}}^{bilinear}(\bar{x},\bar{z}) = D_{\mathrm{Co,Co}}(\bar{x},\bar{z}) \tilde{D}(\bar{x},\bar{z}) = z_{\mathrm{Li}} D_{\mathrm{Co}} (z_{\mathrm{Cl}} D_{\mathrm{Cl}} - z_{\mathrm{Li}} D_{\mathrm{Li}})c_{\mathrm{Li},m}(\bar{x},\bar{z}) + z_{\mathrm{Co}} D_{\mathrm{Co}} D_{\mathrm{Cl}} (z_{\mathrm{Cl}} - z_{\mathrm{Co}})c_{\mathrm{Co},m}(\bar{x},\bar{z}) + z_{\mathrm{Cl}} D_{\mathrm{Co}} D_{\mathrm{Cl}} \chi \qquad \forall \, \bar{x} \neq 0

No applied potential on the system:

.. math:: 0 = \sum_{i \in \mathcal{I}} z_i j_i(\bar{x}) \qquad \forall \, \bar{x} \neq 0

Electroneutrality:

.. math:: 0 = \sum_{i \in \mathcal{I}} z_i c_{i,r}(\bar{x})
.. math:: 0 = \chi + \sum_{i \in \mathcal{I}} z_i c_{i,m}(\bar{x},\bar{z}) \qquad \forall \, \bar{z} \neq 0
.. math:: 0 = \sum_{i \in \mathcal{I}} z_i c_{i,p}(\bar{x})

Partitioning at the retentate-membrane interface:

.. math:: H_{\mathrm{Li}} H_{\mathrm{Cl}} c_{\mathrm{Li},r}(\bar{x}) c_{\mathrm{Cl},r}(\bar{x}) = c_{\mathrm{Li},m}(\bar{x},\bar{z}=0) c_{\mathrm{Cl},m}(\bar{x},\bar{z}=0) \qquad \forall \, \bar{x} \neq 0
.. math:: H_{\mathrm{Co}} H_{\mathrm{Cl}}^{z_{\mathrm{Co}}} c_{\mathrm{Co},r}(\bar{x}) c_{\mathrm{Cl},r}(\bar{x})^{z_{\mathrm{Co}}} =c_{\mathrm{Co},m}(\bar{x},\bar{z}=0) c_{\mathrm{Cl},m}(\bar{x},\bar{z}=0)^{z_{\mathrm{Co}}} \qquad \forall \, \bar{x} \neq 0

Partitioning at the membrane-permeate interface:

.. math:: H_{\mathrm{Li}} H_{\mathrm{Cl}} c_{\mathrm{Li},p}(\bar{x}) c_{\mathrm{Cl},p}(\bar{x}) = c_{\mathrm{Li},m}(\bar{x},\bar{z}=1) c_{\mathrm{Cl},m}(\bar{x},\bar{z}=1) \qquad \forall \, \bar{x} \neq 0
.. math:: H_{\mathrm{Co}} H_{\mathrm{Cl}}^{z_{\mathrm{Co}}} c_{\mathrm{Co},p}(\bar{x}) c_{\mathrm{Cl},p}(\bar{x})^{z_{\mathrm{Co}}} =c_{\mathrm{Co},m}(\bar{x},\bar{z}=1) c_{\mathrm{Cl},m}(\bar{x},\bar{z}=1)^{z_{\mathrm{Co}}} \qquad \forall \, \bar{x} \neq 0

The following boundary conditions complete the model:

.. math:: q_r(\bar{x}=0) = q_f + q_d
.. math:: c_{\mathrm{Li},r}(\bar{x}=0) = \frac{q_f c_{\mathrm{Li},f} + q_d c_{\mathrm{Li},d}}{q_f + q_d}
.. math:: c_{\mathrm{Co},r}(\bar{x}=0) = \frac{q_f c_{\mathrm{Co},f} + q_d c_{\mathrm{Co},d}}{q_f + q_d}

The following boundary conditions (which are expected to be zero) are fixed to improve numerical stability (with the appropriate constraints deactivated as described above):

.. math:: q_p(\bar{x}=0) = \epsilon
.. math:: c_{\mathrm{Li},p}(\bar{x}=0) = \epsilon
.. math:: c_{\mathrm{Co},p}(\bar{x}=0) = \epsilon
.. math:: c_{\mathrm{Li},m} (\bar{x}=0,\bar{z}) = \epsilon \qquad \forall \, \bar{z}
.. math:: c_{\mathrm{Co},m} (\bar{x}=0,\bar{z}) = \epsilon \qquad \forall \, \bar{z}
.. math:: c_{\mathrm{Cl},m} (\bar{x}=0,\bar{z}) = \epsilon \qquad \forall \, \bar{z}
.. math:: \frac{\mathrm{d}q_r(\bar{x})}{\mathrm{d}\bar{x}}(\bar{x}=0)=\epsilon
.. math:: \frac{\mathrm{d}c_{\mathrm{Li},r}(\bar{x})}{\mathrm{d}\bar{x}}(\bar{x}=0)=\epsilon
.. math:: \frac{\mathrm{d}c_{\mathrm{Co},r}(\bar{x})}{\mathrm{d}\bar{x}}(\bar{x}=0)=\epsilon
.. math:: J_w(\bar{x}=0) = \epsilon
.. math:: j_{\mathrm{Li}}(\bar{x}=0) = \epsilon
.. math:: j_{\mathrm{Co}}(\bar{x}=0) = \epsilon
.. math:: j_{\mathrm{Cl}}(\bar{x}=0) = \epsilon
"""

from pyomo.common.config import ConfigBlock, ConfigValue
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.environ import (
    Constraint,
    Param,
    Reference,
    Set,
    Suffix,
    TransformationFactory,
    Var,
    units,
    value,
)
from pyomo.network import Port

from idaes.core import UnitModelBlockData, declare_process_block_class, useDefault
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.constants import Constants


@declare_process_block_class("TwoSaltDiafiltration")
class TwoSaltDiafiltrationData(UnitModelBlockData):
    """
    Two-Salt Diafiltration Unit Model Class.
    """

    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for membrane system",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}
""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {see property package for documentation}
""",
        ),
    )
    CONFIG.declare(
        "NFE_module_length",
        ConfigValue(
            doc="Number of discretization points in the x-direction (across module length)",
        ),
    )
    CONFIG.declare(
        "NFE_membrane_thickness",
        ConfigValue(
            doc="Number of discretization points in the z-direction (across membrane thickness)",
        ),
    )

    def build(self):
        """
        Build method for the two salt diafiltration unit model
        """
        super().build()

        # TODO: generalize to any 2 cations and 1 anion

        self.add_mutable_parameters()
        self.add_variables()
        self.add_constraints()
        self.discretize_model()
        self.fix_boundary_values()
        self.add_scaling_factors()
        self.add_ports()

    def add_mutable_parameters(self):
        """
        Adds default parameters for the two salt diafiltration unit model.

        Values can be changed by the user during implementation.

        Assumes membrane thickness of 100 nm.

        Membrane permeability and fixed charged are estimated from:
        Liu, Xinhong, et al. (2025) https://doi.org/10.1021/acs.iecr.4c04763
        """
        self.numerical_zero_tolerance = Param(
            initialize=1e-10,
            mutable=True,
            doc="Numerical tolerance for zero values in the model",
        )
        self.total_membrane_thickness = Param(
            initialize=1e-7,
            mutable=True,
            units=units.m,
            doc="Thickness of membrane (z-direction)",
        )
        self.membrane_fixed_charge = Param(
            initialize=-140,
            mutable=True,
            units=units.mol / units.m**3,  # mM
            doc="Fixed charge on the membrane",
        )
        self.membrane_permeability = Param(
            initialize=0.01,
            mutable=True,
            units=units.m / units.h / units.bar,
            doc="Hydraulic permeability coefficient",
        )
        self.temperature = Param(
            initialize=298,
            mutable=True,
            units=units.K,
            doc="System temperature",
        )

    def add_variables(self):
        """
        Adds variables for the two salt diafiltration unit model.

        Membrane module dimensions and maximum flowrate (17 m3/h) are
        estimated from NF270-440 modules.

        Assumes 4 modules in series.
        """
        # define length scales
        self.dimensionless_module_length = ContinuousSet(bounds=(0, 1))
        self.dimensionless_membrane_thickness = ContinuousSet(bounds=(0, 1))

        # add a time index since the property package variables are indexed over time
        self.time = Set(initialize=[0])

        # add components
        self.solutes = Set(initialize=["Li", "Co", "Cl"])
        self.cations = Set(initialize=["Li", "Co"])

        # add global variables
        self.total_module_length = Var(
            initialize=4,  # 4 tubes that are ~1m long each (NF270-440)
            units=units.m,
            bounds=[1e-11, None],
            doc="Width of the membrane (x-direction)",
        )
        self.total_membrane_length = Var(
            initialize=41,  # 41 m of length in each tube (NF270-440)
            units=units.m,
            bounds=[1e-11, None],
            doc="Length of the membrane, wound radially",
        )
        self.applied_pressure = Var(
            initialize=10,
            units=units.bar,
            bounds=[1e-11, 41],  # maximum operating presssure (NF270-440)
            doc="Pressure applied to membrane",
        )
        self.feed_flow_volume = Var(
            self.time,
            initialize=12.5,
            units=units.m**3 / units.h,
            bounds=[1e-11, None],
            doc="Volumetric flow rate of the feed",
        )

        def initialize_feed_conc_mol_comp(m, t, j):
            vals = {"Li": 245, "Co": 288, "Cl": 822}
            return vals[j]

        self.feed_conc_mol_comp = Var(
            self.time,
            self.solutes,
            initialize=initialize_feed_conc_mol_comp,
            units=units.mol / units.m**3,  # mM
            bounds=[1e-11, None],
            doc="Mole concentration of solutes in the feed",
        )
        self.diafiltrate_flow_volume = Var(
            self.time,
            initialize=3.75,
            units=units.m**3 / units.h,
            bounds=[1e-11, None],
            doc="Volumetric flow rate of the diafiltrate",
        )

        def initialize_diafiltrate_conc_mol_comp(m, t, j):
            vals = {"Li": 14, "Co": 3, "Cl": 21}
            return vals[j]

        self.diafiltrate_conc_mol_comp = Var(
            self.time,
            self.solutes,
            initialize=initialize_diafiltrate_conc_mol_comp,
            units=units.mol / units.m**3,  # mM
            bounds=[1e-11, None],
            doc="Mole concentration of solutes in the diafiltrate",
        )

        # add variables dependent on dimensionless_module_length
        self.volume_flux_water = Var(
            self.time,
            self.dimensionless_module_length,
            initialize=0.06,
            units=units.m**3 / units.m**2 / units.h,
            bounds=[1e-11, None],
            doc="Volumetric water flux of water across the membrane",
        )

        def initialize_mol_flux_comp(m, t, w, j):
            vals = {"Li": 11, "Co": 13, "Cl": 37}
            return vals[j]

        self.molar_ion_flux = Var(
            self.time,
            self.dimensionless_module_length,
            self.solutes,
            initialize=initialize_mol_flux_comp,
            units=units.mol / units.m**2 / units.h,
            bounds=[1e-11, None],
            doc="Mole flux of solutes across the membrane (z-direction, x-dependent)",
        )
        self.retentate_flow_volume = Var(
            self.time,
            self.dimensionless_module_length,
            initialize=6.7,
            units=units.m**3 / units.h,
            bounds=[1e-11, None],
            doc="Volumetric flow rate of the retentate, x-dependent",
        )

        def initialize_retentate_conc_mol_comp(m, t, w, j):
            vals = {"Li": 198, "Co": 241, "Cl": 680}
            return vals[j]

        self.retentate_conc_mol_comp = Var(
            self.time,
            self.dimensionless_module_length,
            self.solutes,
            initialize=initialize_retentate_conc_mol_comp,
            units=units.mol / units.m**3,  # mM
            bounds=[1e-11, None],
            doc="Mole concentration of solutes in the retentate, x-dependent",
        )
        self.permeate_flow_volume = Var(
            self.time,
            self.dimensionless_module_length,
            initialize=9.5,
            units=units.m**3 / units.h,
            bounds=[1e-11, None],
            doc="Volumetric flow rate of the permeate, x-dependent",
        )

        def initialize_permeate_conc_mol_comp(m, t, w, j):
            vals = {"Li": 191, "Co": 220, "Cl": 632}
            return vals[j]

        self.permeate_conc_mol_comp = Var(
            self.time,
            self.dimensionless_module_length,
            self.solutes,
            initialize=initialize_permeate_conc_mol_comp,
            units=units.mol / units.m**3,  # mM
            bounds=[1e-11, None],
            doc="Mole concentration of solutes in the permeate, x-dependent",
        )
        self.osmotic_pressure = Var(
            self.time,
            self.dimensionless_module_length,
            initialize=4,
            units=units.bar,
            bounds=[1e-11, None],
            doc="Osmostic pressure difference across the membrane",
        )

        # add variables dependent on dimensionless_module_length and dimensionless_membrane_thickness
        def initialize_membrane_conc_mol_comp(m, t, w, l, j):
            vals = {"Li": 109, "Co": 18, "Cl": 4.5}
            return vals[j]

        self.membrane_conc_mol_comp = Var(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            self.solutes,
            initialize=initialize_membrane_conc_mol_comp,
            units=units.mol / units.m**3,  # mM
            bounds=[1e-11, None],
            doc="Mole concentration of solutes in the membrane, x- and z-dependent",
        )
        self.membrane_D_tilde = Var(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            initialize=620,
            units=(units.mm**2 / units.hr) * (units.mol / units.m**3),  # D * c
            doc="Denominator of diffusion and convection coefficients in membrane",
        )

        def initialize_membrane_cross_diffusion_coefficient_bilinear(m, t, w, l, j, k):
            vals = {"Li": {"Li": -3800, "Co": -3800}, "Co": {"Li": -340, "Co": -2500}}
            return vals[j][k]

        self.membrane_cross_diffusion_coefficient_bilinear = Var(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            self.cations,
            self.cations,
            initialize=initialize_membrane_cross_diffusion_coefficient_bilinear,
            units=(units.mm**2 / units.h)
            * (units.mm**2 / units.h * units.mol / units.m**3),  # D * D,tilde
            doc="Bi-linear cross diffusion coefficients for cations in membrane",
        )

        def initialize_membrane_convection_coefficient_bilinear(m, t, w, l, j):
            vals = {"Li": 105, "Co": -115}
            return vals[j]

        self.membrane_convection_coefficient_bilinear = Var(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            self.cations,
            initialize=initialize_membrane_convection_coefficient_bilinear,
            units=(units.mm**2 / units.hr) * (units.mol / units.m**3),  # D,tilde
            doc="Convection coefficients for cations in membrane",
        )

        def initialize_membrane_cross_diffusion_coefficient(m, t, w, l, j, k):
            vals = {"Li": {"Li": -6, "Co": -6}, "Co": {"Li": -0.5, "Co": -4}}
            return vals[j][k]

        self.membrane_cross_diffusion_coefficient = Var(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            self.cations,
            self.cations,
            initialize=initialize_membrane_cross_diffusion_coefficient,
            units=units.mm**2 / units.h,
            doc="Cross diffusion coefficient for cations in membrane",
        )

        def initialize_membrane_convection_coefficient(m, t, w, l, j):
            vals = {"Li": 0.17, "Co": -0.18}
            return vals[j]

        self.membrane_convection_coefficient = Var(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            self.cations,
            initialize=initialize_membrane_convection_coefficient,
            units=units.dimensionless,
            doc="Convection coefficients for cations in membrane",
        )

        # define the (partial) derivative variables
        self.d_retentate_conc_mol_comp_dx = DerivativeVar(
            self.retentate_conc_mol_comp,
            wrt=self.dimensionless_module_length,
            units=units.mol / units.m**3,  # mM
            doc="Solute concentration gradient in the retentate",
        )
        self.d_retentate_flow_volume_dx = DerivativeVar(
            self.retentate_flow_volume,
            wrt=self.dimensionless_module_length,
            units=units.m**3 / units.h,
            doc="Volume flow gradient in the retentate",
        )
        self.d_membrane_conc_mol_comp_dz = DerivativeVar(
            self.membrane_conc_mol_comp,
            wrt=self.dimensionless_membrane_thickness,
            units=units.mol / units.m**3,  # mM
            doc="Solute concentration gradient wrt membrane thickness",
        )

    def add_constraints(self):
        """
        Adds model constraints for the two salt diafiltration unit model.
        """

        # mol balance constraints
        def _overall_mol_balance(blk, x):
            if x == 0:
                return Constraint.Skip
            return blk.d_retentate_flow_volume_dx[0, x] == (
                -blk.volume_flux_water[0, x]
                * blk.total_membrane_length
                * blk.total_module_length
            )

        self.overall_mol_balance = Constraint(
            self.dimensionless_module_length, rule=_overall_mol_balance
        )

        def _lithium_mol_balance(blk, x):
            if x == 0:
                return Constraint.Skip
            return (
                blk.retentate_flow_volume[0, x]
                * blk.d_retentate_conc_mol_comp_dx[0, x, "Li"]
            ) == (
                (
                    blk.volume_flux_water[0, x]
                    * blk.retentate_conc_mol_comp[0, x, "Li"]
                    - blk.molar_ion_flux[0, x, "Li"]
                )
                * blk.total_membrane_length
                * blk.total_module_length
            )

        self.lithium_mol_balance = Constraint(
            self.dimensionless_module_length, rule=_lithium_mol_balance
        )

        def _cobalt_mol_balance(blk, x):
            if x == 0:
                return Constraint.Skip
            return (
                blk.retentate_flow_volume[0, x]
                * blk.d_retentate_conc_mol_comp_dx[0, x, "Co"]
            ) == (
                (
                    blk.volume_flux_water[0, x]
                    * blk.retentate_conc_mol_comp[0, x, "Co"]
                    - blk.molar_ion_flux[0, x, "Co"]
                )
                * blk.total_membrane_length
                * blk.total_module_length
            )

        self.cobalt_mol_balance = Constraint(
            self.dimensionless_module_length, rule=_cobalt_mol_balance
        )

        # bulk flux balance constraints
        def _bulk_flux_equation_overall(blk, x):
            if x == 0:
                return Constraint.Skip
            return (
                blk.permeate_flow_volume[0, x]
                == blk.volume_flux_water[0, x]
                * x
                * blk.total_membrane_length
                * blk.total_module_length
            )

        self.bulk_flux_equation_overall = Constraint(
            self.dimensionless_module_length, rule=_bulk_flux_equation_overall
        )

        def _bulk_flux_equation_lithium(blk, x):
            if x == 0:
                return Constraint.Skip
            return blk.molar_ion_flux[0, x, "Li"] == (
                blk.permeate_conc_mol_comp[0, x, "Li"] * blk.volume_flux_water[0, x]
            )

        self.bulk_flux_equation_lithium = Constraint(
            self.dimensionless_module_length, rule=_bulk_flux_equation_lithium
        )

        def _bulk_flux_equation_cobalt(blk, x):
            if x == 0:
                return Constraint.Skip
            return blk.molar_ion_flux[0, x, "Co"] == (
                blk.permeate_conc_mol_comp[0, x, "Co"] * blk.volume_flux_water[0, x]
            )

        self.bulk_flux_equation_cobalt = Constraint(
            self.dimensionless_module_length, rule=_bulk_flux_equation_cobalt
        )

        # transport constraints (first principles)
        def _lumped_water_flux(blk, x):
            if x == 0:
                return Constraint.Skip
            return blk.volume_flux_water[0, x] == (
                blk.membrane_permeability
                * (blk.applied_pressure - blk.osmotic_pressure[0, x])
            )

        self.lumped_water_flux = Constraint(
            self.dimensionless_module_length, rule=_lumped_water_flux
        )

        def _membrane_D_tilde_calculation(blk, x, z):
            if x == 0:
                return Constraint.Skip
            return blk.membrane_D_tilde[0, x, z] == (
                (
                    (
                        (
                            (blk.config.property_package.charge["Li"] ** 2)
                            * blk.config.property_package.diffusion_coefficient["Li"]
                        )
                        - (
                            blk.config.property_package.charge["Li"]
                            * blk.config.property_package.charge["Cl"]
                            * blk.config.property_package.diffusion_coefficient["Cl"]
                        )
                    )
                    * blk.membrane_conc_mol_comp[0, x, z, "Li"]
                )
                + (
                    (
                        (
                            (blk.config.property_package.charge["Co"] ** 2)
                            * blk.config.property_package.diffusion_coefficient["Co"]
                        )
                        - (
                            blk.config.property_package.charge["Co"]
                            * blk.config.property_package.charge["Cl"]
                            * blk.config.property_package.diffusion_coefficient["Cl"]
                        )
                    )
                    * blk.membrane_conc_mol_comp[0, x, z, "Co"]
                )
                - (
                    blk.config.property_package.charge["Cl"]
                    * blk.config.property_package.diffusion_coefficient["Cl"]
                    * blk.membrane_fixed_charge
                )
            )

        self.membrane_D_tilde_calculation = Constraint(
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            rule=_membrane_D_tilde_calculation,
        )

        def _membrane_diffusion_coefficient_bilinear_calculation(blk, x, z, j, k):
            if x == 0:
                return Constraint.Skip
            return (
                blk.membrane_cross_diffusion_coefficient_bilinear[0, x, z, j, k]
                == blk.membrane_cross_diffusion_coefficient[0, x, z, j, k]
                * blk.membrane_D_tilde[0, x, z]
            )

        self.membrane_diffusion_coefficient_bilinear_calculation = Constraint(
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            self.cations,
            self.cations,
            rule=_membrane_diffusion_coefficient_bilinear_calculation,
        )

        def _membrane_convection_coefficient_bilinear_calculation(blk, x, z, j):
            if x == 0:
                return Constraint.Skip
            return (
                blk.membrane_convection_coefficient_bilinear[0, x, z, j]
                == blk.membrane_convection_coefficient[0, x, z, j]
                * blk.membrane_D_tilde[0, x, z]
            )

        self.membrane_convection_coefficient_bilinear_calculation = Constraint(
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            self.cations,
            rule=_membrane_convection_coefficient_bilinear_calculation,
        )

        def _membrane_diffusion_coefficient_lithium_lithium_calculation(blk, x, z):
            if x == 0:
                return Constraint.Skip
            return blk.membrane_cross_diffusion_coefficient_bilinear[
                0, x, z, "Li", "Li"
            ] == (
                (
                    (
                        blk.config.property_package.charge["Li"]
                        * blk.config.property_package.charge["Cl"]
                        * blk.config.property_package.diffusion_coefficient["Li"]
                        * blk.config.property_package.diffusion_coefficient["Cl"]
                    )
                    - (
                        (blk.config.property_package.charge["Li"] ** 2)
                        * blk.config.property_package.diffusion_coefficient["Li"]
                        * blk.config.property_package.diffusion_coefficient["Cl"]
                    )
                )
                * blk.membrane_conc_mol_comp[0, x, z, "Li"]
                + (
                    (
                        blk.config.property_package.charge["Co"]
                        * blk.config.property_package.charge["Cl"]
                        * blk.config.property_package.diffusion_coefficient["Li"]
                        * blk.config.property_package.diffusion_coefficient["Cl"]
                    )
                    - (
                        (blk.config.property_package.charge["Co"] ** 2)
                        * blk.config.property_package.diffusion_coefficient["Li"]
                        * blk.config.property_package.diffusion_coefficient["Co"]
                    )
                )
                * blk.membrane_conc_mol_comp[0, x, z, "Co"]
                + (
                    blk.config.property_package.charge["Cl"]
                    * blk.config.property_package.diffusion_coefficient["Li"]
                    * blk.config.property_package.diffusion_coefficient["Cl"]
                    * blk.membrane_fixed_charge
                )
            )

        self.membrane_diffusion_coefficient_lithium_lithium_calculation = Constraint(
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            rule=_membrane_diffusion_coefficient_lithium_lithium_calculation,
        )

        def _membrane_diffusion_coefficient_lithium_cobalt_calculation(blk, x, z):
            if x == 0:
                return Constraint.Skip
            return blk.membrane_cross_diffusion_coefficient_bilinear[
                0, x, z, "Li", "Co"
            ] == (
                (
                    (
                        blk.config.property_package.charge["Li"]
                        * blk.config.property_package.charge["Co"]
                        * blk.config.property_package.diffusion_coefficient["Li"]
                        * blk.config.property_package.diffusion_coefficient["Co"]
                    )
                    - (
                        blk.config.property_package.charge["Li"]
                        * blk.config.property_package.charge["Co"]
                        * blk.config.property_package.diffusion_coefficient["Li"]
                        * blk.config.property_package.diffusion_coefficient["Cl"]
                    )
                )
                * blk.membrane_conc_mol_comp[0, x, z, "Li"]
            )

        self.membrane_diffusion_coefficient_lithium_cobalt_calculation = Constraint(
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            rule=_membrane_diffusion_coefficient_lithium_cobalt_calculation,
        )

        def _membrane_diffusion_coefficient_cobalt_lithium_calculation(blk, x, z):
            if x == 0:
                return Constraint.Skip
            return blk.membrane_cross_diffusion_coefficient_bilinear[
                0, x, z, "Co", "Li"
            ] == (
                (
                    (
                        blk.config.property_package.charge["Li"]
                        * blk.config.property_package.charge["Co"]
                        * blk.config.property_package.diffusion_coefficient["Li"]
                        * blk.config.property_package.diffusion_coefficient["Co"]
                    )
                    - (
                        blk.config.property_package.charge["Li"]
                        * blk.config.property_package.charge["Co"]
                        * blk.config.property_package.diffusion_coefficient["Co"]
                        * blk.config.property_package.diffusion_coefficient["Cl"]
                    )
                )
                * blk.membrane_conc_mol_comp[0, x, z, "Co"]
            )

        self.membrane_diffusion_coefficient_cobalt_lithium_calculation = Constraint(
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            rule=_membrane_diffusion_coefficient_cobalt_lithium_calculation,
        )

        def _membrane_diffusion_coefficient_cobalt_cobalt_calculation(blk, x, z):
            if x == 0:
                return Constraint.Skip
            return blk.membrane_cross_diffusion_coefficient_bilinear[
                0, x, z, "Co", "Co"
            ] == (
                (
                    (
                        blk.config.property_package.charge["Li"]
                        * blk.config.property_package.charge["Cl"]
                        * blk.config.property_package.diffusion_coefficient["Co"]
                        * blk.config.property_package.diffusion_coefficient["Cl"]
                    )
                    - (
                        (blk.config.property_package.charge["Li"] ** 2)
                        * blk.config.property_package.diffusion_coefficient["Li"]
                        * blk.config.property_package.diffusion_coefficient["Co"]
                    )
                )
                * blk.membrane_conc_mol_comp[0, x, z, "Li"]
                + (
                    (
                        blk.config.property_package.charge["Co"]
                        * blk.config.property_package.charge["Cl"]
                        * blk.config.property_package.diffusion_coefficient["Co"]
                        * blk.config.property_package.diffusion_coefficient["Cl"]
                    )
                    - (
                        (blk.config.property_package.charge["Co"] ** 2)
                        * blk.config.property_package.diffusion_coefficient["Co"]
                        * blk.config.property_package.diffusion_coefficient["Cl"]
                    )
                )
                * blk.membrane_conc_mol_comp[0, x, z, "Co"]
                + (
                    blk.config.property_package.charge["Cl"]
                    * blk.config.property_package.diffusion_coefficient["Co"]
                    * blk.config.property_package.diffusion_coefficient["Cl"]
                    * blk.membrane_fixed_charge
                )
            )

        self.membrane_diffusion_coefficient_cobalt_cobalt_calculation = Constraint(
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            rule=_membrane_diffusion_coefficient_cobalt_cobalt_calculation,
        )

        def _membrane_convection_coefficient_lithium_calculation(blk, x, z):
            if x == 0:
                return Constraint.Skip
            return blk.membrane_convection_coefficient_bilinear[0, x, z, "Li"] == (
                blk.membrane_D_tilde[0, x, z]
                + (
                    blk.config.property_package.charge["Li"]
                    * blk.config.property_package.diffusion_coefficient["Li"]
                    * blk.membrane_fixed_charge
                )
            )

        self.membrane_convection_coefficient_lithium_calculation = Constraint(
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            rule=_membrane_convection_coefficient_lithium_calculation,
        )

        def _membrane_convection_coefficient_cobalt_calculation(blk, x, z):
            if x == 0:
                return Constraint.Skip
            return blk.membrane_convection_coefficient_bilinear[0, x, z, "Co"] == (
                blk.membrane_D_tilde[0, x, z]
                + (
                    blk.config.property_package.charge["Co"]
                    * blk.config.property_package.diffusion_coefficient["Co"]
                    * blk.membrane_fixed_charge
                )
            )

        self.membrane_convection_coefficient_cobalt_calculation = Constraint(
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            rule=_membrane_convection_coefficient_cobalt_calculation,
        )

        def _lithium_flux_membrane(blk, x, z):
            if x == 0:
                return Constraint.Skip
            return blk.molar_ion_flux[0, x, "Li"] == (
                (
                    blk.membrane_convection_coefficient[0, x, z, "Li"]
                    * blk.membrane_conc_mol_comp[0, x, z, "Li"]
                    * blk.volume_flux_water[0, x]
                )
                + (
                    units.convert(
                        blk.membrane_cross_diffusion_coefficient[0, x, z, "Li", "Li"],
                        to_units=units.m**2 / units.h,
                    )
                    / blk.total_membrane_thickness
                    * blk.d_membrane_conc_mol_comp_dz[0, x, z, "Li"]
                )
                + (
                    units.convert(
                        blk.membrane_cross_diffusion_coefficient[0, x, z, "Li", "Co"],
                        to_units=units.m**2 / units.h,
                    )
                    / blk.total_membrane_thickness
                    * blk.d_membrane_conc_mol_comp_dz[0, x, z, "Co"]
                )
            )

        self.lithium_flux_membrane = Constraint(
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            rule=_lithium_flux_membrane,
        )

        def _cobalt_flux_membrane(blk, x, z):
            if x == 0:
                return Constraint.Skip
            return blk.molar_ion_flux[0, x, "Co"] == (
                (
                    blk.membrane_convection_coefficient[0, x, z, "Co"]
                    * blk.membrane_conc_mol_comp[0, x, z, "Co"]
                    * blk.volume_flux_water[0, x]
                )
                + (
                    units.convert(
                        blk.membrane_cross_diffusion_coefficient[0, x, z, "Co", "Li"],
                        to_units=units.m**2 / units.h,
                    )
                    / blk.total_membrane_thickness
                    * blk.d_membrane_conc_mol_comp_dz[0, x, z, "Li"]
                )
                + (
                    units.convert(
                        blk.membrane_cross_diffusion_coefficient[0, x, z, "Co", "Co"],
                        to_units=units.m**2 / units.h,
                    )
                    / blk.total_membrane_thickness
                    * blk.d_membrane_conc_mol_comp_dz[0, x, z, "Co"]
                )
            )

        self.cobalt_flux_membrane = Constraint(
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            rule=_cobalt_flux_membrane,
        )

        def _chloride_flux_membrane(blk, x):
            if x == 0:
                return Constraint.Skip
            return 0 == (
                (
                    blk.config.property_package.charge["Li"]
                    * blk.molar_ion_flux[0, x, "Li"]
                )
                + (
                    blk.config.property_package.charge["Co"]
                    * blk.molar_ion_flux[0, x, "Co"]
                )
                + (
                    blk.config.property_package.charge["Cl"]
                    * blk.molar_ion_flux[0, x, "Cl"]
                )
            )

        self.chloride_flux_membrane = Constraint(
            self.dimensionless_module_length, rule=_chloride_flux_membrane
        )

        # other physical constraints
        def _osmotic_pressure_calculation(blk, x):
            if x == 0:
                return Constraint.Skip
            return blk.osmotic_pressure[0, x] == units.convert(
                (
                    Constants.gas_constant  # J / mol / K
                    * blk.temperature
                    * (
                        (
                            blk.config.property_package.num_solutes["Li"]
                            * blk.config.property_package.sigma["Li"]
                            * (
                                blk.retentate_conc_mol_comp[0, x, "Li"]
                                - blk.permeate_conc_mol_comp[0, x, "Li"]
                            )
                        )
                        + (
                            blk.config.property_package.num_solutes["Co"]
                            * blk.config.property_package.sigma["Co"]
                            * (
                                blk.retentate_conc_mol_comp[0, x, "Co"]
                                - blk.permeate_conc_mol_comp[0, x, "Co"]
                            )
                        )
                        + (
                            blk.config.property_package.num_solutes["Cl"]
                            * blk.config.property_package.sigma["Cl"]
                            * (
                                blk.retentate_conc_mol_comp[0, x, "Cl"]
                                - blk.permeate_conc_mol_comp[0, x, "Cl"]
                            )
                        )
                    )
                ),
                to_units=units.bar,
            )

        self.osmotic_pressure_calculation = Constraint(
            self.dimensionless_module_length, rule=_osmotic_pressure_calculation
        )

        def _electroneutrality_retentate(blk, x):
            return 0 == (
                blk.config.property_package.charge["Li"]
                * blk.retentate_conc_mol_comp[0, x, "Li"]
                + blk.config.property_package.charge["Co"]
                * blk.retentate_conc_mol_comp[0, x, "Co"]
                + blk.config.property_package.charge["Cl"]
                * blk.retentate_conc_mol_comp[0, x, "Cl"]
            )

        self.electroneutrality_retentate = Constraint(
            self.dimensionless_module_length, rule=_electroneutrality_retentate
        )

        def _electroneutrality_membrane(blk, x, z):
            if x == 0:
                return Constraint.Skip
            return 0 == (
                blk.config.property_package.charge["Li"]
                * blk.membrane_conc_mol_comp[0, x, z, "Li"]
                + blk.config.property_package.charge["Co"]
                * blk.membrane_conc_mol_comp[0, x, z, "Co"]
                + blk.config.property_package.charge["Cl"]
                * blk.membrane_conc_mol_comp[0, x, z, "Cl"]
                + blk.membrane_fixed_charge
            )

        self.electroneutrality_membrane = Constraint(
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            rule=_electroneutrality_membrane,
        )

        def _electroneutrality_permeate(blk, x):
            if x == 0:
                return Constraint.Skip
            return 0 == (
                blk.config.property_package.charge["Li"]
                * blk.permeate_conc_mol_comp[0, x, "Li"]
                + blk.config.property_package.charge["Co"]
                * blk.permeate_conc_mol_comp[0, x, "Co"]
                + blk.config.property_package.charge["Cl"]
                * blk.permeate_conc_mol_comp[0, x, "Cl"]
            )

        self.electroneutrality_permeate = Constraint(
            self.dimensionless_module_length, rule=_electroneutrality_permeate
        )

        # partitioning equations
        def _retentate_membrane_interface_lithium(blk, x):
            if x == 0:
                return Constraint.Skip
            return (
                (
                    blk.config.property_package.partition_coefficient_retentate["Li"]
                    ** (-blk.config.property_package.charge["Cl"])
                )
                * (
                    blk.config.property_package.partition_coefficient_retentate["Cl"]
                    ** blk.config.property_package.charge["Li"]
                )
                * (
                    blk.retentate_conc_mol_comp[0, x, "Li"]
                    ** (-blk.config.property_package.charge["Cl"])
                )
                * (
                    blk.retentate_conc_mol_comp[0, x, "Cl"]
                    ** blk.config.property_package.charge["Li"]
                )
            ) == (
                (
                    blk.membrane_conc_mol_comp[0, x, 0, "Li"]
                    ** (-blk.config.property_package.charge["Cl"])
                )
                * (
                    blk.membrane_conc_mol_comp[0, x, 0, "Cl"]
                    ** blk.config.property_package.charge["Li"]
                )
            )

        self.retentate_membrane_interface_lithium = Constraint(
            self.dimensionless_module_length, rule=_retentate_membrane_interface_lithium
        )

        def _retentate_membrane_interface_cobalt(blk, x):
            if x == 0:
                return Constraint.Skip
            return (
                (
                    blk.config.property_package.partition_coefficient_retentate["Co"]
                    ** (-blk.config.property_package.charge["Cl"])
                )
                * (
                    blk.config.property_package.partition_coefficient_retentate["Cl"]
                    ** blk.config.property_package.charge["Co"]
                )
                * (
                    blk.retentate_conc_mol_comp[0, x, "Co"]
                    ** (-blk.config.property_package.charge["Cl"])
                )
                * (
                    blk.retentate_conc_mol_comp[0, x, "Cl"]
                    ** blk.config.property_package.charge["Co"]
                )
            ) == (
                (
                    blk.membrane_conc_mol_comp[0, x, 0, "Co"]
                    ** (-blk.config.property_package.charge["Cl"])
                )
                * (
                    blk.membrane_conc_mol_comp[0, x, 0, "Cl"]
                    ** blk.config.property_package.charge["Co"]
                )
            )

        self.retentate_membrane_interface_cobalt = Constraint(
            self.dimensionless_module_length, rule=_retentate_membrane_interface_cobalt
        )

        def _membrane_permeate_interface_lithium(blk, x):
            if x == 0:
                return Constraint.Skip
            return (
                (
                    blk.config.property_package.partition_coefficient_permeate["Li"]
                    ** (-blk.config.property_package.charge["Cl"])
                )
                * (
                    blk.config.property_package.partition_coefficient_permeate["Cl"]
                    ** blk.config.property_package.charge["Li"]
                )
                * (
                    blk.permeate_conc_mol_comp[0, x, "Li"]
                    ** (-blk.config.property_package.charge["Cl"])
                )
                * (
                    blk.permeate_conc_mol_comp[0, x, "Cl"]
                    ** blk.config.property_package.charge["Li"]
                )
            ) == (
                (
                    blk.membrane_conc_mol_comp[0, x, 1, "Li"]
                    ** (-blk.config.property_package.charge["Cl"])
                )
                * (
                    blk.membrane_conc_mol_comp[0, x, 1, "Cl"]
                    ** blk.config.property_package.charge["Li"]
                )
            )

        self.membrane_permeate_interface_lithium = Constraint(
            self.dimensionless_module_length, rule=_membrane_permeate_interface_lithium
        )

        def _membrane_permeate_interface_cobalt(blk, x):
            if x == 0:
                return Constraint.Skip
            return (
                (
                    blk.config.property_package.partition_coefficient_permeate["Co"]
                    ** (-blk.config.property_package.charge["Cl"])
                )
                * (
                    blk.config.property_package.partition_coefficient_permeate["Cl"]
                    ** blk.config.property_package.charge["Co"]
                )
                * (
                    blk.permeate_conc_mol_comp[0, x, "Co"]
                    ** (-blk.config.property_package.charge["Cl"])
                )
                * (
                    blk.permeate_conc_mol_comp[0, x, "Cl"]
                    ** blk.config.property_package.charge["Co"]
                )
            ) == (
                (
                    blk.membrane_conc_mol_comp[0, x, 1, "Co"]
                    ** (-blk.config.property_package.charge["Cl"])
                )
                * (
                    blk.membrane_conc_mol_comp[0, x, 1, "Cl"]
                    ** blk.config.property_package.charge["Co"]
                )
            )

        self.membrane_permeate_interface_cobalt = Constraint(
            self.dimensionless_module_length, rule=_membrane_permeate_interface_cobalt
        )

        # boundary conditions
        def _retentate_flow_volume_boundary_condition(blk):
            return (
                blk.retentate_flow_volume[0, 0]
                == blk.feed_flow_volume[0] + blk.diafiltrate_flow_volume[0]
            )

        self.retentate_flow_volume_boundary_condition = Constraint(
            rule=_retentate_flow_volume_boundary_condition
        )

        def _retentate_conc_mol_lithium_boundary_condition(blk):
            return blk.retentate_conc_mol_comp[0, 0, "Li"] == (
                (
                    blk.feed_flow_volume[0] * blk.feed_conc_mol_comp[0, "Li"]
                    + blk.diafiltrate_flow_volume[0]
                    * blk.diafiltrate_conc_mol_comp[0, "Li"]
                )
                / (blk.feed_flow_volume[0] + blk.diafiltrate_flow_volume[0])
            )

        self.retentate_conc_mol_lithium_boundary_condition = Constraint(
            rule=_retentate_conc_mol_lithium_boundary_condition
        )

        def _retentate_conc_mol_cobalt_boundary_condition(blk):
            return blk.retentate_conc_mol_comp[0, 0, "Co"] == (
                (
                    blk.feed_flow_volume[0] * blk.feed_conc_mol_comp[0, "Co"]
                    + blk.diafiltrate_flow_volume[0]
                    * blk.diafiltrate_conc_mol_comp[0, "Co"]
                )
                / (blk.feed_flow_volume[0] + blk.diafiltrate_flow_volume[0])
            )

        self.retentate_conc_mol_cobalt_boundary_condition = Constraint(
            rule=_retentate_conc_mol_cobalt_boundary_condition
        )

    def discretize_model(self):
        discretizer = TransformationFactory("dae.finite_difference")
        discretizer.apply_to(
            self,
            wrt=self.dimensionless_module_length,
            nfe=self.config.NFE_module_length,
            scheme="BACKWARD",
        )
        discretizer.apply_to(
            self,
            wrt=self.dimensionless_membrane_thickness,
            nfe=self.config.NFE_membrane_thickness,
            scheme="BACKWARD",
        )

    def fix_boundary_values(self):
        """
        Fix boundary values for the two salt diafiltration unit model.
        """
        for x in self.dimensionless_module_length:
            # chloride concentration gradient in retentate variable is created by default but
            # is not needed in model; fix to reduce number of variables
            self.d_retentate_conc_mol_comp_dx[0, x, "Cl"].fix(
                value(self.numerical_zero_tolerance)
            )
            # associated discretization equation not needed in model
            if x != 0:
                self.d_retentate_conc_mol_comp_dx_disc_eq[0, x, "Cl"].deactivate()

        # set "zero" boundary values to a sufficiently small value
        # these quantities are known to be zero
        # setting them to (numerical) zero improves numerical stability
        self.permeate_flow_volume[0, 0].fix(value(self.numerical_zero_tolerance))
        self.permeate_conc_mol_comp[0, 0, "Li"].fix(
            value(self.numerical_zero_tolerance)
        )
        self.permeate_conc_mol_comp[0, 0, "Co"].fix(
            value(self.numerical_zero_tolerance)
        )
        self.permeate_conc_mol_comp[0, 0, "Cl"].fix(
            value(self.numerical_zero_tolerance)
        )
        for z in self.dimensionless_membrane_thickness:
            self.membrane_conc_mol_comp[0, 0, z, "Li"].fix(
                value(self.numerical_zero_tolerance)
            )
            self.membrane_conc_mol_comp[0, 0, z, "Co"].fix(
                value(self.numerical_zero_tolerance)
            )
            self.membrane_conc_mol_comp[0, 0, z, "Cl"].fix(
                value(self.numerical_zero_tolerance)
            )
        self.d_retentate_conc_mol_comp_dx[0, 0, "Li"].fix(
            value(self.numerical_zero_tolerance)
        )
        self.d_retentate_conc_mol_comp_dx[0, 0, "Co"].fix(
            value(self.numerical_zero_tolerance)
        )
        self.d_retentate_flow_volume_dx[0, 0].fix(value(self.numerical_zero_tolerance))
        self.volume_flux_water[0, 0].fix(value(self.numerical_zero_tolerance))
        self.molar_ion_flux[0, 0, "Li"].fix(value(self.numerical_zero_tolerance))
        self.molar_ion_flux[0, 0, "Co"].fix(value(self.numerical_zero_tolerance))
        self.molar_ion_flux[0, 0, "Cl"].fix(value(self.numerical_zero_tolerance))

    def add_scaling_factors(self):
        """
        Assigns scaling factors to certain variables and constraints to
        improve solver performance.
        """
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        for x in self.dimensionless_module_length:
            if x != 0:
                self.scaling_factor[self.retentate_membrane_interface_cobalt[x]] = 1e-5
                self.scaling_factor[self.retentate_membrane_interface_lithium[x]] = 1e-3
                self.scaling_factor[self.membrane_permeate_interface_cobalt[x]] = 1e-5
                self.scaling_factor[self.membrane_permeate_interface_lithium[x]] = 1e-3

    def add_ports(self):
        self.feed_inlet = Port(doc="Feed Inlet Port")
        self._feed_flow_volume_ref = Reference(self.feed_flow_volume)
        self.feed_inlet.add(self._feed_flow_volume_ref, "flow_vol")
        self._feed_conc_mol_comp_ref = Reference(self.feed_conc_mol_comp)
        self.feed_inlet.add(self._feed_conc_mol_comp_ref, "conc_mol_comp")

        self.diafiltrate_inlet = Port(doc="Diafiltrate Inlet Port")
        self._diafiltrate_flow_volume_ref = Reference(self.diafiltrate_flow_volume)
        self.diafiltrate_inlet.add(self._diafiltrate_flow_volume_ref, "flow_vol")
        self._diafiltrate_conc_mol_comp_ref = Reference(self.diafiltrate_conc_mol_comp)
        self.diafiltrate_inlet.add(self._diafiltrate_conc_mol_comp_ref, "conc_mol_comp")

        self.retentate_outlet = Port(doc="Retentate Outlet Port")
        self._retentate_flow_volume_ref = Reference(
            self.retentate_flow_volume[:, self.dimensionless_module_length.last()]
        )
        self.retentate_outlet.add(self._retentate_flow_volume_ref, "flow_vol")
        self._retentate_conc_mol_comp_ref = Reference(
            self.retentate_conc_mol_comp[:, self.dimensionless_module_length.last(), :]
        )
        self.retentate_outlet.add(self._retentate_conc_mol_comp_ref, "conc_mol_comp")

        self.permeate_outlet = Port(doc="Permeate Outlet Port")
        self._permeate_flow_volume_ref = Reference(
            self.permeate_flow_volume[:, self.dimensionless_module_length.last()]
        )
        self.permeate_outlet.add(self._permeate_flow_volume_ref, "flow_vol")
        self._permeate_conc_mol_comp_ref = Reference(
            self.permeate_conc_mol_comp[:, self.dimensionless_module_length.last(), :]
        )
        self.permeate_outlet.add(self._permeate_conc_mol_comp_ref, "conc_mol_comp")
