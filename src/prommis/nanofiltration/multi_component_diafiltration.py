#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
r"""
Multi-Component Diafiltration Unit Model
========================================

Author: Molly Dougher

This membrane unit model is for the multi-component diafiltration of a multi-salt system with a common anion. The model can be built with or without the assumption of a boundary layer. Currently, the model and property packages support one, two, and three salt systems; however, the model can be extended to :math:`n` salts by supplying the appropriate properties and arguments (see below). The membrane is designed for use in a diafiltration cascade, i.e., the model represents one spiral-wound membrane module piece within a cascade of several membranes.

Configuration Arguments
-----------------------

The Multi-Component Diafiltration unit model requires a property package that provides the valency (:math:`z_i`), diffusion coefficients (:math:`D_{i,bl}` and :math:`D_{i,m}`) within the boundary layer and membrane, respectively, in :math:`\mathrm{mm}^2 \, \mathrm{h}^{-1}`, thermodynamic reflection coefficient (:math:`\sigma_i`), partition coefficients (:math:`H_{i,r}` and :math:`H_{i,p}`) at the retentate-membrane and membrane-permeate interfaces, and number of dissolved species (:math:`n_i`) for each ion :math:`i` in solution. When used in a flowsheet, the user can provide separate property packages for the feed and product streams.

There are six configuration arguments to create an instance of the Multi-Component Diafiltration Unit Model:

#. ``cation_list``: list of cations present in the system

    ``default=["Li", "Co"]``

#. ``anion_list``: list of anions present in the system

    ``default=["Cl"]``

#. ``include_boundary_layer``: Boolean to specify if the model is to be built with a boundary layer

    ``default=True``

#. ``NFE_module_length``: the desired number of finite elements across the width of the membrane (i.e., the module length)

    ``default=10``

#. ``NFE_boundary_layer_thickness``: the desired number of finite elements across the thickness of the boundary layer

    ``default=5``

#. ``NFE_membrane_thickness``: the desired number of finite elements across the thickness of the membrane

    ``default=5``

Degrees of Freedom
------------------

The Multi-Component Diafiltration unit model has :math:`5+2n` degrees of freedom, where :math:`n` is the number of cations in the system:

#. the length of the membrane module (``total_module_length``)
#. the length of the membrane (``total_membrane_length``)
#. the pressure applied to the membrane system (``applied_pressure``)
#. the volumetric flow rate of the feed (``feed_flow_volume``)
#. the cation concentration in the feed (``feed_conc_mol_comp[t,k]``)
#. the volumetric flow rate of the diafiltrate (``diafiltrate_flow_volume``)
#. the cation concentration in the diafiltrate (``diafiltrate_conc_mol_comp[t,k]``)

Model Structure
---------------

There are (up to) four regions in the Multi-Component Diafiltration model: the retentate, the boundary layer, the membrane, and the permeate. The retentate and the permeate are only discretized with respect to module length (:math:`x`-direction), while the boundary layer and membrane are discretized with respect to both module length (:math:`x`-direction) and thickness (:math:`z_{bl}`-direction and :math:`z_{m}`-direction, respectively). The resulting system of partial differential algebraic equations is solved by discretizing with the backward finite difference method.

A schematic of the Multi-Component Diafiltration model's geometry can be found `here <https://github.com/prommis/prommis/blob/main/src/prommis/nanofiltration/membrane_schematic.png>`_.

Assumptions
-----------

* The membrane module dimensions, maximum applied pressure, and inlet flow rates assume that one tube (one instance of this model) consists of 4 NF270-440 membranes in series.

* The partitioning relationships, which describe how the solutes transition (partition) across the solution-membrane interfaces, are derived assuming Donnan equilibrium. The partitioning coefficients incorporate both steric and Donnan effects.

* The default value for the membrane's surface charge (:math:`-44 \, \mathrm{mM}`), was calculated using zeta potential measurements for NF270 membranes. (See `this reference <https://doi.org/10.1021/acs.iecr.4c04763>`_). Currently, the default property package only supports negatively charged membranes.

* The boundary layer thickness is assumed to be :math:`20 \, \mathrm{\mu m}` and the membrane thickness is assumed to be :math:`100 \, \mathrm{nm}`.

* The default value for the membrane permeability (:math:`0.01 \, \mathrm{m} \, \mathrm{h}^{-1} \, \mathrm{bar}^{-1}`) is based off of parameter estimation results from `this reference <https://doi.org/10.1021/acs.iecr.4c04763>`_ for NF270 membranes.

* The dominating transport mechanism within the bulk/retentate and permeate solutions is convection.

* The transport mechanisms modeled within the both the boundary layer and the membrane are convection, diffusion, and electromigration. 

* The system is uniform with respect to the wound-dimension of the membrane.

* Diffusion that occurs normal to the direction of flux is assumed to be negligible, meaning the concentration gradient of ion :math:`i` only has a :math:`z_{bl}`- or :math:`z_m`-component (perpendicular to the membrane surface).

Sets
----

The Multi-Component Diafiltration model defines the following discrete sets for solutes (:math:`\mathcal{I}`) and cations (:math:`\mathcal{K}`) in the system:

.. math:: \mathcal{I}=\{\mathrm{cation_1, cation_2, ..., cation_n, anion}\}
.. math:: \mathcal{K}=\{\mathrm{cation_1, cation_2, ..., cation_n}\}

where :math:`n` is the desired number of cations.

There are 3 continuous sets for each length dimension: ``dimensionless_module_length`` (in the :math:`x`-direction parallel to the membrane surface), ``dimensionless_boundary_layer_thickness`` (in the :math:`z_{bl}`-direction perpendicular to the membrane surface), and ``dimensionless_membrane_thickness`` (in the :math:`z_m`-direction perpendicular to the membrane surface). The length dimensions, :math:`x`, :math:`z_{bl}`, and :math:`z_m`, are non-dimensionalized as :math:`\bar{x}`, :math:`\bar{z}_{bl}`, and :math:`\bar{z}_m`, respectively, using the module length (:math:`w`), boundary layer thickness (:math:`\delta`), and membrane thickness (:math:`l`), respectively, to improve numerical stability.

.. math:: \bar{x} \in \mathbb{R} \| 0 \leq \bar{x} \leq 1
.. math:: \bar{z}_{bl} \in \mathbb{R} \| 0 \leq \bar{z}_{bl} \leq 1
.. math:: \bar{z}_m \in \mathbb{R} \| 0 \leq \bar{z}_m \leq 1

*Note:* :math:`z_{bl}` *and* :math:`z_m` *point in the same direction (perpendicular to the membrane surface), but are defined as separate length scales to simplify the implementation of the model. The appropriate boundary conditions between* :math:`z_{bl}` *and* :math:`z_m` *are enforced within the model.*

Though this is not implemented as a dynamic model, a set is defined for time.

.. math:: t \in [0]

Default Model Parameters
------------------------

The Multi-Component Diafiltration model has the following parameters.

================ =============================================== ================================== ============= ==========================================================
Parameter        Description                                     Name                               Default Value Units
================ =============================================== ================================== ============= ==========================================================
:math:`\epsilon` numerical tolerance for zero values             ``numerical_zero_tolerance``       :math:`1e-10`
:math:`\delta`   boundary layer thickness                        ``total_boundary_layer_thickness`` :math:`2e-05` :math:`\mathrm{m}`
:math:`l`        membrane thickness                              ``total_membrane_thickness``       :math:`1e-07` :math:`\mathrm{m}`
:math:`L_p`      hydraulic permeability of the membrane          ``membrane_permeability``          :math:`0.01`  :math:`\mathrm{m} \, \mathrm{h}^{-1} \, \mathrm{bar}^{-1}`
:math:`T`        temperature of the system                       ``temperature``                    :math:`298`   :math:`\mathrm{K}`
:math:`\chi`     concentration of surface charge on the membrane ``membrane_fixed_charge``          :math:`-44`   :math:`\mathrm{mol} \, \mathrm{m}^{-3}`
================ =============================================== ================================== ============= ==========================================================

Variables
---------

The Multi-Component Diafiltration model adds the following variables.

=============================== ==================================================================== ======================================================= =========================================================================== ==========================================================================================================
Variable                        Description                                                          Name                                                    Units                                                                       Indexed over
=============================== ==================================================================== ======================================================= =========================================================================== ==========================================================================================================
:math:`c_{i,bl}`                ion concentration in the boundary layer                              ``boundary_layer_conc_mol_comp``                        :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                                     :math:`t`, :math:`\bar{x}`, :math:`\bar{z}_{bl}`, and :math:`i \in \mathcal{I}`
:math:`c_{i,d}`                 ion concentration in the diafiltrate                                 ``diafiltrate_conc_mol_comp``                           :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                                     :math:`t` and :math:`i \in \mathcal{I}`
:math:`c_{i,f}`                 ion concentration in the feed                                        ``feed_conc_mol_comp``                                  :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                                     :math:`t` and :math:`i \in \mathcal{I}`
:math:`c_{i,m}`                 ion concentration in the membrane                                    ``membrane_conc_mol_comp``                              :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                                     :math:`t`, :math:`\bar{x}`, :math:`\bar{z}_m`, and :math:`i \in \mathcal{I}`
:math:`c_{i,p}`                 ion concentration in the permeate                                    ``permeate_conc_mol_comp``                              :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                                     :math:`t`, :math:`\bar{x}`, and :math:`i \in \mathcal{I}`
:math:`c_{i,r}`                 ion concentration in the retentate                                   ``retentate_conc_mol_comp``                             :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                                     :math:`t`, :math:`\bar{x}`, and :math:`i \in \mathcal{I}`
:math:`\tilde{D}_{bl}`          cross-diffusion coefficient denominator in the boundary layer        ``boundary_layer_D_tilde``                              :math:`\mathrm{mm}^2 \, \mathrm{h}^{-1} \, \mathrm{mol} \, \mathrm{m}^{-3}` :math:`t`, :math:`\bar{x}`, and :math:`\bar{z}_{bl}`
:math:`D_{kj,bl}^{bilinear}`    bilinear cross-diffusion coefficient in the membrane                 ``boundary_layer_cross_diffusion_coefficient_bilinear`` :math:`\mathrm{mm}^4 \, \mathrm{h}^{-2} \, \mathrm{mol} \, \mathrm{m}^{-3}` :math:`t`, :math:`\bar{x}`, :math:`\bar{z}_{bl}`, :math:`k \in \mathcal{K}`, and :math:`j \in \mathcal{K}`
:math:`D_{kj,bl}`               cross-diffusion coefficient in the membrane                          ``boundary_layer_cross_diffusion_coefficient``          :math:`\mathrm{mm}^2 \, \mathrm{h}^{-1}`                                    :math:`t`, :math:`\bar{x}`, :math:`\bar{z}_{bl}`, :math:`k \in \mathcal{K}`, and :math:`j \in \mathcal{K}`
:math:`\tilde{D}_m`             cross-diffusion & convection coefficient denominator in the membrane ``membrane_D_tilde``                                    :math:`\mathrm{mm}^2 \, \mathrm{h}^{-1} \, \mathrm{mol} \, \mathrm{m}^{-3}` :math:`t`, :math:`\bar{x}`, and :math:`\bar{z}_m`
:math:`D_{kj,m}^{bilinear}`     bilinear cross-diffusion coefficient in the membrane                 ``membrane_cross_diffusion_coefficient_bilinear``       :math:`\mathrm{mm}^4 \, \mathrm{h}^{-2} \, \mathrm{mol} \, \mathrm{m}^{-3}` :math:`t`, :math:`\bar{x}`, :math:`\bar{z}_m`, :math:`k \in \mathcal{K}`, and :math:`j \in \mathcal{K}`
:math:`\alpha_{k,m}^{bilinear}` bilinear convection coefficient in the membrane                      ``membrane_convection_coefficient_bilinear``            :math:`\mathrm{mm}^2 \, \mathrm{h}^{-1} \, \mathrm{mol} \, \mathrm{m}^{-3}` :math:`t`, :math:`\bar{x}`, :math:`\bar{z}_m`, and :math:`k \in \mathcal{K}`
:math:`D_{kj,m}`                cross-diffusion coefficient in the membrane                          ``membrane_cross_diffusion_coefficient``                :math:`\mathrm{mm}^2 \, \mathrm{h}^{-1}`                                    :math:`t`, :math:`\bar{x}`, :math:`\bar{z}_m`, :math:`k \in \mathcal{K}`, and :math:`j \in \mathcal{K}`
:math:`\alpha_{k,m}`            convection coefficient in the membrane                               ``membrnane_convection_coefficient``                    :math:`\mathrm{dimensionless}`                                              :math:`t`, :math:`\bar{x}`, :math:`\bar{z}_m`, and :math:`k \in \mathcal{K}`
:math:`j_i`                     molar flux of ions through the membrane                              ``molar_ion_flux``                                      :math:`\mathrm{mol} \, \mathrm{m}^{-2} \, \mathrm{h}^{-1}`                  :math:`t`, :math:`\bar{x}`, and :math:`i \in \mathcal{I}`
:math:`J_w`                     water flux across the membrane                                       ``volume_flux_water``                                   :math:`\mathrm{m}^3 \, \mathrm{m}^{-2} \, \mathrm{h}^{-1}`                  :math:`t` and :math:`\bar{x}`
:math:`L`                       length of the membrane                                               ``total_membrane_length``                               :math:`\mathrm{m}`
:math:`\Delta \pi`              osmotic pressure of feed-side fluid                                  ``osmotic_pressure``                                    :math:`\mathrm{bar}`                                                        :math:`t` and :math:`\bar{x}`
:math:`\Delta P`                applied pressure to the membrane                                     ``applied_pressure``                                    :math:`\mathrm{bar}`                                                        :math:`t`
:math:`q_d`                     volumetric flow rate of the diafiltrate                              ``diafiltrate_flow_volume``                             :math:`\mathrm{m}^3 \, \mathrm{h}^{-1}`                                     :math:`t`
:math:`q_f`                     volumetric flow rate of the feed                                     ``feed_flow_volume``                                    :math:`\mathrm{m}^3 \, \mathrm{h}^{-1}`                                     :math:`t`
:math:`q_p`                     volumetric flow rate of the permeate                                 ``permeate_flow_volume``                                :math:`\mathrm{m}^3 \, \mathrm{h}^{-1}`                                     :math:`t` and :math:`\bar{x}`
:math:`q_r`                     volumetric flow rate of the retentate                                ``retentate_flow_volume``                               :math:`\mathrm{m}^3 \, \mathrm{h}^{-1}`                                     :math:`t` and :math:`\bar{x}`
:math:`w`                       length of the membrane module                                        ``total_module_length``                                 :math:`\mathrm{m}`
=============================== ==================================================================== ======================================================= =========================================================================== ==========================================================================================================

Derivative Variables
--------------------

The Multi-Component Diafiltration model adds the following derivative variables.

======================================================= ================================================ ===================================== ======================================= ===============================================================================
Variable                                                Description                                      Name                                  Units                                   Indexed over
======================================================= ================================================ ===================================== ======================================= ===============================================================================
:math:`\frac{\mathrm{d}c_{k,r}}{\mathrm{d}\bar{x}}`     ion concentration gradient in the retentate      ``d_retentate_conc_mol_comp_dx``      :math:`\mathrm{mol} \, \mathrm{m}^{-3}` :math:`t`, :math:`\bar{x}`, and :math:`k \in \mathcal{K}`
:math:`\frac{\mathrm{d}q_r}{\mathrm{d}\bar{x}}`         retentate flow rate gradient                     ``d_retentate_flow_volume_dx``        :math:`\mathrm{m}^3 \, \mathrm{h}^{-1}` :math:`t` and :math:`\bar{x}`
:math:`\frac{\partial c_{k,bl}}{\partial \bar{z}_{bl}}` ion concentration gradient in the boundary layer ``d_boundary_layer_conc_mol_comp_dz`` :math:`\mathrm{mol} \, \mathrm{m}^{-3}` :math:`t`, :math:`\bar{x}`, :math:`\bar{z}_{bl}`, and :math:`k \in \mathcal{K}`
:math:`\frac{\partial c_{k,m}}{\partial \bar{z}_m}`     ion concentration gradient in the membrane       ``d_membrane_conc_mol_comp_dz``       :math:`\mathrm{mol} \, \mathrm{m}^{-3}` :math:`t`, :math:`\bar{x}`, :math:`\bar{z}_m`, and :math:`k \in \mathcal{K}`
======================================================= ================================================ ===================================== ======================================= ===============================================================================

Constraints
-----------

**Differential mole balances:**

.. math:: \frac{\mathrm{d}q_r(\bar{x})}{\mathrm{d}\bar{x}} = - J_w(\bar{x}) wL  \qquad \forall \, \bar{x} \in (0, 1]
.. math:: q_r(\bar{x}) \frac{\mathrm{d}c_{k,r}(\bar{x})}{\mathrm{d}\bar{x}} = wL (J_w(\bar{x}) c_{k,r}(\bar{x}) - j_{k}(\bar{x}))  \qquad \forall \, \bar{x} \in (0, 1], \, k \in \mathcal{K}

**Bulk flux balances:**

.. math:: q_p(\bar{x}) = \bar{x} wL J_w(\bar{x}) \qquad \forall \, \bar{x} \in (0, 1]
.. math:: j_{k}(\bar{x}) = c_{k,p}(\bar{x}) J_w(\bar{x}) \qquad \forall \, \bar{x} \in (0, 1], \, k \in \mathcal{K}

**Overall water flux through the membrane:**

.. math:: J_w (\bar{x}) = L_p (\Delta P - \Delta \pi (\bar{x})) \qquad \forall \, \bar{x} \in (0, 1]

*without a boundary layer:*

.. math:: \Delta \pi (\bar{x}) = \mathrm{R} \mathrm{T} \sum_{i \in \mathcal{I}} n_i \sigma_i (c_{i,r}(\bar{x})-c_{i,p}(\bar{x})) \qquad \forall \, \bar{x} \in (0, 1]

*with a boundary layer:*

.. math:: \Delta \pi (\bar{x}) = \mathrm{R} \mathrm{T} \sum_{i \in \mathcal{I}} n_i \sigma_i (c_{i,bl}(\bar{x}, \bar{z}_{bl}=1)-c_{i,p}(\bar{x})) \qquad \forall \, \bar{x} \in (0, 1]

**Cation flux through the boundary layer and membrane:**

*Derived from the extended Nernst-Planck equation*

.. math:: j_k(\bar{x}) = c_{k,bl}(\bar{x},\bar{z}_{bl}) J_w(\bar{x}) + \frac{1}{\delta} \sum_{j \in \mathcal{K}} \left(D_{kj,bl} (\bar{x},\bar{z}_{bl}) \nabla c_{j,bl} (\bar{x},\bar{z}_{bl}) \right) \qquad \forall \, \bar{x} \in (0, 1], \, k \in \mathcal{K}
.. math:: j_k(\bar{x}) = \alpha_{k,m}(\bar{x},\bar{z}_m) c_{k,m}(\bar{x},\bar{z}_m) J_w(\bar{x}) + \frac{1}{l} \sum_{j \in \mathcal{K}} \left(D_{kj,m} (\bar{x},\bar{z}_m) \nabla c_{j,m} (\bar{x},\bar{z}_m) \right) \qquad \forall \, \bar{x} \in (0, 1], \, k \in \mathcal{K}

where

.. math:: 
    \alpha_{k,h}(\bar{x},\bar{z}_h) = 
    \begin{cases}
        1,& \text{if } h = bl \\
        1 + \dfrac{z_k D_{k,h} \chi}{\tilde{D}_h (\bar{x},\bar{z}_h)},& \text{if } h = m
    \end{cases}
.. math:: 
    D_{kj,h}(\bar{x},\bar{z}_h) = 
    \begin{cases}
        \dfrac{(z_k z_j D_{k,h} D_{j,h} - z_k z_j D_{k,h} D_{a,h})c_{k,h} (\bar{x},\bar{z}_h)}{\tilde{D}_h (\bar{x},\bar{z}_h)},& \text{if } k \neq j, h \in \{bl, m\} \\
        \dfrac{\sum_{t \in \mathcal{C}} \left((z_t z_a D_{k,h} D_{a,h} - \beta_{kt,h})c_{t,h} (\bar{x},\bar{z}_h) \right)}{\tilde{D}_h (\bar{x},\bar{z}_h)} ,& \text{if } k=j, h=bl \\
        \dfrac{\sum_{t \in \mathcal{C}} \left((z_t z_a D_{k,h} D_{a,h} - \beta_{kt,h})c_{t,h} (\bar{x},\bar{z}_h) \right) + z_a D_{k,h} D_{a,h} \chi}{\tilde{D}_h (\bar{x},\bar{z}_h)} ,& \text{if } k=j, h=m \\
    \end{cases}
.. math::
    \beta_{kt,h} = 
    \begin{cases}
        z_t^2 D_{t,h} D_{k,h} ,& \text{if } k\neq t \\
        z_t^2 D_{t,h} D_{a,h} ,& \text{if } k=t \\
    \end{cases}
.. math:: 
    \tilde{D}_h (\bar{x},\bar{z}_h) = 
    \begin{cases}
        \sum_{j \in \mathcal{K}} \left((z_j^2 D_{j,h} - z_j z_a D_{a,h})c_{j,h} (\bar{x},\bar{z}_h) \right),& \text{if } h = bl \\
        \sum_{j \in \mathcal{K}} \left((z_j^2 D_{j,h} - z_j z_a D_{a,h})c_{j,h} (\bar{x},\bar{z}_h) \right) - z_a D_{a,h} \chi,& \text{if } h = m \\
    \end{cases}
.. math:: \nabla c_{k,h} (\bar{x},\bar{z}_h)= \dfrac{\partial c_{k,h}(\bar{x},\bar{z}_h)}{\partial \bar{z}_h}

and the subscript :math:`a` represents the anion in solution.

The diffusion and convection coefficients are reformulated to bilinear constraints:

.. math:: \alpha_{k,m}^{bilinear}(\bar{x},\bar{z}_m) = \alpha_{k,m}(\bar{x},\bar{z}_m) \tilde{D}_m(\bar{x},\bar{z}_m) = \tilde{D}_m(\bar{x},\bar{z}_m) + z_k D_{k,m} \chi \qquad \forall \, \bar{x} \in (0, 1]
.. math:: D_{kj,h}^{bilinear}(\bar{x},\bar{z}_h) = D_{kj,h}(\bar{x},\bar{z}_h) \tilde{D}_h(\bar{x},\bar{z}_h) \qquad \forall \, h \in \{bl, m\}, \, \bar{x} \in (0, 1]

*Note that the single solute diffusion coefficients are provided in* :math:`\mathrm{mm}^2\ \, \mathrm{h}^{-1}` *to improve numerical stability. When used in the Nernst-Planck equations, the diffusion coefficients are converted to* :math:`\mathrm{m}^2\ \, \mathrm{h}^{-1}`.

**No applied potential on the system:**

.. math:: 0 = \sum_{i \in \mathcal{I}} z_i j_i(\bar{x}) \qquad \forall \, \bar{x} \in (0, 1]

**Electroneutrality:**

.. math:: 0 = \sum_{i \in \mathcal{I}} z_i c_{i,r}(\bar{x})
.. math:: 0 = \sum_{i \in \mathcal{I}} z_i c_{i,bl}(\bar{x},\bar{z}_{bl}) \qquad \forall \, \bar{x} \in (0, 1]
.. math:: 0 = \chi + \sum_{i \in \mathcal{I}} z_i c_{i,m}(\bar{x},\bar{z}_m) \qquad \forall \, \bar{x} \in (0, 1]
.. math:: 0 = \sum_{i \in \mathcal{I}} z_i c_{i,p}(\bar{x}) \qquad \forall \, \bar{x} \in (0, 1]

**Partitioning:**

At the feed-side solution-membrane interface:

*without a boundary layer:*

.. math:: H_{k,r}^{-z_a} H_{a,r}^{z_k} = \left(\frac{c_{k,m} (\bar{x},\bar{z}=0)}{c_{k,r} (\bar{x})}\right)^{-z_a} \left(\frac{c_{a,m} (\bar{x},\bar{z}=0)}{c_{a,r}(\bar{x})}\right)^{z_k} \qquad \forall \, \bar{x} \in (0, 1], \, k \in \mathcal{K}

*with a boundary layer:*

.. math:: c_{k,r} (\bar{x}) = c_{k,bl} (\bar{x},\bar{z}_{bl}=0) \qquad \forall \, \bar{x} \in (0, 1], \, k \in \mathcal{K}
.. math:: H_{k,r}^{-z_a} H_{a,r}^{z_k} = \left(\frac{c_{k,m} (\bar{x},\bar{z}_m=0)} {c_{k,bl} (\bar{x},\bar{z}_{bl}=1)}\right)^{-z_a} \left(\frac{c_{a,m} (\bar{x},\bar{z}_m=0)}{c_{a,bl}(\bar{x},\bar{z}_{bl}=1)}\right)^{z_k} \qquad \forall \, \bar{x} \in (0, 1], \, k \in \mathcal{K}


At the membrane-permeate interface:

.. math:: H_{k,p}^{-z_a} H_{a,p}^{z_k} = \left(\frac{c_{k,m} (\bar{x},\bar{z}=1)}{c_{k,p} (\bar{x})}\right)^{-z_a} \left(\frac{c_{a,m} (\bar{x},\bar{z}=1)}{c_{a,p}(\bar{x})}\right)^{z_k} \qquad \forall \, \bar{x} \in (0, 1], \, k \in \mathcal{K}

**Boundary conditions:**

.. math:: q_r(\bar{x}=0) = q_f + q_d
.. math:: c_{k,r}(\bar{x}=0) = \frac{q_f c_{k,f} + q_d c_{k,d}}{q_f + q_d} \qquad \forall \, k \in \mathcal{K}
.. math:: c_{k,bl} (\bar{x}=0,\bar{z}_{bl}) = \epsilon \qquad \forall \, \bar{z}_{bl}, \, k \in \mathcal{K}
.. math:: c_{k,m} (\bar{x}=0,\bar{z}_m) = \epsilon \qquad \forall \, \bar{z}_m, \, k \in \mathcal{K}

The following constraints (which are expected to be zero) are enforced to improve numerical stability (with the appropriate constraints deactivated as described above):

.. math:: q_p(\bar{x}=0) = \epsilon
.. math:: c_{i,p}(\bar{x}=0) = \epsilon \qquad \forall \, i \in \mathcal{I}
.. math:: \frac{\mathrm{d}q_r(\bar{x})}{\mathrm{d}\bar{x}}(\bar{x}=0)=\epsilon
.. math:: \frac{\mathrm{d}c_{k,r}(\bar{x})}{\mathrm{d}\bar{x}}(\bar{x}=0)=\epsilon \qquad \forall \, k \in \mathcal{K}
.. math:: J_w(\bar{x}=0) = \epsilon
.. math:: j_i(\bar{x}=0) = \epsilon \qquad \forall \, i \in \mathcal{I}
"""

from pyomo.common.config import ConfigBlock, ConfigValue, ListOf
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.environ import (
    Constraint,
    Expression,
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
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core import UnitModelBlockData, declare_process_block_class, useDefault
from idaes.core.initialization import BlockTriangularizationInitializer
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.constants import Constants
from idaes.core.util.exceptions import ConfigurationError


class MultiComponentDiafiltrationInitializer(BlockTriangularizationInitializer):
    """
    Multi-Component Diafiltration Initializer Class.
    """

    def initialization_routine(self, model):
        """
        Initializes the retentate and permeate streams, membrane and boundary
        layer concentrations, and un-initialized derivative variables.

        Method then calls the block triangularization initializer method.
        """

        for t in model.time:
            for x in model.dimensionless_module_length:
                model.retentate_flow_volume[t, x].set_value(
                    value(model.feed_flow_volume[t]) * 1 / 3
                )
                model.d_retentate_flow_volume_dx[t, x].set_value(-10)
                model.permeate_flow_volume[t, x].set_value(
                    value(model.feed_flow_volume[t]) * 2 / 3
                )
                for j in model.solutes:
                    model.retentate_conc_mol_comp[t, x, j].set_value(
                        value(model.feed_conc_mol_comp[t, j]) * 0.95
                    )
                    if len(model.config.cation_list) == 1:
                        model.d_retentate_conc_mol_comp_dx[t, x, j].set_value(1)
                    else:
                        model.d_retentate_conc_mol_comp_dx[t, x, j].set_value(10)
                    model.permeate_conc_mol_comp[t, x, j].set_value(
                        value(model.feed_conc_mol_comp[t, j]) * 0.8
                    )
                if model.config.include_boundary_layer:
                    for z in model.dimensionless_boundary_layer_thickness:
                        for j in model.solutes:
                            model.boundary_layer_conc_mol_comp[t, x, z, j].set_value(
                                value(model.feed_conc_mol_comp[t, j]) * 0.75
                            )
                            model.d_boundary_layer_conc_mol_comp_dz[
                                t, x, z, j
                            ].set_value(10)
                        # update diffusion coefficients
                        if x != 0:
                            calculate_variable_from_constraint(
                                model.boundary_layer_D_tilde[t, x, z],
                                model.boundary_layer_D_tilde_calculation[t, x, z],
                            )
                            for k in model.cations:
                                for j in model.cations:
                                    calculate_variable_from_constraint(
                                        model.boundary_layer_cross_diffusion_coefficient_bilinear[
                                            t, x, z, k, j
                                        ],
                                        model.boundary_layer_cross_diffusion_coefficient_calculation[
                                            t, x, z, k, j
                                        ],
                                    )
                                    calculate_variable_from_constraint(
                                        model.boundary_layer_cross_diffusion_coefficient[
                                            t, x, z, k, j
                                        ],
                                        model.boundary_layer_cross_diffusion_coefficient_bilinear_calculation[
                                            t, x, z, k, j
                                        ],
                                    )

                for z in model.dimensionless_membrane_thickness:
                    for j in model.solutes:
                        # adjust membrane concentration based on charge for 3 salt system
                        if len(model.config.cation_list) == 1:
                            model.membrane_conc_mol_comp[t, x, z, j].set_value(
                                value(model.feed_conc_mol_comp[t, j]) * 0.1
                            )
                        else:
                            if value(model.config.property_package.charge[j]) == 1:
                                model.membrane_conc_mol_comp[t, x, z, j].set_value(
                                    value(model.feed_conc_mol_comp[t, j]) * 0.2
                                )
                            elif value(model.config.property_package.charge[j]) >= 2:
                                model.membrane_conc_mol_comp[t, x, z, j].set_value(
                                    value(model.feed_conc_mol_comp[t, j]) * 1e-2
                                )
                            # update anion concentration to consider fixed membrane charge
                            if x != 0:
                                calculate_variable_from_constraint(
                                    model.membrane_conc_mol_comp[
                                        t, x, z, model.config.anion_list[0]
                                    ],
                                    model.electroneutrality_membrane[t, x, z],
                                )
                        # Note: this threshold is not rigorously tested
                        if value(model.feed_ionic_strength[t]) < 800:
                            model.d_membrane_conc_mol_comp_dz[t, x, z, j].set_value(1)
                        else:
                            model.d_membrane_conc_mol_comp_dz[t, x, z, j].set_value(0.1)

                    # update diffusion and convection coefficients
                    # improves numerics for multi-salt systems
                    if len(model.config.cation_list) >= 3:
                        if x != 0:
                            calculate_variable_from_constraint(
                                model.membrane_D_tilde[t, x, z],
                                model.membrane_D_tilde_calculation[t, x, z],
                            )
                            for k in model.cations:
                                calculate_variable_from_constraint(
                                    model.membrane_convection_coefficient_bilinear[
                                        t, x, z, k
                                    ],
                                    model.membrane_convection_coefficient_calculation[
                                        t, x, z, k
                                    ],
                                )
                                calculate_variable_from_constraint(
                                    model.membrane_convection_coefficient[t, x, z, k],
                                    model.membrane_convection_coefficient_bilinear_calculation[
                                        t, x, z, k
                                    ],
                                )
                                for j in model.cations:
                                    calculate_variable_from_constraint(
                                        model.membrane_cross_diffusion_coefficient_bilinear[
                                            t, x, z, k, j
                                        ],
                                        model.membrane_cross_diffusion_coefficient_calculation[
                                            t, x, z, k, j
                                        ],
                                    )
                                    calculate_variable_from_constraint(
                                        model.membrane_cross_diffusion_coefficient[
                                            t, x, z, k, j
                                        ],
                                        model.membrane_cross_diffusion_coefficient_bilinear_calculation[
                                            t, x, z, k, j
                                        ],
                                    )

        super().initialization_routine(model)


@declare_process_block_class("MultiComponentDiafiltration")
class MultiComponentDiafiltrationData(UnitModelBlockData):
    """
    Multi-Component Diafiltration Unit Model Class.
    """

    # Set default initializer
    default_initializer = MultiComponentDiafiltrationInitializer

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
        "cation_list",
        ConfigValue(
            domain=ListOf(str),
            default=["Li", "Co"],
            doc="List of cations present in the system",
        ),
    )
    CONFIG.declare(
        "anion_list",
        ConfigValue(
            domain=ListOf(str),
            default=["Cl"],
            doc="List of anions present in the system",
        ),
    )
    CONFIG.declare(
        "include_boundary_layer",
        ConfigValue(
            default=True,
            doc="Boolean to specify if the model is to be built with a boundary layer",
        ),
    )
    CONFIG.declare(
        "NFE_module_length",
        ConfigValue(
            default=10,
            doc="Number of finite elements across module length (in the x-direction)",
        ),
    )
    CONFIG.declare(
        "NFE_boundary_layer_thickness",
        ConfigValue(
            default=5,
            doc="Number of finite elements across the boundary layer (in the z-direction)",
        ),
    )
    CONFIG.declare(
        "NFE_membrane_thickness",
        ConfigValue(
            default=5,
            doc="Number of finite elements across the membrane thickness (in the z-direction)",
        ),
    )

    def build(self):
        """
        Build method for the multi-component diafiltration unit model.
        """
        super().build()

        if len(self.config.anion_list) > 1:
            raise ConfigurationError(
                "The multi-component diafiltration unit model only supports systems with a common anion"
            )

        self.add_mutable_parameters()
        self.add_variables()
        self.add_constraints()
        self.discretize_model()
        self.deactivate_unnecessary_objects()
        self.add_scaling_factors()
        self.add_ports()
        self.add_helpful_expressions()

    def add_mutable_parameters(self):
        """
        Adds default parameters for the multi-component diafiltration unit model.

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
        if self.config.include_boundary_layer:
            self.total_boundary_layer_thickness = Param(
                initialize=2e-5,  # Baker, Chapter 4, page 176
                mutable=True,
                units=units.m,
                doc="Thickness of boundary layer (z-direction)",
            )
        self.total_membrane_thickness = Param(
            initialize=1e-7,
            mutable=True,
            units=units.m,
            doc="Thickness of membrane (z-direction)",
        )
        self.membrane_fixed_charge = Param(
            initialize=-44,
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
        Adds variables for the multi-component diafiltration unit model.

        Membrane module dimensions and maximum flowrate (17 m3/h) are
        estimated from NF270-440 modules.

        Assumes 4 modules in series.
        """
        # define length scales
        self.dimensionless_module_length = ContinuousSet(bounds=(0, 1))
        if self.config.include_boundary_layer:
            self.dimensionless_boundary_layer_thickness = ContinuousSet(bounds=(0, 1))
        self.dimensionless_membrane_thickness = ContinuousSet(bounds=(0, 1))

        # add a time index since the property package variables are indexed over time
        self.time = Set(initialize=[0])

        # add components
        self.solutes = Set(initialize=self.config.cation_list + self.config.anion_list)
        self.cations = Set(initialize=self.config.cation_list)

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
            self.time,
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
            vals = {k: 200 for k in self.config.cation_list}
            vals.update({self.config.anion_list[0]: 600})
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
            vals = {k: 10 for k in self.config.cation_list}
            vals.update({self.config.anion_list[0]: 30})
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

        def initialize_molar_ion_flux(m, t, w, j):
            vals = {k: 10 for k in self.config.cation_list}
            vals.update({self.config.anion_list[0]: 30})
            return vals[j]

        self.molar_ion_flux = Var(
            self.time,
            self.dimensionless_module_length,
            self.solutes,
            initialize=initialize_molar_ion_flux,
            units=units.mol / units.m**2 / units.h,
            bounds=[1e-11, None],
            doc="Mole flux of solutes across the membrane (z-direction, x-dependent)",
        )
        self.retentate_flow_volume = Var(
            self.time,
            self.dimensionless_module_length,
            initialize=6.75,
            units=units.m**3 / units.h,
            bounds=[1e-11, None],
            doc="Volumetric flow rate of the retentate, x-dependent",
        )

        def initialize_retentate_conc_mol_comp(m, t, w, j):
            vals = {
                i: 0.95 * initialize_feed_conc_mol_comp(m, t, i) for i in self.solutes
            }
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
            initialize=10,
            units=units.m**3 / units.h,
            bounds=[1e-11, None],
            doc="Volumetric flow rate of the permeate, x-dependent",
        )

        def initialize_permeate_conc_mol_comp(m, t, w, j):
            vals = {
                i: 0.8 * initialize_feed_conc_mol_comp(m, t, i) for i in self.solutes
            }
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
        if self.config.include_boundary_layer:

            def initialize_boundary_layer_conc_mol_comp(m, t, w, l, j):
                vals = {
                    i: 0.5 * initialize_feed_conc_mol_comp(m, t, i)
                    for i in self.solutes
                }
                return vals[j]

            self.boundary_layer_conc_mol_comp = Var(
                self.time,
                self.dimensionless_module_length,
                self.dimensionless_boundary_layer_thickness,
                self.solutes,
                initialize=initialize_boundary_layer_conc_mol_comp,
                units=units.mol / units.m**3,  # mM
                bounds=[1e-11, None],
                doc="Mole concentration of solutes in the boundary layer, x- and z-dependent",
            )

            self.boundary_layer_D_tilde = Var(
                self.time,
                self.dimensionless_module_length,
                self.dimensionless_boundary_layer_thickness,
                initialize=1000,
                units=(units.mm**2 / units.hr) * (units.mol / units.m**3),  # D * c
                doc="Denominator of diffusion and convection coefficients in boundary layer",
            )

            def initialize_boundary_layer_cross_diffusion_coefficient_bilinear(
                m, t, w, l, j, k
            ):
                vals = {
                    k: {j: -3000 for j in self.config.cation_list}
                    for k in self.config.cation_list
                }
                return vals[j][k]

            self.boundary_layer_cross_diffusion_coefficient_bilinear = Var(
                self.time,
                self.dimensionless_module_length,
                self.dimensionless_boundary_layer_thickness,
                self.cations,
                self.cations,
                initialize=initialize_boundary_layer_cross_diffusion_coefficient_bilinear,
                units=(units.mm**2 / units.h)
                * (units.mm**2 / units.h * units.mol / units.m**3),  # D * D,tilde
                doc="Bi-linear cross diffusion coefficient for cations in boundary layer",
            )

            def initialize_boundary_layer_cross_diffusion_coefficient(m, t, w, l, j, k):
                vals = {
                    k: {j: -5 for j in self.config.cation_list}
                    for k in self.config.cation_list
                }
                return vals[j][k]

            self.boundary_layer_cross_diffusion_coefficient = Var(
                self.time,
                self.dimensionless_module_length,
                self.dimensionless_boundary_layer_thickness,
                self.cations,
                self.cations,
                initialize=initialize_boundary_layer_cross_diffusion_coefficient,
                units=units.mm**2 / units.h,
                doc="Cross diffusion coefficient for cations in boundary layer",
            )

        def initialize_membrane_conc_mol_comp(m, t, w, l, j):
            vals = {
                i: 0.1 * initialize_feed_conc_mol_comp(m, t, i) for i in self.solutes
            }
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
            vals = {
                k: {j: -3000 for j in self.config.cation_list}
                for k in self.config.cation_list
            }
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
            doc="Bi-linear cross diffusion coefficient for cations in membrane",
        )

        def initialize_membrane_convection_coefficient_bilinear(m, t, w, l, j):
            vals = {k: 100 for k in self.config.cation_list}
            return vals[j]

        self.membrane_convection_coefficient_bilinear = Var(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            self.cations,
            initialize=initialize_membrane_convection_coefficient_bilinear,
            units=(units.mm**2 / units.hr) * (units.mol / units.m**3),  # D,tilde
            doc="Convection coefficient for cations in membrane",
        )

        def initialize_membrane_cross_diffusion_coefficient(m, t, w, l, j, k):
            vals = {
                k: {j: -5 for j in self.config.cation_list}
                for k in self.config.cation_list
            }
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
            vals = {k: 0.2 for k in self.config.cation_list}
            return vals[j]

        self.membrane_convection_coefficient = Var(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            self.cations,
            initialize=initialize_membrane_convection_coefficient,
            units=units.dimensionless,
            doc="Convection coefficient for cations in membrane",
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
        if self.config.include_boundary_layer:
            self.d_boundary_layer_conc_mol_comp_dz = DerivativeVar(
                self.boundary_layer_conc_mol_comp,
                wrt=self.dimensionless_boundary_layer_thickness,
                units=units.mol / units.m**3,  # mM
                doc="Solute concentration gradient wrt z in the boundary layer",
            )
        self.d_membrane_conc_mol_comp_dz = DerivativeVar(
            self.membrane_conc_mol_comp,
            wrt=self.dimensionless_membrane_thickness,
            units=units.mol / units.m**3,  # mM
            doc="Solute concentration gradient wrt membrane thickness",
        )

    def add_constraints(self):
        """
        Adds model constraints for the multi-component diafiltration unit model.
        """

        # mol balance constraints
        def _overall_mol_balance(blk, t, x):
            if x == 0:
                return Constraint.Skip
            return blk.d_retentate_flow_volume_dx[t, x] == (
                -blk.volume_flux_water[t, x]
                * blk.total_membrane_length
                * blk.total_module_length
            )

        self.overall_mol_balance = Constraint(
            self.time, self.dimensionless_module_length, rule=_overall_mol_balance
        )

        def _cation_mol_balance(blk, t, x, k):
            if x == 0:
                return Constraint.Skip
            return (
                blk.retentate_flow_volume[t, x]
                * blk.d_retentate_conc_mol_comp_dx[t, x, k]
            ) == (
                (
                    blk.volume_flux_water[t, x] * blk.retentate_conc_mol_comp[t, x, k]
                    - blk.molar_ion_flux[t, x, k]
                )
                * blk.total_membrane_length
                * blk.total_module_length
            )

        self.cation_mol_balance = Constraint(
            self.time,
            self.dimensionless_module_length,
            self.cations,
            rule=_cation_mol_balance,
        )

        # bulk flux balance constraints
        def _overall_bulk_flux_equation(blk, t, x):
            if x == 0:
                return Constraint.Skip
            return (
                blk.permeate_flow_volume[t, x]
                == blk.volume_flux_water[t, x]
                * x
                * blk.total_membrane_length
                * blk.total_module_length
            )

        self.overall_bulk_flux_equation = Constraint(
            self.time,
            self.dimensionless_module_length,
            rule=_overall_bulk_flux_equation,
        )

        def _cation_bulk_flux_equation(blk, t, x, k):
            if x == 0:
                return Constraint.Skip
            return blk.molar_ion_flux[t, x, k] == (
                blk.permeate_conc_mol_comp[t, x, k] * blk.volume_flux_water[t, x]
            )

        self.cation_bulk_flux_equation = Constraint(
            self.time,
            self.dimensionless_module_length,
            self.cations,
            rule=_cation_bulk_flux_equation,
        )

        # transport constraints (first principles)
        def _lumped_water_flux(blk, t, x):
            if x == 0:
                return Constraint.Skip
            return blk.volume_flux_water[t, x] == (
                blk.membrane_permeability
                * (blk.applied_pressure[t] - blk.osmotic_pressure[t, x])
            )

        self.lumped_water_flux = Constraint(
            self.time, self.dimensionless_module_length, rule=_lumped_water_flux
        )

        if self.config.include_boundary_layer:

            def _boundary_layer_D_tilde_calculation(blk, t, x, z):
                if x == 0:
                    return Constraint.Skip
                a0 = self.config.anion_list[0]
                charge = blk.config.property_package.charge
                conc_bl = blk.boundary_layer_conc_mol_comp
                D_bl = blk.config.property_package.boundary_layer_diffusion_coefficient
                return blk.boundary_layer_D_tilde[t, x, z] == sum(
                    (
                        (
                            ((charge[k] ** 2) * D_bl[k])
                            - (charge[k] * charge[a0] * D_bl[a0])
                        )
                        * conc_bl[t, x, z, k]
                    )
                    for k in self.cations
                )

            self.boundary_layer_D_tilde_calculation = Constraint(
                self.time,
                self.dimensionless_module_length,
                self.dimensionless_boundary_layer_thickness,
                rule=_boundary_layer_D_tilde_calculation,
            )

            def _boundary_layer_cross_diffusion_coefficient_bilinear_calculation(
                blk, t, x, z, k, j
            ):
                if x == 0:
                    return Constraint.Skip
                return (
                    blk.boundary_layer_cross_diffusion_coefficient_bilinear[
                        t, x, z, k, j
                    ]
                    == blk.boundary_layer_cross_diffusion_coefficient[t, x, z, k, j]
                    * blk.boundary_layer_D_tilde[t, x, z]
                )

            self.boundary_layer_cross_diffusion_coefficient_bilinear_calculation = Constraint(
                self.time,
                self.dimensionless_module_length,
                self.dimensionless_boundary_layer_thickness,
                self.cations,
                self.cations,
                rule=_boundary_layer_cross_diffusion_coefficient_bilinear_calculation,
            )

            def _boundary_layer_cross_diffusion_coefficient_calculation(
                blk, t, x, z, k, j
            ):
                if x == 0:
                    return Constraint.Skip
                a0 = self.config.anion_list[0]
                charge = blk.config.property_package.charge
                conc_bl = blk.boundary_layer_conc_mol_comp
                D_bl = blk.config.property_package.boundary_layer_diffusion_coefficient
                # off-diagonal
                if k != j:
                    return blk.boundary_layer_cross_diffusion_coefficient_bilinear[
                        t, x, z, k, j
                    ] == (
                        (
                            (charge[k] * charge[j] * D_bl[k] * D_bl[j])
                            - (charge[k] * charge[j] * D_bl[k] * D_bl[a0])
                        )
                        * conc_bl[t, x, z, k]
                    )
                # diagonal
                if k == j:
                    return blk.boundary_layer_cross_diffusion_coefficient_bilinear[
                        t, x, z, k, j
                    ] == (
                        sum(
                            (
                                (
                                    (charge[i] * charge[a0] * D_bl[k] * D_bl[a0])
                                    - (charge[i] ** 2 * D_bl[i] * D_bl[k])
                                )
                                * conc_bl[t, x, z, i]
                            )
                            for i in blk.cations
                            if k != i
                        )
                        + sum(
                            (
                                (
                                    (charge[i] * charge[a0] * D_bl[k] * D_bl[a0])
                                    - (charge[i] ** 2 * D_bl[i] * D_bl[a0])
                                )
                                * conc_bl[t, x, z, i]
                            )
                            for i in blk.cations
                            if k == i
                        )
                    )

            self.boundary_layer_cross_diffusion_coefficient_calculation = Constraint(
                self.time,
                self.dimensionless_module_length,
                self.dimensionless_boundary_layer_thickness,
                self.cations,
                self.cations,
                rule=_boundary_layer_cross_diffusion_coefficient_calculation,
            )

        def _membrane_D_tilde_calculation(blk, t, x, z):
            if x == 0:
                return Constraint.Skip
            a0 = self.config.anion_list[0]
            charge = blk.config.property_package.charge
            chi = blk.membrane_fixed_charge
            conc_mem = blk.membrane_conc_mol_comp
            D_mem = blk.config.property_package.membrane_diffusion_coefficient
            return blk.membrane_D_tilde[t, x, z] == (
                sum(
                    (
                        (
                            ((charge[k] ** 2) * D_mem[k])
                            - (charge[k] * charge[a0] * D_mem[a0])
                        )
                        * conc_mem[t, x, z, k]
                    )
                    for k in blk.cations
                )
                - (charge[a0] * D_mem[a0] * chi)
            )

        self.membrane_D_tilde_calculation = Constraint(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            rule=_membrane_D_tilde_calculation,
        )

        def _membrane_cross_diffusion_coefficient_bilinear_calculation(
            blk, t, x, z, k, j
        ):
            if x == 0:
                return Constraint.Skip
            return (
                blk.membrane_cross_diffusion_coefficient_bilinear[t, x, z, k, j]
                == blk.membrane_cross_diffusion_coefficient[t, x, z, k, j]
                * blk.membrane_D_tilde[t, x, z]
            )

        self.membrane_cross_diffusion_coefficient_bilinear_calculation = Constraint(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            self.cations,
            self.cations,
            rule=_membrane_cross_diffusion_coefficient_bilinear_calculation,
        )

        def _membrane_convection_coefficient_bilinear_calculation(blk, t, x, z, k):
            if x == 0:
                return Constraint.Skip
            return (
                blk.membrane_convection_coefficient_bilinear[t, x, z, k]
                == blk.membrane_convection_coefficient[t, x, z, k]
                * blk.membrane_D_tilde[t, x, z]
            )

        self.membrane_convection_coefficient_bilinear_calculation = Constraint(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            self.cations,
            rule=_membrane_convection_coefficient_bilinear_calculation,
        )

        def _membrane_cross_diffusion_coefficient_calculation(blk, t, x, z, k, j):
            if x == 0:
                return Constraint.Skip
            a0 = self.config.anion_list[0]
            charge = blk.config.property_package.charge
            chi = blk.membrane_fixed_charge
            conc_mem = blk.membrane_conc_mol_comp
            D_mem = blk.config.property_package.membrane_diffusion_coefficient
            # off-diagonal
            if k != j:
                return blk.membrane_cross_diffusion_coefficient_bilinear[
                    t, x, z, k, j
                ] == (
                    (
                        (charge[k] * charge[j] * D_mem[k] * D_mem[j])
                        - (charge[k] * charge[j] * D_mem[k] * D_mem[a0])
                    )
                    * conc_mem[t, x, z, k]
                )
            # diagonal
            if k == j:
                return blk.membrane_cross_diffusion_coefficient_bilinear[
                    t, x, z, k, j
                ] == (
                    sum(
                        (
                            (
                                (charge[i] * charge[a0] * D_mem[k] * D_mem[a0])
                                - (charge[i] ** 2 * D_mem[i] * D_mem[k])
                            )
                            * conc_mem[t, x, z, i]
                        )
                        for i in blk.cations
                        if k != i
                    )
                    + sum(
                        (
                            (
                                (charge[i] * charge[a0] * D_mem[k] * D_mem[a0])
                                - (charge[i] ** 2 * D_mem[i] * D_mem[a0])
                            )
                            * conc_mem[t, x, z, i]
                        )
                        for i in blk.cations
                        if k == i
                    )
                    + charge[a0] * D_mem[k] * D_mem[a0] * chi
                )

        self.membrane_cross_diffusion_coefficient_calculation = Constraint(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            self.cations,
            self.cations,
            rule=_membrane_cross_diffusion_coefficient_calculation,
        )

        def _membrane_convection_coefficient_calculation(blk, t, x, z, k):
            if x == 0:
                return Constraint.Skip
            charge = blk.config.property_package.charge
            chi = blk.membrane_fixed_charge
            D_mem = blk.config.property_package.membrane_diffusion_coefficient
            return blk.membrane_convection_coefficient_bilinear[t, x, z, k] == (
                blk.membrane_D_tilde[t, x, z] + (charge[k] * D_mem[k] * chi)
            )

        self.membrane_convection_coefficient_calculation = Constraint(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            self.cations,
            rule=_membrane_convection_coefficient_calculation,
        )

        if self.config.include_boundary_layer:

            def _cation_flux_boundary_layer(blk, t, x, z, k):
                if x == 0 or z == 0:
                    return Constraint.Skip
                return blk.molar_ion_flux[t, x, k] == (
                    (
                        blk.boundary_layer_conc_mol_comp[t, x, z, k]
                        * blk.volume_flux_water[t, x]
                    )
                    + sum(
                        (
                            units.convert(
                                blk.boundary_layer_cross_diffusion_coefficient[
                                    t, x, z, k, i
                                ],
                                to_units=units.m**2 / units.h,
                            )
                            / (blk.total_boundary_layer_thickness)
                            * blk.d_boundary_layer_conc_mol_comp_dz[t, x, z, i]
                        )
                        for i in self.cations
                    )
                )

            self.cation_flux_boundary_layer = Constraint(
                self.time,
                self.dimensionless_module_length,
                self.dimensionless_boundary_layer_thickness,
                self.cations,
                rule=_cation_flux_boundary_layer,
            )

        def _cation_flux_membrane(blk, t, x, z, k):
            if x == 0:
                return Constraint.Skip
            return blk.molar_ion_flux[t, x, k] == (
                (
                    blk.membrane_convection_coefficient[t, x, z, k]
                    * blk.membrane_conc_mol_comp[t, x, z, k]
                    * blk.volume_flux_water[t, x]
                )
                + sum(
                    (
                        units.convert(
                            blk.membrane_cross_diffusion_coefficient[t, x, z, k, i],
                            to_units=units.m**2 / units.h,
                        )
                        / blk.total_membrane_thickness
                        * blk.d_membrane_conc_mol_comp_dz[t, x, z, i]
                    )
                    for i in blk.cations
                )
            )

        self.cation_flux_membrane = Constraint(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            self.cations,
            rule=_cation_flux_membrane,
        )

        def _anion_flux_membrane(blk, t, x):
            if x == 0:
                return Constraint.Skip
            charge = blk.config.property_package.charge
            return 0 == sum(
                charge[j] * blk.molar_ion_flux[t, x, j] for j in blk.solutes
            )

        self.anion_flux_membrane = Constraint(
            self.time, self.dimensionless_module_length, rule=_anion_flux_membrane
        )

        # other physical constraints
        def _osmotic_pressure_calculation(blk, t, x):
            if x == 0:
                return Constraint.Skip
            conc_p = blk.permeate_conc_mol_comp
            conc_r = blk.retentate_conc_mol_comp
            n = blk.config.property_package.num_solutes
            R = Constants.gas_constant  # J / mol / K
            sigma = blk.config.property_package.sigma
            T = blk.temperature
            if self.config.include_boundary_layer:
                conc_bl = blk.boundary_layer_conc_mol_comp
                return blk.osmotic_pressure[t, x] == units.convert(
                    (
                        R
                        * T
                        * sum(
                            (n[j] * sigma[j] * (conc_bl[t, x, 1, j] - conc_p[t, x, j]))
                            for j in blk.solutes
                        )
                    ),
                    to_units=units.bar,
                )
            else:
                return blk.osmotic_pressure[t, x] == units.convert(
                    (
                        R
                        * T
                        * sum(
                            (n[j] * sigma[j] * (conc_r[t, x, j] - conc_p[t, x, j]))
                            for j in blk.solutes
                        )
                    ),
                    to_units=units.bar,
                )

        self.osmotic_pressure_calculation = Constraint(
            self.time,
            self.dimensionless_module_length,
            rule=_osmotic_pressure_calculation,
        )

        def _electroneutrality_retentate(blk, t, x):
            charge = blk.config.property_package.charge
            conc_r = blk.retentate_conc_mol_comp
            return 0 == sum(charge[j] * conc_r[t, x, j] for j in blk.solutes)

        self.electroneutrality_retentate = Constraint(
            self.time,
            self.dimensionless_module_length,
            rule=_electroneutrality_retentate,
        )

        if self.config.include_boundary_layer:

            def _electroneutrality_boundary_layer(blk, t, x, z):
                if x == 0:
                    return Constraint.Skip
                charge = blk.config.property_package.charge
                conc_bl = blk.boundary_layer_conc_mol_comp
                return 0 == sum(charge[j] * conc_bl[t, x, z, j] for j in blk.solutes)

            self.electroneutrality_boundary_layer = Constraint(
                self.time,
                self.dimensionless_module_length,
                self.dimensionless_boundary_layer_thickness,
                rule=_electroneutrality_boundary_layer,
            )

        def _electroneutrality_membrane(blk, t, x, z):
            if x == 0:
                return Constraint.Skip
            charge = blk.config.property_package.charge
            chi = blk.membrane_fixed_charge
            conc_mem = blk.membrane_conc_mol_comp
            return 0 == (
                sum(charge[j] * conc_mem[t, x, z, j] for j in blk.solutes) + chi
            )

        self.electroneutrality_membrane = Constraint(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            rule=_electroneutrality_membrane,
        )

        def _electroneutrality_permeate(blk, t, x):
            if x == 0:
                return Constraint.Skip
            charge = blk.config.property_package.charge
            conc_p = blk.permeate_conc_mol_comp
            return 0 == sum(charge[j] * conc_p[t, x, j] for j in blk.solutes)

        self.electroneutrality_permeate = Constraint(
            self.time,
            self.dimensionless_module_length,
            rule=_electroneutrality_permeate,
        )

        # partitioning equations
        if self.config.include_boundary_layer:

            def _retentate_boundary_layer_interface(blk, t, x, k):
                if x == 0:
                    return Constraint.Skip
                return (
                    blk.retentate_conc_mol_comp[t, x, k]
                    == blk.boundary_layer_conc_mol_comp[t, x, 0, k]
                )

            self.retentate_boundary_layer_interface = Constraint(
                self.time,
                self.dimensionless_module_length,
                self.cations,
                rule=_retentate_boundary_layer_interface,
            )

            def _cation_equilibrium_boundary_layer_membrane_interface(blk, t, x, k):
                if x == 0:
                    return Constraint.Skip
                a0 = self.config.anion_list[0]
                charge = blk.config.property_package.charge
                conc_bl = blk.boundary_layer_conc_mol_comp
                conc_mem = blk.membrane_conc_mol_comp
                H_r = blk.config.property_package.partition_coefficient_retentate
                return (
                    (H_r[k] ** (-charge[a0]))
                    * (H_r[a0] ** charge[k])
                    * (conc_bl[t, x, 1, k] ** (-charge[a0]))
                    * (conc_bl[t, x, 1, a0] ** charge[k])
                ) == (
                    (conc_mem[t, x, 0, k] ** (-charge[a0]))
                    * (conc_mem[t, x, 0, a0] ** charge[k])
                )

            self.cation_equilibrium_boundary_layer_membrane_interface = Constraint(
                self.time,
                self.dimensionless_module_length,
                self.cations,
                rule=_cation_equilibrium_boundary_layer_membrane_interface,
            )
        else:

            def _cation_equilibrium_retentate_membrane_interface(blk, t, x, k):
                if x == 0:
                    return Constraint.Skip
                a0 = self.config.anion_list[0]
                charge = blk.config.property_package.charge
                conc_mem = blk.membrane_conc_mol_comp
                conc_r = blk.retentate_conc_mol_comp
                H_r = blk.config.property_package.partition_coefficient_retentate
                return (
                    (H_r[k] ** (-charge[a0]))
                    * (H_r[a0] ** charge[k])
                    * (conc_r[t, x, k] ** (-charge[a0]))
                    * (conc_r[t, x, a0] ** charge[k])
                ) == (
                    (conc_mem[t, x, 0, k] ** (-charge[a0]))
                    * (conc_mem[t, x, 0, a0] ** charge[k])
                )

            self.cation_equilibrium_retentate_membrane_interface = Constraint(
                self.time,
                self.dimensionless_module_length,
                self.cations,
                rule=_cation_equilibrium_retentate_membrane_interface,
            )

        def _cation_equilibrium_membrane_permeate_interface(blk, t, x, k):
            if x == 0:
                return Constraint.Skip
            a0 = self.config.anion_list[0]
            charge = blk.config.property_package.charge
            conc_mem = blk.membrane_conc_mol_comp
            conc_p = blk.permeate_conc_mol_comp
            H_p = blk.config.property_package.partition_coefficient_permeate
            return (
                (H_p[k] ** (-charge[a0]))
                * (H_p[a0] ** charge[k])
                * (conc_p[t, x, k] ** (-charge[a0]))
                * (conc_p[t, x, a0] ** charge[k])
            ) == (
                (conc_mem[t, x, 1, k] ** (-charge[a0]))
                * (conc_mem[t, x, 1, a0] ** charge[k])
            )

        self.cation_equilibrium_membrane_permeate_interface = Constraint(
            self.time,
            self.dimensionless_module_length,
            self.cations,
            rule=_cation_equilibrium_membrane_permeate_interface,
        )

        # boundary conditions
        def _retentate_flow_volume_boundary_condition(blk, t):
            return (
                blk.retentate_flow_volume[t, 0]
                == blk.feed_flow_volume[t] + blk.diafiltrate_flow_volume[t]
            )

        self.retentate_flow_volume_boundary_condition = Constraint(
            self.time, rule=_retentate_flow_volume_boundary_condition
        )

        def _retentate_conc_mol_comp_boundary_condition(blk, t, k):
            return blk.retentate_conc_mol_comp[t, 0, k] == (
                (
                    blk.feed_flow_volume[t] * blk.feed_conc_mol_comp[t, k]
                    + blk.diafiltrate_flow_volume[t]
                    * blk.diafiltrate_conc_mol_comp[t, k]
                )
                / (blk.feed_flow_volume[t] + blk.diafiltrate_flow_volume[t])
            )

        self.retentate_conc_mol_comp_boundary_condition = Constraint(
            self.time, self.cations, rule=_retentate_conc_mol_comp_boundary_condition
        )

        if self.config.include_boundary_layer:

            def _boundary_layer_conc_mol_comp_boundary_condition(blk, t, z, k):
                return (
                    blk.boundary_layer_conc_mol_comp[t, 0, z, k]
                    == self.numerical_zero_tolerance * units.mol / units.m**3
                )

            self.boundary_layer_conc_mol_comp_boundary_condition = Constraint(
                self.time,
                self.dimensionless_boundary_layer_thickness,
                self.cations,
                rule=_boundary_layer_conc_mol_comp_boundary_condition,
            )

        def _membrane_conc_mol_comp_boundary_condition(blk, t, z, k):
            return (
                blk.membrane_conc_mol_comp[t, 0, z, k]
                == self.numerical_zero_tolerance * units.mol / units.m**3
            )

        self.membrane_conc_mol_comp_boundary_condition = Constraint(
            self.time,
            self.dimensionless_membrane_thickness,
            self.cations,
            rule=_membrane_conc_mol_comp_boundary_condition,
        )

        # constraints to improve numerical stability
        def _permeate_flow_volume_boundary_condition(blk, t):
            return (
                blk.permeate_flow_volume[t, 0]
                == self.numerical_zero_tolerance * units.m**3 / units.h
            )

        self.permeate_flow_volume_boundary_condition = Constraint(
            self.time, rule=_permeate_flow_volume_boundary_condition
        )

        def _permeate_conc_mol_comp_boundary_condition(blk, t, j):
            return (
                blk.permeate_conc_mol_comp[t, 0, j]
                == self.numerical_zero_tolerance * units.mol / units.m**3
            )

        self.permeate_conc_mol_comp_boundary_condition = Constraint(
            self.time, self.solutes, rule=_permeate_conc_mol_comp_boundary_condition
        )

        def _d_retentate_flow_volume_dx_boundary_condition(blk, t):
            return (
                blk.d_retentate_flow_volume_dx[t, 0]
                == self.numerical_zero_tolerance * units.m**3 / units.h
            )

        self.d_retentate_flow_volume_dx_boundary_condition = Constraint(
            self.time, rule=_d_retentate_flow_volume_dx_boundary_condition
        )

        def _d_retentate_conc_mol_comp_dx_boundary_condition(blk, t, k):
            return (
                blk.d_retentate_conc_mol_comp_dx[t, 0, k]
                == self.numerical_zero_tolerance * units.mol / units.m**3
            )

        self.d_retentate_conc_mol_comp_dx_boundary_condition = Constraint(
            self.time,
            self.cations,
            rule=_d_retentate_conc_mol_comp_dx_boundary_condition,
        )

        def _volume_flux_water_boundary_condition(blk, t):
            return (
                blk.volume_flux_water[t, 0]
                == self.numerical_zero_tolerance * units.m / units.h
            )

        self.volume_flux_water_boundary_condition = Constraint(
            self.time, rule=_volume_flux_water_boundary_condition
        )

        def _molar_ion_flux_boundary_condition(blk, t, j):
            return (
                blk.molar_ion_flux[t, 0, j]
                == self.numerical_zero_tolerance * units.mol / units.m**2 / units.h
            )

        self.molar_ion_flux_boundary_condition = Constraint(
            self.time, self.solutes, rule=_molar_ion_flux_boundary_condition
        )

    def discretize_model(self):
        discretizer = TransformationFactory("dae.finite_difference")
        discretizer.apply_to(
            self,
            wrt=self.dimensionless_module_length,
            nfe=self.config.NFE_module_length,
            scheme="BACKWARD",
        )
        if self.config.include_boundary_layer:
            discretizer.apply_to(
                self,
                wrt=self.dimensionless_boundary_layer_thickness,
                nfe=self.config.NFE_boundary_layer_thickness,
                scheme="BACKWARD",
            )
        discretizer.apply_to(
            self,
            wrt=self.dimensionless_membrane_thickness,
            nfe=self.config.NFE_membrane_thickness,
            scheme="BACKWARD",
        )

    def deactivate_unnecessary_objects(self):
        """
        Deactivates variables and constraints not needed in the multi-component
        diafiltration unit model.
        """
        a0 = self.config.anion_list[0]
        for t in self.time:
            for x in self.dimensionless_module_length:
                # anion concentration gradient in retentate variable is created by default but
                # is not needed in model; fix to reduce number of variables
                self.d_retentate_conc_mol_comp_dx[t, x, a0].fix(
                    value(self.numerical_zero_tolerance)
                )
                # associated discretization equation not needed in model
                if x != 0:
                    self.d_retentate_conc_mol_comp_dx_disc_eq[t, x, a0].deactivate()

                if self.config.include_boundary_layer:
                    for z in self.dimensionless_boundary_layer_thickness:
                        # anion concentration gradient in boundary layer variable is created by default but
                        # is not needed in model; fix to reduce number of variables
                        self.d_boundary_layer_conc_mol_comp_dz[t, x, z, a0].fix(
                            value(self.numerical_zero_tolerance)
                        )
                        # associated discretization equation not needed in model
                        if z != 0:
                            self.d_boundary_layer_conc_mol_comp_dz_disc_eq[
                                t, x, z, a0
                            ].deactivate()
                for z in self.dimensionless_membrane_thickness:
                    # anion concentration gradient in membrane variable is created by default but
                    # is not needed in model; fix to reduce number of variables
                    self.d_membrane_conc_mol_comp_dz[t, x, z, a0].fix(
                        value(self.numerical_zero_tolerance)
                    )
                    # associated discretization equation not needed in model
                    if z != 0:
                        self.d_membrane_conc_mol_comp_dz_disc_eq[
                            t, x, z, a0
                        ].deactivate()

    def add_scaling_factors(self):
        """
        Assigns scaling factors to certain variables and constraints to
        improve solver performance.
        """
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        self.scaling_factor[self.volume_flux_water] = 1e2
        if self.config.include_boundary_layer:
            self.scaling_factor[self.boundary_layer_D_tilde] = 1e-3
            self.scaling_factor[
                self.boundary_layer_cross_diffusion_coefficient_bilinear
            ] = 1e-4
            self.scaling_factor[self.boundary_layer_cross_diffusion_coefficient] = 1e1
        self.scaling_factor[self.membrane_D_tilde] = 1e-1
        self.scaling_factor[self.membrane_cross_diffusion_coefficient_bilinear] = 1e-2
        self.scaling_factor[self.membrane_convection_coefficient_bilinear] = 1e-1
        self.scaling_factor[self.membrane_cross_diffusion_coefficient] = 1e1
        self.scaling_factor[self.membrane_convection_coefficient] = 1e1

        if len(self.config.cation_list) >= 2:
            for t in self.time:
                for x in self.dimensionless_module_length:
                    if x != 0:
                        self.scaling_factor[self.lumped_water_flux[t, x]] = 1e3

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

    def add_helpful_expressions(self):
        def _feed_ionic_strength(
            blk,
            t,
        ):
            charge = blk.config.property_package.charge
            return 0.5 * sum(
                (
                    (
                        (
                            blk.feed_flow_volume[t] * blk.feed_conc_mol_comp[t, j]
                            + blk.diafiltrate_flow_volume[t]
                            * blk.diafiltrate_conc_mol_comp[t, j]
                        )
                        / (blk.feed_flow_volume[t] + blk.diafiltrate_flow_volume[t])
                    )
                    * charge[j] ** 2
                )
                for j in blk.solutes
            )

        self.feed_ionic_strength = Expression(self.time, rule=_feed_ionic_strength)
