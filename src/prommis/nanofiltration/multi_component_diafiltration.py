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
    exp,
    units,
    value,
)
from pyomo.network import Port
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core import UnitModelBlockData, declare_process_block_class, useDefault
from idaes.core.initialization import InitializerBase
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.constants import Constants
from idaes.core.util.exceptions import ConfigurationError


class MultiComponentDiafiltrationInitializer(InitializerBase):
    """
    Multi-Component Diafiltration Initializer Class.
    """

    CONFIG = InitializerBase.CONFIG()

    CONFIG.declare(
        "multiplier_H_feed",
        ConfigValue(
            default=1,
            doc="Multiplicative factor to adjust H_feed guess",
        ),
    )
    CONFIG.declare(
        "multiplier_H_perm",
        ConfigValue(
            default=1,
            doc="Multiplicative factor to adjust H_perm guess",
        ),
    )

    def initialize(self, model):
        """
        Uses a combination of heuristics and data-based surrogates for
        important concentration-based parameters (the partition and
        sieving coefficients) to populate the model with a
        reasonable initial point.
        """
        # define helpful aliases
        alpha_mem_bi = model.membrane_convection_coefficient_bilinear
        alpha_mem_bi_calc = model.membrane_convection_coefficient_bilinear_calculation
        alpha_mem = model.membrane_convection_coefficient
        alpha_mem_calc = model.membrane_convection_coefficient_calculation
        a0 = model.config.anion_list[0]
        boundary_layer = model.config.include_boundary_layer
        charge = model.config.property_package.charge
        conc_f_tot = model.total_feed_conc_mol_comp
        conc_ret = model.retentate_conc_mol_comp
        conc_perm = model.permeate_conc_mol_comp
        conc_mem = model.membrane_conc_mol_comp
        d_conc_ret_dx = model.d_retentate_conc_mol_comp_dx
        d_conc_mem_dz = model.d_membrane_conc_mol_comp_dz
        d_q_r_dx = model.d_retentate_flow_volume_dx
        D_mem_bi = model.membrane_cross_diffusion_coefficient_bilinear
        D_mem_bi_calc = model.membrane_cross_diffusion_coefficient_bilinear_calculation
        D_mem = model.membrane_cross_diffusion_coefficient
        D_mem_calc = model.membrane_cross_diffusion_coefficient_calculation
        q_f_tot = model.total_feed_flow_volume
        q_ret = model.retentate_flow_volume
        q_perm = model.permeate_flow_volume
        zero_val = value(model.numerical_zero_tolerance)
        if boundary_layer:
            conc_bl = model.boundary_layer_conc_mol_comp
            d_conc_bl_dz = model.d_boundary_layer_conc_mol_comp_dz
            D_bl_bi = model.boundary_layer_cross_diffusion_coefficient_bilinear
            D_bl_bi_calc = (
                model.boundary_layer_cross_diffusion_coefficient_bilinear_calculation
            )
            D_bl = model.boundary_layer_cross_diffusion_coefficient
            D_bl_calc = model.boundary_layer_cross_diffusion_coefficient_calculation

        for t in model.time:
            # set initial conditions
            q_ret[t, 0].set_value(value(q_f_tot[t]))
            d_q_r_dx[t, 0].set_value(zero_val)
            q_perm[t, 0].set_value(zero_val)
            model.volume_flux_water[t, 0].set_value(zero_val)

            for j in model.solutes:
                conc_ret[t, 0, j].set_value(value(conc_f_tot[t, j]))
                d_conc_ret_dx[t, 0, j].set_value(zero_val)
                conc_perm[t, 0, j].set_value(zero_val)
                if boundary_layer:
                    for z_bl in model.dimensionless_boundary_layer_thickness:
                        conc_bl[t, 0, z_bl, j].set_value(zero_val)
                        d_conc_bl_dz[t, 0, z_bl, j].set_value(zero_val)
                for z_m in model.dimensionless_membrane_thickness:
                    conc_mem[t, 0, z_m, j].set_value(zero_val)
                    d_conc_mem_dz[t, 0, z_m, j].set_value(zero_val)
                model.molar_ion_flux[t, 0, j].set_value(zero_val)

            x_prev = 0
            for x in model.dimensionless_module_length:
                if x != 0:
                    # temporary retentate conditions
                    q_ret[t, x].set_value(value(q_ret[t, x_prev]))
                    for j in model.solutes:
                        conc_ret[t, x, j].set_value(value(conc_ret[t, x_prev, j]))

                    # guess permeate concentration with observed sieving coefficient surrogate
                    # S_i_obs = c_{i,p} / c_{i,r}
                    # a, b, and c are fitted from successful "brute force" solves
                    for k in model.cations:
                        if value(charge[k]) == 3:
                            params = {"a": -0.09, "b": -0.021, "c": 0.14}
                        elif value(charge[k]) == 2:
                            params = {"a": -0.33, "b": -0.045, "c": 0.51}
                        elif value(charge[k]) == 1:
                            params = {"a": -0.71, "b": -0.031, "c": 0.69}
                        S_guess = (
                            params["a"] * exp(params["b"] * value(conc_f_tot[t, k]))
                            + params["c"]
                        )
                        conc_perm[t, x, k].set_value(
                            value(conc_ret[t, x, k]) * round(S_guess, 1)
                        )
                    calculate_variable_from_constraint(
                        conc_perm[t, x, a0],
                        model.electroneutrality_permeate[t, x],
                    )
                    if boundary_layer:
                        # geuss interface concentration with constant CP modulus = 1.01
                        for k in model.cations:
                            conc_bl[t, x, 1, k].set_value(
                                value(conc_ret[t, x, k]) * 1.01
                            )
                        calculate_variable_from_constraint(
                            conc_bl[t, x, 1, a0],
                            model.electroneutrality_boundary_layer[t, x, 1],
                        )
                    # calculate osmotic pressure and fluxes
                    calculate_variable_from_constraint(
                        model.osmotic_pressure[t, x],
                        model.osmotic_pressure_calculation[t, x],
                    )
                    calculate_variable_from_constraint(
                        model.volume_flux_water[t, x], model.lumped_water_flux[t, x]
                    )
                    for k in model.cations:
                        # calculate molar flux assuming convection only
                        model.molar_ion_flux[t, x, k].set_value(
                            value(conc_perm[t, x, k])
                            * value(model.volume_flux_water[t, x])
                        )
                    calculate_variable_from_constraint(
                        model.molar_ion_flux[t, x, a0], model.anion_flux_membrane[t, x]
                    )
                    # calculate flow rates
                    calculate_variable_from_constraint(
                        q_perm[t, x],
                        model.overall_mass_balance[t, x],
                    )
                    q_ret[t, x].set_value(value(q_f_tot[t]) - value(q_perm[t, x]))
                    # calculate derivatives
                    calculate_variable_from_constraint(
                        d_q_r_dx[t, x],
                        model.overall_mass_balance[t, x],
                    )
                    # d_qr / d_x = (qr(x) - qr(x_prev)) / (x - x_prev)
                    # (d_qr / d_x)*(x - x_prev) + qr(x_prev) = qr(x)
                    q_ret[t, x].set_value(
                        (value(d_q_r_dx[t, x]) * (x - x_prev)) + value(q_ret[t, x_prev])
                    )
                    for k in model.cations:
                        calculate_variable_from_constraint(
                            d_conc_ret_dx[t, x, k],
                            model.cation_mol_balance[t, x, k],
                        )
                        # d_cr / d_x = (cr(x) - cr(x_prev)) / (x - x_prev)
                        # (d_cr / d_x)*(x - x_prev) + cr(x_prev) = cr(x)
                        conc_ret[t, x, k].set_value(
                            (value(d_conc_ret_dx[t, x, k]) * (x - x_prev))
                            + value(conc_ret[t, x_prev, k])
                        )
                        calculate_variable_from_constraint(
                            conc_ret[t, x, a0],
                            model.electroneutrality_retentate[t, x],
                        )

                    if boundary_layer:
                        z_bl_prev = 0
                        for z_bl in model.dimensionless_boundary_layer_thickness:
                            if z_bl == 0:
                                for j in model.solutes:
                                    conc_bl[t, x, z_bl, j].set_value(
                                        value(conc_ret[t, x, j])
                                    )
                            else:
                                # guess boundary layer concentrations with linear profile
                                # slope = (c_int - c_r) / (1 - 0)
                                # c_bl = slope * z_bl + c_r
                                for k in model.cations:
                                    slope = value(conc_bl[t, x, 1, k]) - value(
                                        conc_bl[t, x, 0, k]
                                    )
                                    conc_bl[t, x, z_bl, k].set_value(
                                        slope * z_bl + value(conc_bl[t, x, 0, k])
                                    )
                                calculate_variable_from_constraint(
                                    conc_bl[t, x, z_bl, a0],
                                    model.electroneutrality_boundary_layer[t, x, z_bl],
                                )
                                for j in model.solutes:
                                    # d_c_bl / d_z_bl = (c_bl(z_bl) - c_bl(z_bl_prev)) / (z_bl - z_bl_prev)
                                    d_conc_bl_dz[t, x, z_bl, j].set_value(
                                        (
                                            value(conc_bl[t, x, z_bl, j])
                                            - value(conc_bl[t, x, z_bl_prev, j])
                                        )
                                        / (z_bl - z_bl_prev)
                                    )
                                    if z_bl_prev == 0:
                                        d_conc_bl_dz[t, x, z_bl_prev, j].set_value(
                                            value(d_conc_bl_dz[t, x, z_bl, j])
                                        )
                            # update diffusion coefficients
                            calculate_variable_from_constraint(
                                model.boundary_layer_D_tilde[t, x, z_bl],
                                model.boundary_layer_D_tilde_calculation[t, x, z_bl],
                            )
                            for k in model.cations:
                                for j in model.cations:
                                    calculate_variable_from_constraint(
                                        D_bl_bi[t, x, z_bl, k, j],
                                        D_bl_calc[t, x, z_bl, k, j],
                                    )
                                    calculate_variable_from_constraint(
                                        D_bl[t, x, z_bl, k, j],
                                        D_bl_bi_calc[t, x, z_bl, k, j],
                                    )
                            z_bl_prev = z_bl
                    z_m_prev = 0
                    # guess interface concentration with partition coefficient surrogate
                    # H_i = c_{i,m} / c_{i,s}
                    # a, b, and c are fitted from successful "brute force" solves
                    for k in model.cations:
                        if value(charge[k]) == 3:
                            feed_params = {"a": 2, "b": -0.07, "c": 0.2}
                            perm_params = {"a": 50, "b": -0.09, "c": 10}
                        elif value(charge[k]) == 2:
                            feed_params = {"a": 3, "b": -0.07, "c": 0.3}
                            perm_params = {"a": 20, "b": -0.09, "c": 0.5}
                        elif value(charge[k]) == 1:
                            feed_params = {"a": 6, "b": -0.07, "c": 0.5}
                            perm_params = {"a": 80, "b": -0.1, "c": 0.9}

                        H_feed = (
                            feed_params["a"]
                            * exp(feed_params["b"] * value(conc_f_tot[t, k]))
                            + feed_params["c"]
                        )
                        H_feed = value(self.config.multiplier_H_feed) * H_feed

                        H_perm = (
                            perm_params["a"]
                            * exp(perm_params["b"] * value(conc_f_tot[t, k]))
                            + perm_params["c"]
                        )
                        H_perm = value(self.config.multiplier_H_perm) * H_perm

                        conc_mem[t, x, 0, k].set_value(
                            round(H_feed, 1) * value(conc_ret[t, x, k])
                        )
                        conc_mem[t, x, 1, k].set_value(
                            round(H_perm, 1) * value(conc_perm[t, x, k])
                        )
                    calculate_variable_from_constraint(
                        conc_mem[t, x, 0, a0],
                        model.electroneutrality_membrane[t, x, 0],
                    )

                    for z_m in model.dimensionless_membrane_thickness:
                        if z_m != 0:
                            # guess boundary layer concentrations with linear profile
                            # slope = (c_m(1) - c_m(0)) / (1 - 0)
                            # c_m(x) = slope * z_m + c_m(0)
                            for k in model.cations:
                                slope = value(conc_mem[t, x, 1, k]) - value(
                                    conc_mem[t, x, 0, k]
                                )
                                conc_mem[t, x, z_m, k].set_value(
                                    slope * z_m + value(conc_mem[t, x, 0, k])
                                )
                            calculate_variable_from_constraint(
                                conc_mem[t, x, z_m, a0],
                                model.electroneutrality_membrane[t, x, z_m],
                            )
                            for j in model.solutes:
                                # d_cm / d_z_m = (c_m(z_m) - c_m(z_m_prev)) / (z_m - z_m_prev)
                                d_conc_mem_dz[t, x, z_m, j].set_value(
                                    (
                                        value(conc_mem[t, x, z_m, j])
                                        - value(conc_mem[t, x, z_m_prev, j])
                                    )
                                    / (z_m - z_m_prev)
                                )
                                if z_m_prev == 0:
                                    d_conc_mem_dz[t, x, z_m_prev, j].set_value(
                                        value(d_conc_mem_dz[t, x, z_m, j])
                                    )

                        # update diffusion and convection coefficients
                        calculate_variable_from_constraint(
                            model.membrane_D_tilde[t, x, z_m],
                            model.membrane_D_tilde_calculation[t, x, z_m],
                        )
                        for k in model.cations:
                            calculate_variable_from_constraint(
                                alpha_mem_bi[t, x, z_m, k],
                                alpha_mem_calc[t, x, z_m, k],
                            )
                            calculate_variable_from_constraint(
                                alpha_mem[t, x, z_m, k],
                                alpha_mem_bi_calc[t, x, z_m, k],
                            )
                            for j in model.cations:
                                calculate_variable_from_constraint(
                                    D_mem_bi[t, x, z_m, k, j],
                                    D_mem_calc[t, x, z_m, k, j],
                                )
                                calculate_variable_from_constraint(
                                    D_mem[t, x, z_m, k, j],
                                    D_mem_bi_calc[t, x, z_m, k, j],
                                )
                        z_m_prev = z_m
                x_prev = x


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
        "total_membrane_thickness",
        ConfigValue(
            default=1e-7,
            doc="Total membrane thickness in m (in the z-direction)",
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
            initialize=self.config.total_membrane_thickness,
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
            initialize=11,
            mutable=True,
            units=units.L / units.m**2 / units.h / units.bar,
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
            initialize=20,
            units=units.bar,
            bounds=[1e-11, 41],  # maximum operating presssure (NF270-440)
            doc="Pressure applied to membrane",
        )
        self.feed_flow_volume = Var(
            self.time,
            units=units.m**3 / units.h,
            bounds=[1e-11, None],
            doc="Volumetric flow rate of the feed",
        )
        self.feed_conc_mol_comp = Var(
            self.time,
            self.solutes,
            units=units.mol / units.m**3,  # mM
            bounds=[1e-11, None],
            doc="Mole concentration of solutes in the feed",
        )
        self.diafiltrate_flow_volume = Var(
            self.time,
            units=units.m**3 / units.h,
            bounds=[1e-11, None],
            doc="Volumetric flow rate of the diafiltrate",
        )
        self.diafiltrate_conc_mol_comp = Var(
            self.time,
            self.solutes,
            units=units.mol / units.m**3,  # mM
            bounds=[1e-11, None],
            doc="Mole concentration of solutes in the diafiltrate",
        )

        # add variables dependent on dimensionless_module_length
        self.volume_flux_water = Var(
            self.time,
            self.dimensionless_module_length,
            units=units.m**3 / units.m**2 / units.h,
            bounds=[1e-11, None],
            doc="Volumetric water flux of water across the membrane",
        )
        self.molar_ion_flux = Var(
            self.time,
            self.dimensionless_module_length,
            self.solutes,
            units=units.mol / units.m**2 / units.h,
            bounds=[1e-11, None],
            doc="Mole flux of solutes across the membrane (z-direction, x-dependent)",
        )
        self.retentate_flow_volume = Var(
            self.time,
            self.dimensionless_module_length,
            units=units.m**3 / units.h,
            bounds=[1e-11, None],
            doc="Volumetric flow rate of the retentate, x-dependent",
        )
        self.retentate_conc_mol_comp = Var(
            self.time,
            self.dimensionless_module_length,
            self.solutes,
            units=units.mol / units.m**3,  # mM
            bounds=[1e-11, None],
            doc="Mole concentration of solutes in the retentate, x-dependent",
        )
        self.permeate_flow_volume = Var(
            self.time,
            self.dimensionless_module_length,
            units=units.m**3 / units.h,
            bounds=[1e-11, None],
            doc="Volumetric flow rate of the permeate, x-dependent",
        )
        self.permeate_conc_mol_comp = Var(
            self.time,
            self.dimensionless_module_length,
            self.solutes,
            units=units.mol / units.m**3,  # mM
            bounds=[1e-11, None],
            doc="Mole concentration of solutes in the permeate, x-dependent",
        )
        self.osmotic_pressure = Var(
            self.time,
            self.dimensionless_module_length,
            units=units.bar,
            bounds=[1e-11, None],
            doc="Osmostic pressure difference across the membrane",
        )
        self.Donnan_potential_feed_side = Var(
            self.time,
            self.dimensionless_module_length,
            units=units.dimensionless,
            doc="Dimensionless Donnan potential (feed-side)",
        )
        self.Donnan_potential_permeate_side = Var(
            self.time,
            self.dimensionless_module_length,
            units=units.dimensionless,
            doc="Dimensionless Donnan potential (permeate-side)",
        )
        self.partitioning_term_bilinear_feed = Var(
            self.time,
            self.dimensionless_module_length,
            self.solutes,
            units=units.mol / units.m**3,  # mM
            bounds=[1e-11, None],
            doc="Bi-linear partitioning term for the feed-side interface",
        )
        self.partitioning_term_bilinear_permeate = Var(
            self.time,
            self.dimensionless_module_length,
            self.solutes,
            units=units.mol / units.m**3,  # mM
            bounds=[1e-11, None],
            doc="Bi-linear partitioning term for the permeate-side interface",
        )

        # add variables dependent on dimensionless_module_length and dimensionless_membrane_thickness
        if self.config.include_boundary_layer:
            self.boundary_layer_conc_mol_comp = Var(
                self.time,
                self.dimensionless_module_length,
                self.dimensionless_boundary_layer_thickness,
                self.solutes,
                units=units.mol / units.m**3,  # mM
                bounds=[1e-11, None],
                doc="Mole concentration of solutes in the boundary layer, x- and z-dependent",
            )
            self.boundary_layer_D_tilde = Var(
                self.time,
                self.dimensionless_module_length,
                self.dimensionless_boundary_layer_thickness,
                units=(units.mm**2 / units.hr) * (units.mol / units.m**3),  # D * c
                doc="Denominator of diffusion and convection coefficients in boundary layer",
            )
            self.boundary_layer_cross_diffusion_coefficient_bilinear = Var(
                self.time,
                self.dimensionless_module_length,
                self.dimensionless_boundary_layer_thickness,
                self.cations,
                self.cations,
                units=(units.mm**2 / units.h)
                * (units.mm**2 / units.h * units.mol / units.m**3),  # D * D,tilde
                doc="Bi-linear cross diffusion coefficient for cations in boundary layer",
            )
            self.boundary_layer_cross_diffusion_coefficient = Var(
                self.time,
                self.dimensionless_module_length,
                self.dimensionless_boundary_layer_thickness,
                self.cations,
                self.cations,
                units=units.mm**2 / units.h,
                doc="Cross diffusion coefficient for cations in boundary layer",
            )

        self.membrane_conc_mol_comp = Var(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            self.solutes,
            units=units.mol / units.m**3,  # mM
            bounds=[1e-11, None],
            doc="Mole concentration of solutes in the membrane, x- and z-dependent",
        )
        self.membrane_D_tilde = Var(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            units=(units.mm**2 / units.hr) * (units.mol / units.m**3),  # D * c
            doc="Denominator of diffusion and convection coefficients in membrane",
        )
        self.membrane_cross_diffusion_coefficient_bilinear = Var(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            self.cations,
            self.cations,
            units=(units.mm**2 / units.h)
            * (units.mm**2 / units.h * units.mol / units.m**3),  # D * D,tilde
            doc="Bi-linear cross diffusion coefficient for cations in membrane",
        )
        self.membrane_convection_coefficient_bilinear = Var(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            self.cations,
            units=(units.mm**2 / units.hr) * (units.mol / units.m**3),  # D,tilde
            doc="Convection coefficient for cations in membrane",
        )
        self.membrane_cross_diffusion_coefficient = Var(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            self.cations,
            self.cations,
            units=units.mm**2 / units.h,
            doc="Cross diffusion coefficient for cations in membrane",
        )
        self.membrane_convection_coefficient = Var(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            self.cations,
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
        def _differential_overall_mass_balance(blk, t, x):
            if x == 0:
                return Constraint.Skip
            return blk.d_retentate_flow_volume_dx[t, x] == (
                -blk.volume_flux_water[t, x]
                * blk.total_membrane_length
                * blk.total_module_length
            )

        self.differential_overall_mass_balance = Constraint(
            self.time,
            self.dimensionless_module_length,
            rule=_differential_overall_mass_balance,
        )

        def _differential_cation_mol_balance(blk, t, x, k):
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

        self.differential_cation_mol_balance = Constraint(
            self.time,
            self.dimensionless_module_length,
            self.cations,
            rule=_differential_cation_mol_balance,
        )

        def _overall_mass_balance(blk, t, x):
            return (
                blk.retentate_flow_volume[t, x] + blk.permeate_flow_volume[t, x]
                == blk.feed_flow_volume[t] + blk.diafiltrate_flow_volume[t]
            )

        self.overall_mass_balance = Constraint(
            self.time, self.dimensionless_module_length, rule=_overall_mass_balance
        )

        def _cation_mol_balance(blk, t, x, k):
            return (
                blk.retentate_flow_volume[t, x] * blk.retentate_conc_mol_comp[t, x, k]
            ) + (
                blk.permeate_flow_volume[t, x] * blk.permeate_conc_mol_comp[t, x, k]
            ) == (
                blk.feed_flow_volume[t] * blk.feed_conc_mol_comp[t, k]
            ) + (
                blk.diafiltrate_flow_volume[t] * blk.diafiltrate_conc_mol_comp[t, k]
            )

        self.cation_mol_balance = Constraint(
            self.time,
            self.dimensionless_module_length,
            self.cations,
            rule=_cation_mol_balance,
        )

        # transport constraints (first principles)
        def _lumped_water_flux(blk, t, x):
            if x == 0:
                return Constraint.Skip
            return blk.volume_flux_water[t, x] == (
                units.convert(
                    blk.membrane_permeability, to_units=units.m / units.h / units.bar
                )
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

            def _boundary_layer_membrane_interface(blk, t, x, j):
                if x == 0:
                    return Constraint.Skip
                return blk.membrane_conc_mol_comp[t, x, 0, j] == (
                    blk.partitioning_term_bilinear_feed[t, x, j]
                )

            self.boundary_layer_membrane_interface = Constraint(
                self.time,
                self.dimensionless_module_length,
                self.solutes,
                rule=_boundary_layer_membrane_interface,
            )
        else:

            def _retentate_membrane_interface(blk, t, x, j):
                if x == 0:
                    return Constraint.Skip
                return (
                    blk.membrane_conc_mol_comp[t, x, 0, j]
                    == blk.partitioning_term_bilinear_feed[t, x, j]
                )

            self.retentate_membrane_interface = Constraint(
                self.time,
                self.dimensionless_module_length,
                self.solutes,
                rule=_retentate_membrane_interface,
            )

        def _partitioning_term_bilinear_feed_constraint(blk, t, x, j):
            if x == 0:
                return Constraint.Skip
            charge = blk.config.property_package.charge
            conc_r = blk.retentate_conc_mol_comp

            H_nonDonnan = blk.config.property_package.non_Donnan_partition_coefficient
            if self.config.include_boundary_layer:
                conc_bl = blk.boundary_layer_conc_mol_comp
                return blk.partitioning_term_bilinear_feed[t, x, j] == (
                    conc_bl[t, x, 1, j]
                    * H_nonDonnan[j]
                    * exp(-charge[j] * blk.Donnan_potential_feed_side[t, x])
                )
            else:
                return blk.partitioning_term_bilinear_feed[t, x, j] == (
                    conc_r[t, x, j]
                    * H_nonDonnan[j]
                    * exp(-charge[j] * blk.Donnan_potential_feed_side[t, x])
                )

        self.partitioning_term_bilinear_feed_constraint = Constraint(
            self.time,
            self.dimensionless_module_length,
            self.solutes,
            rule=_partitioning_term_bilinear_feed_constraint,
        )

        def _membrane_permeate_interface(blk, t, x, j):
            if x == 0:
                return Constraint.Skip
            return blk.membrane_conc_mol_comp[t, x, 1, j] == (
                blk.partitioning_term_bilinear_permeate[t, x, j]
            )

        self.membrane_permeate_interface = Constraint(
            self.time,
            self.dimensionless_module_length,
            self.solutes,
            rule=_membrane_permeate_interface,
        )

        def _partitioning_term_bilinear_permeate_constraint(blk, t, x, j):
            if x == 0:
                return Constraint.Skip
            charge = blk.config.property_package.charge
            conc_p = blk.permeate_conc_mol_comp
            H_nonDonnan = blk.config.property_package.non_Donnan_partition_coefficient
            return blk.partitioning_term_bilinear_permeate[t, x, j] == (
                conc_p[t, x, j]
                * H_nonDonnan[j]
                * exp(-charge[j] * blk.Donnan_potential_permeate_side[t, x])
            )

        self.partitioning_term_bilinear_permeate_constraint = Constraint(
            self.time,
            self.dimensionless_module_length,
            self.solutes,
            rule=_partitioning_term_bilinear_permeate_constraint,
        )

        # boundary conditions and constraints to improve numerical stability
        if self.config.include_boundary_layer:

            def _boundary_layer_conc_mol_comp_boundary_condition(blk, t, z, j):
                return (
                    blk.boundary_layer_conc_mol_comp[t, 0, z, j]
                    == self.numerical_zero_tolerance * units.mol / units.m**3
                )

            self.boundary_layer_conc_mol_comp_boundary_condition = Constraint(
                self.time,
                self.dimensionless_boundary_layer_thickness,
                self.solutes,
                rule=_boundary_layer_conc_mol_comp_boundary_condition,
            )

        def _membrane_conc_mol_comp_boundary_condition(blk, t, z, j):
            return (
                blk.membrane_conc_mol_comp[t, 0, z, j]
                == self.numerical_zero_tolerance * units.mol / units.m**3
            )

        self.membrane_conc_mol_comp_boundary_condition = Constraint(
            self.time,
            self.dimensionless_membrane_thickness,
            self.solutes,
            rule=_membrane_conc_mol_comp_boundary_condition,
        )

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

    def add_scaling_factors(self):
        """
        Assigns scaling factors to certain variables and constraints to
        improve solver performance.
        """
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        self.scaling_factor[self.volume_flux_water] = 1e2
        self.scaling_factor[self.lumped_water_flux] = 1e3

        if self.config.include_boundary_layer:
            self.scaling_factor[self.boundary_layer_D_tilde] = 1e-2
            self.scaling_factor[
                self.boundary_layer_cross_diffusion_coefficient_bilinear
            ] = 1e-3
            self.scaling_factor[
                self.boundary_layer_cross_diffusion_coefficient_bilinear_calculation
            ] = 1e-2
            self.scaling_factor[
                self.boundary_layer_cross_diffusion_coefficient_calculation
            ] = 1e-2

        self.scaling_factor[self.membrane_D_tilde] = 1e1
        self.scaling_factor[self.membrane_cross_diffusion_coefficient_bilinear] = 1e3
        self.scaling_factor[
            self.membrane_cross_diffusion_coefficient_bilinear_calculation
        ] = 1e3
        self.scaling_factor[self.membrane_convection_coefficient_bilinear] = 1e2
        self.scaling_factor[self.membrane_cross_diffusion_coefficient] = 1e5
        self.scaling_factor[self.membrane_cross_diffusion_coefficient_calculation] = 1e5
        self.scaling_factor[self.membrane_convection_coefficient] = 1e3

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
        def _feed_ionic_strength(blk, t):
            charge = blk.config.property_package.charge
            q_feed = blk.feed_flow_volume
            q_diaf = blk.diafiltrate_flow_volume
            conc_feed = blk.feed_conc_mol_comp
            conc_diaf = blk.diafiltrate_conc_mol_comp
            return 0.5 * sum(
                (
                    (
                        (q_feed[t] * conc_feed[t, j] + q_diaf[t] * conc_diaf[t, j])
                        / (q_feed[t] + q_diaf[t])
                    )
                    * charge[j] ** 2
                )
                for j in blk.solutes
            )

        self.feed_ionic_strength = Expression(self.time, rule=_feed_ionic_strength)

        def _total_feed_flow_volume(
            blk,
            t,
        ):
            return blk.feed_flow_volume[t] + blk.diafiltrate_flow_volume[t]

        self.total_feed_flow_volume = Expression(
            self.time, rule=_total_feed_flow_volume
        )

        def _total_feed_conc_mol_comp(blk, t, j):
            q_feed = blk.feed_flow_volume
            q_diaf = blk.diafiltrate_flow_volume
            conc_feed = blk.feed_conc_mol_comp
            conc_diaf = blk.diafiltrate_conc_mol_comp
            return (q_feed[t] * conc_feed[t, j] + q_diaf[t] * conc_diaf[t, j]) / (
                q_feed[t] + q_diaf[t]
            )

        self.total_feed_conc_mol_comp = Expression(
            self.time, self.solutes, rule=_total_feed_conc_mol_comp
        )

        def _overall_partition_coefficient_feed_side(blk, t, x, j):
            if self.config.include_boundary_layer:
                return (
                    blk.membrane_conc_mol_comp[t, x, 0, j]
                    / blk.boundary_layer_conc_mol_comp[t, x, 1, j]
                )
            else:
                return (
                    blk.membrane_conc_mol_comp[t, x, 0, j]
                    / blk.retentate_conc_mol_comp[t, x, j]
                )

        self.overall_partition_coefficient_feed_side = Expression(
            self.time,
            self.dimensionless_module_length,
            self.solutes,
            rule=_overall_partition_coefficient_feed_side,
        )

        def _overall_partition_coefficient_permeate_side(blk, t, x, j):
            return (
                blk.membrane_conc_mol_comp[t, x, 1, j]
                / blk.permeate_conc_mol_comp[t, x, j]
            )

        self.overall_partition_coefficient_permeate_side = Expression(
            self.time,
            self.dimensionless_module_length,
            self.solutes,
            rule=_overall_partition_coefficient_permeate_side,
        )

        def _observed_rejection_percent(blk, t, x, j):
            return (
                1
                - (
                    blk.permeate_conc_mol_comp[t, x, j]
                    / blk.retentate_conc_mol_comp[t, x, j]
                )
            ) * 100

        self.observed_rejection_percent = Expression(
            self.time,
            self.dimensionless_module_length,
            self.solutes,
            rule=_observed_rejection_percent,
        )

        if self.config.include_boundary_layer:

            def _actual_rejection_percent(blk, t, x, j):
                return (
                    1
                    - (
                        blk.permeate_conc_mol_comp[t, x, j]
                        / blk.boundary_layer_conc_mol_comp[t, x, 1, j]
                    )
                ) * 100

            self.actual_rejection_percent = Expression(
                self.time,
                self.dimensionless_module_length,
                self.solutes,
                rule=_actual_rejection_percent,
            )

        def _observed_sieving_coefficient(blk, t, x, j):
            return (
                blk.permeate_conc_mol_comp[t, x, j]
                / blk.retentate_conc_mol_comp[t, x, j]
            )

        self.observed_sieving_coefficient = Expression(
            self.time,
            self.dimensionless_module_length,
            self.solutes,
            rule=_observed_sieving_coefficient,
        )

        if self.config.include_boundary_layer:

            def _actual_sieving_coefficient(blk, t, x, j):
                return (
                    blk.permeate_conc_mol_comp[t, x, j]
                    / blk.boundary_layer_conc_mol_comp[t, x, 1, j]
                )

            self.actual_sieving_coefficient = Expression(
                self.time,
                self.dimensionless_module_length,
                self.solutes,
                rule=_actual_sieving_coefficient,
            )

            def _boundary_layer_convective_flux(blk, t, x, z, j):
                if x == 0:
                    return Constraint.Skip
                return units.convert(
                    (
                        blk.boundary_layer_conc_mol_comp[t, x, z, j]
                        * blk.volume_flux_water[t, x]
                    ),
                    to_units=units.mol / units.m**2 / units.h,
                )

            self.boundary_layer_convective_flux = Expression(
                self.time,
                self.dimensionless_module_length,
                self.dimensionless_boundary_layer_thickness,
                self.solutes,
                rule=_boundary_layer_convective_flux,
            )

            def _boundary_layer_diffusive_flux(blk, t, x, z, j):
                if x == 0:
                    return Constraint.Skip
                return units.convert(
                    (
                        -blk.config.property_package.boundary_layer_diffusion_coefficient[
                            j
                        ]
                        * blk.d_boundary_layer_conc_mol_comp_dz[t, x, z, j]
                        / blk.total_boundary_layer_thickness
                    ),
                    to_units=units.mol / units.m**2 / units.h,
                )

            self.boundary_layer_diffusive_flux = Expression(
                self.time,
                self.dimensionless_module_length,
                self.dimensionless_boundary_layer_thickness,
                self.solutes,
                rule=_boundary_layer_diffusive_flux,
            )

            def _boundary_layer_electric_potential_gradient(blk, t, x, z):
                if x == 0:
                    return Constraint.Skip
                R = Constants.gas_constant  # J / mol / K
                T = blk.temperature
                F = Constants.faraday_constant  # C / mol
                charge = blk.config.property_package.charge
                D_bl = blk.config.property_package.boundary_layer_diffusion_coefficient
                d_conc_bl = blk.d_boundary_layer_conc_mol_comp_dz
                conc_bl = blk.boundary_layer_conc_mol_comp
                return (-R * T / F) * (
                    (
                        sum(
                            charge[j]
                            * D_bl[j]
                            * (
                                d_conc_bl[t, x, z, j]
                                / blk.total_boundary_layer_thickness
                            )
                            for j in blk.solutes
                        )
                    )
                    / (
                        sum(
                            charge[j] ** 2 * D_bl[j] * conc_bl[t, x, z, j]
                            for j in blk.solutes
                        )
                    )
                )

            self.boundary_layer_electric_potential_gradient = Expression(
                self.time,
                self.dimensionless_module_length,
                self.dimensionless_boundary_layer_thickness,
                rule=_boundary_layer_electric_potential_gradient,
            )

            def _boundary_layer_electromigrative_flux(blk, t, x, z, j):
                if x == 0:
                    return Constraint.Skip
                R = Constants.gas_constant  # J / mol / K
                T = blk.temperature
                F = Constants.faraday_constant  # C / mol
                charge = blk.config.property_package.charge
                D_bl = blk.config.property_package.boundary_layer_diffusion_coefficient
                conc_bl = blk.boundary_layer_conc_mol_comp
                return units.convert(
                    (
                        -(charge[j] * D_bl[j] * F)
                        / (R * T)
                        * conc_bl[t, x, z, j]
                        * blk.boundary_layer_electric_potential_gradient[t, x, z]
                    ),
                    to_units=units.mol / units.m**2 / units.h,
                )

            self.boundary_layer_electromigrative_flux = Expression(
                self.time,
                self.dimensionless_module_length,
                self.dimensionless_boundary_layer_thickness,
                self.solutes,
                rule=_boundary_layer_electromigrative_flux,
            )

        def _membrane_convective_flux(blk, t, x, z, j):
            if x == 0:
                return Constraint.Skip
            return units.convert(
                (blk.membrane_conc_mol_comp[t, x, z, j] * blk.volume_flux_water[t, x]),
                to_units=units.mol / units.m**2 / units.h,
            )

        self.membrane_convective_flux = Expression(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            self.solutes,
            rule=_membrane_convective_flux,
        )

        def _membrane_diffusive_flux(blk, t, x, z, j):
            if x == 0:
                return Constraint.Skip
            return units.convert(
                (
                    -blk.config.property_package.membrane_diffusion_coefficient[j]
                    * blk.d_membrane_conc_mol_comp_dz[t, x, z, j]
                    / blk.total_membrane_thickness
                ),
                to_units=units.mol / units.m**2 / units.h,
            )

        self.membrane_diffusive_flux = Expression(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            self.solutes,
            rule=_membrane_diffusive_flux,
        )

        def _membrane_electric_potential_gradient(blk, t, x, z):
            if x == 0:
                return Constraint.Skip
            R = Constants.gas_constant  # J / mol / K
            T = blk.temperature
            F = Constants.faraday_constant  # C / mol
            charge = blk.config.property_package.charge
            D_mem = blk.config.property_package.membrane_diffusion_coefficient
            d_con_mem = blk.d_membrane_conc_mol_comp_dz
            con_mem = blk.membrane_conc_mol_comp
            J_w = blk.volume_flux_water
            chi = blk.membrane_fixed_charge
            return (-R * T / F) * (
                (
                    units.convert(
                        sum(
                            charge[j]
                            * D_mem[j]
                            * (d_con_mem[t, x, z, j] / blk.total_membrane_thickness)
                            for j in blk.solutes
                        ),
                        to_units=units.mol / units.m**2 / units.h,
                    )
                    + (J_w[t, x] * chi)
                )
                / (
                    sum(
                        charge[j] ** 2 * D_mem[j] * con_mem[t, x, z, j]
                        for j in blk.solutes
                    )
                )
            )

        self.membrane_electric_potential_gradient = Expression(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            rule=_membrane_electric_potential_gradient,
        )

        def _membrane_electromigrative_flux(blk, t, x, z, j):
            if x == 0:
                return Constraint.Skip
            R = Constants.gas_constant  # J / mol / K
            T = blk.temperature
            F = Constants.faraday_constant  # C / mol
            charge = blk.config.property_package.charge
            D_mem = blk.config.property_package.membrane_diffusion_coefficient
            conc_mem = blk.membrane_conc_mol_comp
            return units.convert(
                (
                    -(charge[j] * D_mem[j] * F)
                    / (R * T)
                    * conc_mem[t, x, z, j]
                    * blk.membrane_electric_potential_gradient[t, x, z]
                ),
                to_units=units.mol / units.m**2 / units.h,
            )

        self.membrane_electromigrative_flux = Expression(
            self.time,
            self.dimensionless_module_length,
            self.dimensionless_membrane_thickness,
            self.solutes,
            rule=_membrane_electromigrative_flux,
        )
