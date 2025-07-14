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

This membrane unit model is for the multi-component diafiltration of a two-salt system with a common anion. The membrane is designed for use in a diafiltration cascade, i.e., represents a spiral-wound membrane module.

Configuration Arguments
-----------------------

The Two-Salt Diafiltration model requires a property package that provides the valency (:math:`z_i`), reflection coefficient (:math:`\sigma_i`), partition coefficient (:math:`H_i`), and number of dissolved species (:math:`n_i`) for each ion :math:`i` in solution. When used in a flowsheet, the user can provide separate property packages for the feed and product streams. 

Additionally, there are two required arguments, ``NFEx`` and ``NFEz``, to specfiy the desired number of finite elements across the width and thickness of the membrane, respectively.

Degrees of Freedom
------------------

The Two-Salt Diafiltration unit model has three degrees of freedom (variable names and default values are provided in parentheses):

#. the width of the membrane, i.e., the length of the membrane module (``membrane_width``; :math:`4 \, \mathrm{m}`)
#. the length of the membrane (``membrane_length``; :math:`40 \, \mathrm{m}`)
#. the pressure applied to the membrane system (``applied_pressure``; :math:`15 \, \mathrm{bar}`)

To run a simulation (with zero degrees of freedom) in a flowsheet, the additional following variables must be fixed to obtain zero degrees of freedom (variable names and default values are provided in parentheses):

#. the volumetric flow rate of the feed (``feed_flow_volume``; :math:`100 \, \mathrm{m}^3 \, \mathrm{h}^{-1}`)
#. the lithium concentration in the feed (``feed_conc_mol_comp[t,"Li"]``; :math:`245 \, \mathrm{mol} \, \mathrm{m}^{-3}`)
#. the cobalt concentration in the feed (``feed_conc_mol_comp[t,"Co"]``; :math:`288 \, \mathrm{mol} \, \mathrm{m}^{-3}`)
#. the volumetric flow rate of the diafiltrate (``diafiltrate_flow_volume``; :math:`30 \, \mathrm{m}^3 \, \mathrm{h}^{-1}`)
#. the lithium concentration in the diafiltrate (``diafiltrate_conc_mol_comp[0,"Li"]``; :math:`14 \, \mathrm{mol} \, \mathrm{m}^{-3}`)
#. the cobalt concentration in the diafiltrate (``diafiltrate_conc_mol_comp[0,"Co"]``; :math:`3 \, \mathrm{mol} \, \mathrm{m}^{-3}`)

Model Structure
---------------

There are three phases in the Two-Salt Diafiltration model: the retentate, the membrane, and the permeate. The retentate and the permeate are only discretized with respect to :math:`x`, while the membrane is discretized with respect to both :math:`x` and :math:`z`. The resulting system of partial differential algebraic equations is solved by discretizing with the backward finite element method.

Assumptions
-----------

The partition coefficients, which describe how the solutes partition across the solution-membrane interfaces, are derived considering steric and electrostatic hindrance.

The membrane has a negative surface charge (:math:`-140 \, \mathrm{mM}`), which was calculated using zeta potential measurements for NF270 membranes. (See `this reference <https://doi.org/10.1021/acs.iecr.4c04763>`_).

The formation of a boundary layer at the membrane surface due to concentration polarization is neglected for mathematical simplicity.

The dominating transport mechanism within the bulk/retentate solution is convection in the :math:`x`-direction (parallel to the membrane surface). The dominating transport mechanism within the permeate solution is convection in the :math:`z`-direction (perpendicular to the membrane surface).

The transport mechanisms modeled within the membrane are convection, diffusion, and electromigration. Diffusion within the membrane that is normal to the pore walls is ignored, meaning the concentration gradient of ion :math:`i` within the membrane only has a :math:`z`-component (perpendicular to the membrane surface).

Sets
----

The Two-Salt Diafiltration model defines the following discrete sets for ions in the system.

.. math:: \mathcal{I}=\{\mathrm{Li^+,Co^{2+},Cl^-}\}

There are 2 continuous sets for each length dimension: the :math:`x`-direction parallel to the membrane surface and the :math:`z`-direction perpendicular to the membrane surface. :math:`x` and :math:`z` are non-dimensionalized with the module length or membrane width (:math:`w`) and membrane thickness (:math:`l`), respectively, to improve numerics.

.. math:: \bar{x} \in \mathbb{R} \| 0 \leq \bar{x} \leq 1
.. math:: \bar{z} \in \mathbb{R} \| 0 \leq \bar{z} \leq 1

Some variables are discretized over time to be compatible with the property package, even though this is not a dynamic model. Thus, the following set is defined for time.

.. math:: t \in [0]

Default Model Parameters
------------------------

The Two-Salt Diafiltration model has the following parameters.

================ =============================================== ============================ ============= ==========================================================
Parameter        Description                                     Name                         Default Value Units
================ =============================================== ============================ ============= ==========================================================
:math:`\epsilon` numerical tolerance for zero values             ``numerical_zero_tolerance`` 1e-10
:math:`l`        thickness of the membrane                       ``membrane_thickness``       1e-07         :math:`\mathrm{m}`
:math:`L_p`      hydraulic permeability of the membrane          ``membrane_permeability``    0.03          :math:`\mathrm{m} \, \mathrm{h}^{-1} \, \mathrm{bar}^{-1}`
:math:`T`        temperature of the system                       ``temperature``              298           :math:`\mathrm{K}`
:math:`\chi`     concentration of surface charge on the membrane ``membrane_fixed_charge``    -140          :math:`\mathrm{mol} \, \mathrm{m}^{-3}`
================ =============================================== ============================ ============= ==========================================================

Variables
---------

The Two-Salt Diafiltration model adds the following variables.

==================================== =============================================== ====================================================== ========================================================== ====================================================
Variable                             Description                                     Name                                                   Units                                                      Notes
==================================== =============================================== ====================================================== ========================================================== ====================================================
:math:`c_{\mathrm{Co^{2+}},d}`       concentration of cobalt ion in the diafiltrate  ``diafiltrate_conc_mol_comp[t, "Co"]``                 :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                    discretized over :math:`t`
:math:`c_{\mathrm{Cl^-},d}`          concentration of chloride in the diafiltrate    ``diafiltrate_conc_mol_comp[t, "Cl"]``                 :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                    discretized over :math:`t`
:math:`c_{\mathrm{Li^+},d}`          concentration of lithium ion in the diafiltrate ``diafiltrate_conc_mol_comp[t, "Li"]``                 :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                    discretized over :math:`t`
:math:`c_{\mathrm{Co^{2+}},f}`       concentration of cobalt ion in the feed         ``feed_conc_mol_comp[t, "Co"]``                        :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                    discretized over :math:`t`
:math:`c_{\mathrm{Cl^-},f}`          concentration of chloride in the feed           ``feed_conc_mol_comp[t, "Cl"]``                        :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                    discretized over :math:`t`
:math:`c_{\mathrm{Li^+},f}`          concentration of lithium ion in the feed        ``feed_conc_mol_comp[t, "Li"]``                        :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                    discretized over :math:`t`
:math:`c_{\mathrm{Cl^-},m}`          concentration of chloride in the membrane       ``membrane_conc_mol_chloride``                         :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                    discretized over :math:`\bar{x}` and :math:`\bar{z}`
:math:`c_{\mathrm{Co^{2+}},m}`       concentration of cobalt ion in the membrane     ``membrane_conc_mol_cobalt``                           :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                    discretized over :math:`\bar{x}` and :math:`\bar{z}`
:math:`c_{\mathrm{Li^+},m}`          concentration of lithium ion in the membrane    ``membrane_conc_mol_lithium``                          :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                    discretized over :math:`\bar{x}` and :math:`\bar{z}`
:math:`c_{\mathrm{Cl^-},p}`          concentration of chloride in the membrane       ``permeate_conc_mol_comp[t, "Cl", :math:`\bar{x}`]``   :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                    discretized over :math:`t` and :math:`\bar{x}`
:math:`c_{\mathrm{Co^{2+}},p}`       concentration of cobalt ion in the permeate     ``permeate_conc_mol_comp[t, "Co", :math:`\bar{x}`]``   :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                    discretized over :math:`t` and :math:`\bar{x}`
:math:`c_{\mathrm{Li^+},p}`          concentration of lithium ion in the permeate    ``permeate_conc_mol_comp[t, "Li", :math:`\bar{x}`]``   :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                    discretized over :math:`t` and :math:`\bar{x}`
:math:`c_{\mathrm{Cl^-},r}`          concentration of chloride in the retentate      ``retentate_conc_mol_comp[t, "Cl", :math:`\bar{x}`]``  :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                    discretized over :math:`t` and :math:`\bar{x}`
:math:`c_{\mathrm{Co^{2+}},r}`       concentration of cobalt ion in the retentate    ``retentate_conc_mol_comp[t, "Co", :math:`\bar{x}`]``  :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                    discretized over :math:`t` and :math:`\bar{x}`
:math:`c_{\mathrm{Li^+},r}`          concentration of lithium ion in the retentate   ``retentate_conc_mol_comp[t, "Li", :math:`\bar{x}`]``  :math:`\mathrm{mol} \, \mathrm{m}^{-3}`                    discretized over :math:`t` and :math:`\bar{x}`
:math:`\alpha_{Li^+}`                linearized convection coefficient (Li)          ``convection_coefficient_lithium``                     :math:`\mathrm{dimensionless}`                             discretized over :math:`\bar{x}` and :math:`\bar{z}`
:math:`\alpha_{Co^{2+}}`             linearized convection coefficient (Co)          ``convection_coefficient_cobalt``                      :math:`\mathrm{dimensionless}`                             discretized over :math:`\bar{x}` and :math:`\bar{z}`
:math:`D_{\mathrm{Li^+,Li^+}}`       linearized cross-diffusion coefficient (Li,Li)  ``D_lithium_lithium``                                  :math:`\mathrm{m}^2 \, \mathrm{h}^{-1}`                    discretized over :math:`\bar{x}` and :math:`\bar{z}`
:math:`D_{\mathrm{Li^+,Co^{2+}}}`    linearized cross-diffusion coefficient (Li,Co)  ``D_lithium_cobalt``                                   :math:`\mathrm{m}^2 \, \mathrm{h}^{-1}`                    discretized over :math:`\bar{x}` and :math:`\bar{z}`
:math:`D_{\mathrm{Co^{2+},Li^+}}`    linearized cross-diffusion coefficient (Co,Li)  ``D_cobalt_lithium``                                   :math:`\mathrm{m}^2 \, \mathrm{h}^{-1}`                    discretized over :math:`\bar{x}` and :math:`\bar{z}`
:math:`D_{\mathrm{Co^{2+},Co^{2+}}}` linearized cross-diffusion coefficient (Co,Co)  ``D_cobalt_cobalt``                                    :math:`\mathrm{m}^2 \, \mathrm{h}^{-1}`                    discretized over :math:`\bar{x}` and :math:`\bar{z}`
:math:`j_{\mathrm{Cl^-}}`            molar flux of chloride across the membrane      ``mol_flux_chloride``                                  :math:`\mathrm{mol} \, \mathrm{m}^{-2} \, \mathrm{h}^{-1}` discretized over :math:`\bar{x}`
:math:`j_{\mathrm{Co^{2+}}}`         molar flux of cobalt ion across the membrane    ``mol_flux_cobalt``                                    :math:`\mathrm{mol} \, \mathrm{m}^{-2} \, \mathrm{h}^{-1}` discretized over :math:`\bar{x}`
:math:`j_{\mathrm{Li^+}}`            molar flux of lithium ion across the membrane   ``mol_flux_lithium``                                   :math:`\mathrm{mol} \, \mathrm{m}^{-2} \, \mathrm{h}^{-1}` discretized over :math:`\bar{x}`
:math:`J_w`                          water flux across the membrane                  ``volume_flux_water``                                  :math:`\mathrm{m}^3 \, \mathrm{m}^{-2} \, \mathrm{h}^{-1}` discretized over :math:`\bar{x}`
:math:`L`                            length of the membrane                          ``membrane_length``                                    :math:`\mathrm{m}`
:math:`\Delta \pi`                   osmotic pressure of feed-side fluid             ``osmotic_pressure``                                   :math:`\mathrm{bar}`                                       discretized over :math:`\bar{x}`
:math:`\Delta P`                     applied pressure to the membrane                ``applied_pressure``                                   :math:`\mathrm{bar}`
:math:`q_d`                          volumetic flow rate of the diafiltrate          ``diafiltrate_flow_volume``                            :math:`\mathrm{m}^3 \, \mathrm{h}^{-1}`                    discretized over :math:`t`
:math:`q_f`                          volumetic flow rate of the feed                 ``feed_flow_volume``                                   :math:`\mathrm{m}^3 \, \mathrm{h}^{-1}`                    discretized over :math:`t`
:math:`q_p`                          volumetic flow rate of the permeate             ``permeate_flow_volume``                               :math:`\mathrm{m}^3 \, \mathrm{h}^{-1}`                    discretized over :math:`t` and :math:`\bar{x}`
:math:`q_r`                          volumetic flow rate of the retentate            ``retentate_flow_volume``                              :math:`\mathrm{m}^3 \, \mathrm{h}^{-1}`                    discretized over :math:`t` and :math:`\bar{x}`
:math:`w`                            width of the membrane                           ``membrane_width``                                     :math:`\mathrm{m}`
==================================== =============================================== ====================================================== ========================================================== ====================================================

Derivative Variables
--------------------

The Two-Salt Diafiltration model adds the following derivative variables.

================================================================== =============================================== =========================================================== ======================================= ====================================================
Variable                                                           Description                                     Name                                                        Units                                   Notes
================================================================== =============================================== =========================================================== ======================================= ====================================================
:math:`\frac{\mathrm{d}c_{\mathrm{Co^{2+}},r}}{\mathrm{d}\bar{x}}` cobalt ion concentration gradient in retentate  ``d_retentate_conc_mass_comp_dx[t, "Co", :math:`\bar{x}`]`` :math:`\mathrm{kg} \, \mathrm{m}^{-3}`  discretized over :math:`t` and :math:`\bar{x}`
:math:`\frac{\mathrm{d}c_{\mathrm{Li^+},r}}{\mathrm{d}\bar{x}}`    lithium ion concentration gradient in retentate ``d_retentate_conc_mass_comp_dx[t, "Li", :math:`\bar{x}`]`` :math:`\mathrm{kg} \, \mathrm{m}^{-3}`  discretized over :math:`t` and :math:`\bar{x}`
:math:`\frac{\mathrm{d}q_r}{\mathrm{d}\bar{x}}`                    retentate flow rate gradient                    ``d_retentate_flow_volume_dx``                              :math:`\mathrm{m}^3 \, \mathrm{h}^{-1}` discretized over :math:`t` and :math:`\bar{x}`
:math:`\frac{\partial c_{\mathrm{Co^{2+}},m}}{\partial \bar{z}}`   cobalt ion concentration gradient in membrane   ``d_membrane_conc_mass_cobalt_dz``                          :math:`\mathrm{kg} \, \mathrm{m}^{-3}`  discretized over :math:`\bar{x}` and :math:`\bar{z}`
:math:`\frac{\partial c_{\mathrm{Li^+},m}}{\partial \bar{z}}`      lithium ion concentration gradient in membrane  ``d_membrane_conc_mass_lithium_dz``                         :math:`\mathrm{kg} \, \mathrm{m}^{-3}`  discretized over :math:`\bar{x}` and :math:`\bar{z}`
================================================================== =============================================== =========================================================== ======================================= ====================================================

Constraints
-----------

Differential mass balances in the retentate:

.. math:: \frac{\mathrm{d}q_r(\bar{x})}{\mathrm{d}\bar{x}} = - J_w(\bar{x}) wL  \qquad \forall \, \bar{x} \neq 0
.. math:: q_r(\bar{x}) \frac{\mathrm{d}c_{\mathrm{Li^+},r}(\bar{x})}{\mathrm{d}\bar{x}} = wL (J_w(\bar{x}) c_{\mathrm{Li^+},r}(\bar{x}) - j_{\mathrm{Li^+}}(\bar{x}))  \qquad \forall \, \bar{x} \neq 0
.. math:: q_r(\bar{x}) \frac{\mathrm{d}c_{\mathrm{Co^{2+}},r}(\bar{x})}{\mathrm{d}\bar{x}} = wL (J_w(\bar{x}) c_{\mathrm{Co^{2+}},r}(\bar{x}) - j_{\mathrm{Co^{2+}}}(\bar{x}))  \qquad \forall \, \bar{x} \neq 0

Electroneutrality in the retentate:

.. math:: 0 = z_{\mathrm{Li^+}} c_{\mathrm{Li^+},r}(\bar{x}) + z_{\mathrm{Co^{2+}}} c_{\mathrm{Co^{2+}},r}(\bar{x}) + z_{\mathrm{Cl^-}} c_{\mathrm{Cl^-},r}(\bar{x})

Overall water flux through the membrane:

.. math:: J_w (\bar{x}) = L_p (\Delta P - \Delta \pi (\bar{x})) \qquad \forall \, \bar{x} \neq 0

Osmotic pressure:

.. math:: \Delta \pi (\bar{x}) = \mathrm{R} \mathrm{T} \sum_{i \in \mathcal{I}} n_i \sigma_i (c_{i,r}(\bar{x})-c_{i,p}(\bar{x}))

Nernst-Plank equations for the ion flux through the membrane:

.. math:: j_{\mathrm{Li^+}}(\bar{x}) = \alpha_{Li^+}(\bar{x},\bar{z}) c_{\mathrm{Li^+},m}(\bar{x},\bar{z}) J_w(\bar{x}) - \frac{D_{\mathrm{Li^+,Li^+}}}{l} \frac{\partial c_{\mathrm{Li^+},m}(\bar{x},\bar{z})}{\partial \bar{z}} - \frac{D_{\mathrm{Li^+,Co^{2+}}}}{l} \frac{\partial c_{\mathrm{Co^{2+}},m}(\bar{x},\bar{z})}{\partial \bar{z}} \qquad \forall \, \bar{z} \neq 0
.. math:: j_{\mathrm{Co^{2+}}}(\bar{x}) = \alpha_{Li^+}(\bar{x},\bar{z}) c_{\mathrm{Co^{2+}},m}(\bar{x},\bar{z}) J_w(\bar{x}) - \frac{D_{\mathrm{Co^{2+},Li^+}}}{l} \frac{\partial c_{\mathrm{Li^+},m}(\bar{x},\bar{z})}{\partial \bar{z}} - \frac{D_{\mathrm{Co^{2+},Co^{2+}}}}{l} \frac{\partial c_{\mathrm{Co^{2+}},m}(\bar{x},\bar{z})}{\partial \bar{z}} \qquad \forall \, \bar{z} \neq 0

with linearized convection and cross-diffusion coefficients:

.. math:: \alpha_{Li^+}(\bar{x},\bar{z}) = \omega_0 + \omega_1 c_{\mathrm{Li^+},m}(\bar{x},\bar{z}) + \omega_2 c_{\mathrm{Co^{2+}},m}(\bar{x},\bar{z})
.. math:: \alpha_{Co^{2+}}(\bar{x},\bar{z}) = \omega_3 + \omega_4 c_{\mathrm{Li^+},m}(\bar{x},\bar{z}) + \omega_5 c_{\mathrm{Co^{2+}},m}(\bar{x},\bar{z})

.. math:: - D_{\mathrm{Li^+,Li^+}}(\bar{x},\bar{z}) = \beta_0 + \beta_1 c_{\mathrm{Li^+},m}(\bar{x},\bar{z}) + \beta_2 c_{\mathrm{Co^{2+}},m}(\bar{x},\bar{z})
.. math:: - D_{\mathrm{Li^+,Co^{2+}}}(\bar{x},\bar{z}) = \beta_3 + \beta_4 c_{\mathrm{Li^+},m}(\bar{x},\bar{z}) + \beta_5 c_{\mathrm{Co^{2+}},m}(\bar{x},\bar{z})
.. math:: - D_{\mathrm{Co^{2+},Li^+}}(\bar{x},\bar{z}) = \beta_6 + \beta_7 c_{\mathrm{Li^+},m}(\bar{x},\bar{z}) + \beta_8 c_{\mathrm{Co^{2+}},m}(\bar{x},\bar{z})
.. math:: - D_{\mathrm{Co^{2+},Co^{2+}}}(\bar{x},\bar{z}) = \beta_9 + \beta_{10} c_{\mathrm{Li^+},m}(\bar{x},\bar{z}) + \beta_{11} c_{\mathrm{Co^{2+}},m}(\bar{x},\bar{z})

that have the following regressed parameter values when :math:`\chi = 0 mM` (valid for membrane lithium and cobalt concentration range of :math:`50-80 mM`) and when :math:`\chi = -140 mM` (valid for membrane lithium concentration range of :math:`50-80 mM` and membrane cobalt concentration range of :math:`80-110 mM`):

================== =========================== ============================== ============================================================
Parameter          Value (:math:`\chi = 0 mM`) Value (:math:`\chi = -140 mM`) Units
================== =========================== ============================== ============================================================
:math:`\omega_0`   :math:`1`                   :math:`0.260`                  :math:`\mathrm{dimensionless}`
:math:`\omega_1`   :math:`0`                   :math:`0.00139`                :math:`\mathrm{mol} \, \mathrm{m}^{-3}`
:math:`\omega_2`   :math:`0`                   :math:`0.00314`                :math:`\mathrm{mol} \, \mathrm{m}^{-3}`
:math:`\omega_3`   :math:`1`                   :math:`0.0860`                 :math:`\mathrm{dimensionless}`
:math:`\omega_4`   :math:`0`                   :math:`0.00198`                :math:`\mathrm{mol} \, \mathrm{m}^{-3}`
:math:`\omega_5`   :math:`0`                   :math:`0.00448`                :math:`\mathrm{mol} \, \mathrm{m}^{-3}`
:math:`\beta_0`    :math:`-4.07e-06`           :math:`-4.33e-06`              :math:`\mathrm{m}^2 \, \mathrm{h}^{-1}`
:math:`\beta_1`    :math:`-3.99e-09`           :math:`-4.25e-09`              :math:`\mathrm{m}^5 \, \mathrm{mol}^{-1} \, \mathrm{h}^{-1}`
:math:`\beta_2`    :math:`4.01e-09`            :math:`5.14e-09`               :math:`\mathrm{m}^5 \, \mathrm{mol}^{-1} \, \mathrm{h}^{-1}`
:math:`\beta_3`    :math:`-9.62e-07`           :math:`-1.63e-06`              :math:`\mathrm{m}^2 \, \mathrm{h}^{-1}`
:math:`\beta_4`    :math:`-1.03e-08`           :math:`-1.10e-08`              :math:`\mathrm{m}^5 \, \mathrm{mol}^{-1} \, \mathrm{h}^{-1}`
:math:`\beta_5`    :math:`1.04e-08`            :math:`1.33e-08`               :math:`\mathrm{m}^5 \, \mathrm{mol}^{-1} \, \mathrm{h}^{-1}`
:math:`\beta_6`    :math:`-5.24e-07`           :math:`-1.32e-06`              :math:`\mathrm{m}^2 \, \mathrm{h}^{-1}`
:math:`\beta_7`    :math:`2.49e-09`            :math:`4.72e-09`               :math:`\mathrm{m}^5 \, \mathrm{mol}^{-1} \, \mathrm{h}^{-1}`
:math:`\beta_8`    :math:`-2.50e-09`           :math:`1.47e-09`               :math:`\mathrm{m}^5 \, \mathrm{mol}^{-1} \, \mathrm{h}^{-1}`
:math:`\beta_9`    :math:`-4.00e-06`           :math:`-6.05e-06`              :math:`\mathrm{m}^2 \, \mathrm{h}^{-1}`
:math:`\beta_{10}` :math:`6.45e-09`            :math:`1.22e-08`               :math:`\mathrm{m}^5 \, \mathrm{mol}^{-1} \, \mathrm{h}^{-1}`
:math:`\beta_{11}` :math:`-6.48e-09`           :math:`3.81e-09`               :math:`\mathrm{m}^5 \, \mathrm{mol}^{-1} \, \mathrm{h}^{-1}`
================== =========================== ============================== ============================================================

No applied potential on the system:

.. math:: 0 = z_{\mathrm{Li^+}} j_{\mathrm{Li^+}}(\bar{x}) + z_{\mathrm{Co^{2+}}} j_{\mathrm{Co^{2+}}}(\bar{x}) + z_{\mathrm{Cl^-}} j_{\mathrm{Cl^-}}(\bar{x}) \qquad \forall \, \bar{x} \neq 0

Electroneutrality in the membrane:

.. math:: 0 = z_{\mathrm{Li^+}} c_{\mathrm{Li^+},m}(\bar{x},\bar{z}) + z_{\mathrm{Co^{2+}}} c_{\mathrm{Co^{2+}},m}(\bar{x},\bar{z}) + z_{\mathrm{Cl^-}} c_{\mathrm{Cl^-},m}(\bar{x},\bar{z}) + \chi \qquad \forall \, \bar{z} \neq 0

Mass balance (via convection) on the permeate outlet:

.. math:: q_p(\bar{x}) = \bar{x} wL J_w(\bar{x}) \qquad \forall \, \bar{x} \neq 0
.. math:: j_{\mathrm{Li^+}}(\bar{x}) = c_{\mathrm{Li^+},p}(\bar{x}) J_w(\bar{x}) \qquad \forall \, \bar{x} \neq 0
.. math:: j_{\mathrm{Co^{2+}}}(\bar{x}) = c_{\mathrm{Co^{2+}},p}(\bar{x}) J_w(\bar{x}) \qquad \forall \, \bar{x} \neq 0

Electroneutrality in the retentate:

.. math:: 0 = z_{\mathrm{Li^+}} c_{\mathrm{Li^+},p}(\bar{x}) + z_{\mathrm{Co^{2+}}} c_{\mathrm{Co^{2+}},p}(\bar{x}) + z_{\mathrm{Cl^-}} c_{\mathrm{Cl^-},p}(\bar{x})

Partitioning at the retentate-membrane interface:

.. math:: H_{Li^+} H_{Cl^-} c_{\mathrm{Li^+},r}(\bar{x}) c_{\mathrm{Cl^-},r}(\bar{x}) = c_{\mathrm{Li^+},m}(\bar{x},\bar{z}=0) c_{\mathrm{Cl^-},m}(\bar{x},\bar{z}=0) \qquad \forall \, \bar{x} \neq 0
.. math:: H_{Co^{2+}} H_{Cl^-}^{z_{Co^{2+}}} c_{\mathrm{Co^{2+}},r}(\bar{x}) c_{\mathrm{Cl^-},r}(\bar{x})^{z_{Co^{2+}}} =c_{\mathrm{Co^{2+}},m}(\bar{x},\bar{z}=0) c_{\mathrm{Cl^-},m}(\bar{x},\bar{z}=0)^{z_{Co^{2+}}} \qquad \forall \, \bar{x} \neq 0

Partitioning at the membrane-permeate interface:

.. math:: H_{Li^+} H_{Cl^-} c_{\mathrm{Li^+},p}(\bar{x}) c_{\mathrm{Cl^-},p}(\bar{x}) = c_{\mathrm{Li^+},m}(\bar{x},\bar{z}=1) c_{\mathrm{Cl^-},m}(\bar{x},\bar{z}=1) \qquad \forall \, \bar{x} \neq 0
.. math:: H_{Co^{2+}} H_{Cl^-}^{z_{Co^{2+}}} c_{\mathrm{Co^{2+}},p}(\bar{x}) c_{\mathrm{Cl^-},p}(\bar{x})^{z_{Co^{2+}}} =c_{\mathrm{Co^{2+}},m}(\bar{x},\bar{z}=1) c_{\mathrm{Cl^-},m}(\bar{x},\bar{z}=1)^{z_{Co^{2+}}} \qquad \forall \, \bar{x} \neq 0

The following initial conditions are fixed to complete the model:

.. math:: q_r(\bar{x}=0) = q_f + q_d
.. math:: c_{\mathrm{Li^+},r}(\bar{x}=0) = \frac{q_f c_{\mathrm{Li^+},f} + q_d c_{\mathrm{Li^+},d}}{q_f + q_d}
.. math:: c_{\mathrm{Co^{2+}},r}(\bar{x}=0) = \frac{q_f c_{\mathrm{Co^{2+}},f} + q_d c_{\mathrm{Co^{2+}},d}}{q_f + q_d}
.. math:: q_p(\bar{x}=0) = \epsilon
.. math:: c_{\mathrm{Li^+},p}(\bar{x}=0) = \epsilon
.. math:: c_{\mathrm{Co^{2+}},p}(\bar{x}=0) = \epsilon
.. math:: c_{\mathrm{Li^+},m} (\bar{x}=0,\bar{z}) = \epsilon \qquad \forall \, \bar{z}
.. math:: c_{\mathrm{Co^{2+}},m} (\bar{x}=0,\bar{z}) = \epsilon \qquad \forall \, \bar{z}
.. math:: c_{\mathrm{Cl^-},m} (\bar{x}=0,\bar{z}) = \epsilon \qquad \forall \, \bar{z}
.. math:: \frac{\mathrm{d}q_r(\bar{x})}{\mathrm{d}\bar{x}}(\bar{x}=0)=\epsilon
.. math:: \frac{\mathrm{d}c_{\mathrm{Li^+},r}(\bar{x})}{\mathrm{d}\bar{x}}(\bar{x}=0)=\epsilon
.. math:: \frac{\mathrm{d}c_{\mathrm{Co^{2+}},r}(\bar{x})}{\mathrm{d}\bar{x}}(\bar{x}=0)=\epsilon

The following initial conditions are fixed to improve numerical stability (with the appropriate constraints deactivated as described above):

.. math:: J_w(\bar{x}=0) = \epsilon
.. math:: j_{\mathrm{Li^+}}(\bar{x}=0) = \epsilon
.. math:: j_{\mathrm{Co^{2+}}}(\bar{x}=0) = \epsilon
.. math:: j_{\mathrm{Cl^-}}(\bar{x}=0) = \epsilon

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
from idaes.core.util.scaling import constraint_autoscale_large_jac


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
        "NFEx",
        ConfigValue(
            doc="Number of discretization points in the x-direction",
        ),
    )
    CONFIG.declare(
        "NFEz",
        ConfigValue(
            doc="Number of discretization points in the z-direction",
        ),
    )
    CONFIG.declare(
        "charged_membrane",
        ConfigValue(
            default=True, doc="Boolean argument if the membrane has a surface charge"
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
        self.fix_initial_values()
        self.add_scaling_factors()
        self.add_ports()

    def add_mutable_parameters(self):
        """
        Adds parameters for the two salt diafiltration unit model.

        Assigns default values that can be changed by the user during implementation.
        """
        self.numerical_zero_tolerance = Param(
            initialize=1e-10,
            mutable=True,
            doc="Numerical tolerance for zero values in the model",
        )
        self.membrane_thickness = Param(
            initialize=1e-7,
            mutable=True,
            units=units.m,
            doc="Thickness of membrane (z-direction)",
        )
        if self.config.charged_membrane:
            self.membrane_fixed_charge = Param(
                initialize=-140,
                mutable=True,
                units=units.mol / units.m**3,  # mM
                doc="Fixed charge on the membrane",
            )
        else:
            self.membrane_fixed_charge = Param(
                initialize=0,
                mutable=True,
                units=units.mol / units.m**3,  # mM
                doc="Fixed charge on the membrane",
            )
        self.membrane_permeability = Param(
            initialize=0.03,
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
        """
        # define length scales
        self.x_bar = ContinuousSet(bounds=(0, 1))
        self.z_bar = ContinuousSet(bounds=(0, 1))

        # add a time index since the property package variables are indexed over time
        self.time = Set(initialize=[0])

        # add components
        self.solutes = Set(initialize=["Li", "Co", "Cl"])

        # add global variables
        self.membrane_width = Var(
            initialize=4,
            units=units.m,
            bounds=[1e-11, None],
            doc="Width of the membrane (x-direction)",
        )
        self.membrane_length = Var(
            initialize=40,
            units=units.m,
            bounds=[1e-11, None],
            doc="Length of the membrane, wound radially",
        )
        self.applied_pressure = Var(
            initialize=15,
            units=units.bar,
            bounds=[1e-11, None],
            doc="Pressure applied to membrane",
        )
        self.feed_flow_volume = Var(
            self.time,
            initialize=100,
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
            initialize=30,
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

        # add variables dependent on x_bar
        self.volume_flux_water = Var(
            self.x_bar,
            initialize=0.3,
            units=units.m**3 / units.m**2 / units.h,
            bounds=[1e-11, None],
            doc="Volumetric water flux of water across the membrane",
        )
        self.mol_flux_lithium = Var(
            self.x_bar,
            initialize=50,
            units=units.mol / units.m**2 / units.h,
            bounds=[1e-11, None],
            doc="Mole flux of lithium across the membrane (z-direction, x-dependent)",
        )
        self.mol_flux_cobalt = Var(
            self.x_bar,
            initialize=60,
            units=units.mol / units.m**2 / units.h,
            bounds=[1e-11, None],
            doc="Mole flux of cobalt across the membrane (z-direction, x-dependent)",
        )
        self.mol_flux_chloride = Var(
            self.x_bar,
            initialize=170,
            units=units.mol / units.m**2 / units.h,
            bounds=[1e-11, None],
            doc="Mole flux of chloride across the membrane (z-direction, x-dependent)",
        )
        self.retentate_flow_volume = Var(
            self.time,
            self.x_bar,
            initialize=100,
            units=units.m**3 / units.h,
            bounds=[1e-11, None],
            doc="Volumetric flow rate of the retentate, x-dependent",
        )

        def initialize_retentate_conc_mol_comp(m, t, p, j):
            vals = {"Li": 48, "Co": 222, "Cl": 492}
            return vals[j]

        self.retentate_conc_mol_comp = Var(
            self.time,
            self.x_bar,
            self.solutes,
            initialize=initialize_retentate_conc_mol_comp,
            units=units.mol / units.m**3,  # mM
            bounds=[1e-11, None],
            doc="Mole concentration of solutes in the retentate, x-dependent",
        )
        self.permeate_flow_volume = Var(
            self.time,
            self.x_bar,
            initialize=30,
            units=units.m**3 / units.h,
            bounds=[1e-11, None],
            doc="Volumetric flow rate of the permeate, x-dependent",
        )

        def initialize_permeate_conc_mol_comp(m, t, p, j):
            vals = {
                "Li": 48,
                "Co": 222,
                "Cl": 492,
            }
            return vals[j]

        self.permeate_conc_mol_comp = Var(
            self.time,
            self.x_bar,
            self.solutes,
            initialize=initialize_permeate_conc_mol_comp,
            units=units.mol / units.m**3,  # mM
            bounds=[1e-11, None],
            doc="Mole concentration of solutes in the permeate, x-dependent",
        )
        self.osmotic_pressure = Var(
            self.x_bar,
            initialize=1,
            units=units.bar,
            bounds=[1e-11, None],
            doc="Osmostic pressure of the feed-side fluid",
        )

        # add variables dependent on x_bar and z_bar
        self.membrane_conc_mol_lithium = Var(
            self.x_bar,
            self.z_bar,
            initialize=60,
            units=units.mol / units.m**3,  # mM
            bounds=[1e-11, None],
            doc="Mole concentration of lithium in the membrane, x- and z-dependent",
        )
        self.membrane_conc_mol_cobalt = Var(
            self.x_bar,
            self.z_bar,
            initialize=60,
            units=units.mol / units.m**3,  # mM
            bounds=[1e-11, None],
            doc="Mole concentration of cobalt in the membrane, x- and z-dependent",
        )
        self.membrane_conc_mol_chloride = Var(
            self.x_bar,
            self.z_bar,
            initialize=180,
            units=units.mol / units.m**3,  # mM
            bounds=[1e-11, None],
            doc="Mole concentration of chloride in the membrane, x- and z-dependent",
        )
        self.D_lithium_lithium = Var(
            self.x_bar,
            self.z_bar,
            initialize=1e-6,
            units=units.m**2 / units.h,
            bounds=[1e-11, None],
            doc="Linearized cross diffusion coefficient for lithium-lithium",
        )
        self.D_lithium_cobalt = Var(
            self.x_bar,
            self.z_bar,
            initialize=1e-7,
            units=units.m**2 / units.h,
            bounds=[1e-11, None],
            doc="Linearized cross diffusion coefficient for lithium-cobalt",
        )
        self.D_cobalt_lithium = Var(
            self.x_bar,
            self.z_bar,
            initialize=1e-7,
            units=units.m**2 / units.h,
            bounds=[1e-11, None],
            doc="Linearized cross diffusion coefficient for cobalt-lithium",
        )
        self.D_cobalt_cobalt = Var(
            self.x_bar,
            self.z_bar,
            initialize=1e-6,
            units=units.m**2 / units.h,
            bounds=[1e-11, None],
            doc="Linearized cross diffusion coefficient for cobalt-cobalt",
        )
        self.convection_coefficient_lithium = Var(
            self.x_bar,
            self.z_bar,
            initialize=1,
            units=units.dimensionless,
            doc="Linearized convection coefficient for lithium",
        )
        self.convection_coefficient_cobalt = Var(
            self.x_bar,
            self.z_bar,
            initialize=1,
            units=units.dimensionless,
            doc="Linearized convection coefficient for cobalt",
        )

        # define the (partial) derivative variables
        self.d_retentate_conc_mol_comp_dx = DerivativeVar(
            self.retentate_conc_mol_comp,
            wrt=self.x_bar,
            units=units.mol / units.m**3,  # mM
            doc="Solute concentration gradient in the retentate",
        )
        self.d_retentate_flow_volume_dx = DerivativeVar(
            self.retentate_flow_volume,
            wrt=self.x_bar,
            units=units.m**3 / units.h,
            doc="Volume flow gradient in the retentate",
        )
        self.d_membrane_conc_mol_lithium_dz = DerivativeVar(
            self.membrane_conc_mol_lithium,
            wrt=self.z_bar,
            units=units.mol / units.m**3,  # mM
            doc="Lithium concentration gradient wrt z in the membrane",
        )
        self.d_membrane_conc_mol_cobalt_dz = DerivativeVar(
            self.membrane_conc_mol_cobalt,
            wrt=self.z_bar,
            units=units.mol / units.m**3,  # mM
            doc="Cobalt concentration gradient wrt z in the membrane",
        )

    def add_constraints(self):
        """
        Adds model constraints for the two salt diafiltration unit model.
        """

        # mol balance constraints
        def _overall_mol_balance(self, x):
            if x == 0:
                return Constraint.Skip
            return self.d_retentate_flow_volume_dx[0, x] == (
                -self.volume_flux_water[x] * self.membrane_length * self.membrane_width
            )

        self.overall_mol_balance = Constraint(self.x_bar, rule=_overall_mol_balance)

        def _lithium_mol_balance(self, x):
            if x == 0:
                return Constraint.Skip
            return (
                self.retentate_flow_volume[0, x]
                * self.d_retentate_conc_mol_comp_dx[0, x, "Li"]
            ) == (
                (
                    self.volume_flux_water[x] * self.retentate_conc_mol_comp[0, x, "Li"]
                    - self.mol_flux_lithium[x]
                )
                * self.membrane_length
                * self.membrane_width
            )

        self.lithium_mol_balance = Constraint(self.x_bar, rule=_lithium_mol_balance)

        def _cobalt_mol_balance(self, x):
            if x == 0:
                return Constraint.Skip
            return (
                self.retentate_flow_volume[0, x]
                * self.d_retentate_conc_mol_comp_dx[0, x, "Co"]
            ) == (
                (
                    self.volume_flux_water[x] * self.retentate_conc_mol_comp[0, x, "Co"]
                    - self.mol_flux_cobalt[x]
                )
                * self.membrane_length
                * self.membrane_width
            )

        self.cobalt_mol_balance = Constraint(self.x_bar, rule=_cobalt_mol_balance)

        # transport constraints (geometric)
        def _geometric_flux_equation_overall(self, x):
            if x == 0:
                return Constraint.Skip
            return (
                self.permeate_flow_volume[0, x]
                == self.volume_flux_water[x]
                * x
                * self.membrane_length
                * self.membrane_width
            )

        self.geometric_flux_equation_overall = Constraint(
            self.x_bar, rule=_geometric_flux_equation_overall
        )

        def _geometric_flux_equation_lithium(self, x):
            if x == 0:
                return Constraint.Skip
            return self.mol_flux_lithium[x] == (
                self.permeate_conc_mol_comp[0, x, "Li"] * self.volume_flux_water[x]
            )

        self.geometric_flux_equation_lithium = Constraint(
            self.x_bar, rule=_geometric_flux_equation_lithium
        )

        def _geometric_flux_equation_cobalt(self, x):
            if x == 0:
                return Constraint.Skip
            return self.mol_flux_cobalt[x] == (
                self.permeate_conc_mol_comp[0, x, "Co"] * self.volume_flux_water[x]
            )

        self.geometric_flux_equation_cobalt = Constraint(
            self.x_bar, rule=_geometric_flux_equation_cobalt
        )

        # transport constraints (first principles)
        def _lumped_water_flux(self, x):
            if x == 0:
                return Constraint.Skip
            return self.volume_flux_water[x] == (
                self.membrane_permeability
                * (self.applied_pressure - self.osmotic_pressure[x])
            )

        self.lumped_water_flux = Constraint(self.x_bar, rule=_lumped_water_flux)

        # TODO: write exceptions for the valid membrane concentration ranges
        def _D_lithium_lithium_calculation(self, x, z):
            if value(self.membrane_fixed_charge) == 0:
                # these params assume chi=0 mM and lithium & cobalt mem_conc range 50-80 mM
                self.diffusion_params_lithium_lithium = [-4.07e-06, -3.99e-09, 4.01e-09]
            if value(self.membrane_fixed_charge) == -140:
                # these params assume chi=-140 mM and lithium mem_conc range 50-80 mM & cobalt mem_conc range 80-110 mM
                self.diffusion_params_lithium_lithium = [-4.33e-06, -4.25e-09, 5.14e-09]
            return self.D_lithium_lithium[x, z] == -(
                (self.diffusion_params_lithium_lithium[0] * units.m**2 / units.h)
                + (
                    (
                        self.diffusion_params_lithium_lithium[1]
                        * units.m**5
                        / units.mol
                        / units.h
                    )
                    * (self.membrane_conc_mol_lithium[x, z])
                )
                + (
                    (
                        self.diffusion_params_lithium_lithium[2]
                        * units.m**5
                        / units.mol
                        / units.h
                    )
                    * (self.membrane_conc_mol_cobalt[x, z])
                )
            )

        self.D_lithium_lithium_calculation = Constraint(
            self.x_bar, self.z_bar, rule=_D_lithium_lithium_calculation
        )

        def _D_lithium_cobalt_calculation(self, x, z):
            if value(self.membrane_fixed_charge) == 0:
                # these params assume chi=0 mM and lithium & cobalt mem_conc range 50-80 mM
                self.diffusion_params_lithium_cobalt = [-9.62e-07, -1.03e-08, 1.04e-08]
            if value(self.membrane_fixed_charge) == -140:
                # these params assume chi=-140 mM and lithium mem_conc range 50-80 mM & cobalt mem_conc range 80-110 mM
                self.diffusion_params_lithium_cobalt = [-1.63e-06, -1.10e-08, 1.33e-08]
            return self.D_lithium_cobalt[x, z] == -(
                (self.diffusion_params_lithium_cobalt[0] * units.m**2 / units.h)
                + (
                    (
                        self.diffusion_params_lithium_cobalt[1]
                        * units.m**5
                        / units.mol
                        / units.h
                    )
                    * (self.membrane_conc_mol_lithium[x, z])
                )
                + (
                    (
                        self.diffusion_params_lithium_cobalt[2]
                        * units.m**5
                        / units.mol
                        / units.h
                    )
                    * (self.membrane_conc_mol_cobalt[x, z])
                )
            )

        self.D_lithium_cobalt_calculation = Constraint(
            self.x_bar, self.z_bar, rule=_D_lithium_cobalt_calculation
        )

        def _D_cobalt_lithium_calculation(self, x, z):
            if value(self.membrane_fixed_charge) == 0:
                # these params assume chi=0 mM and lithium & cobalt mem_conc range 50-80 mM
                self.diffusion_param_cobalt_lithium = [-5.24e-07, 2.49e-09, -2.50e-09]
            if value(self.membrane_fixed_charge) == -140:
                # these params assume chi=-140 mM and lithium mem_conc range 50-80 mM & cobalt mem_conc range 80-110 mM
                self.diffusion_param_cobalt_lithium = [-1.32e-06, 4.72e-09, 1.47e-09]
            return self.D_cobalt_lithium[x, z] == -(
                (self.diffusion_param_cobalt_lithium[0] * units.m**2 / units.h)
                + (
                    (
                        self.diffusion_param_cobalt_lithium[1]
                        * units.m**5
                        / units.mol
                        / units.h
                    )
                    * (self.membrane_conc_mol_lithium[x, z])
                )
                + (
                    (
                        self.diffusion_param_cobalt_lithium[2]
                        * units.m**5
                        / units.mol
                        / units.h
                    )
                    * (self.membrane_conc_mol_cobalt[x, z])
                )
            )

        self.D_cobalt_lithium_calculation = Constraint(
            self.x_bar, self.z_bar, rule=_D_cobalt_lithium_calculation
        )

        def _D_cobalt_cobalt_calculation(self, x, z):
            if value(self.membrane_fixed_charge) == 0:
                # these params assume chi=0 mM and lithium & cobalt mem_conc range 50-80 mM
                self.diffusion_params_cobalt_cobalt = [-4.00e-06, 6.45e-09, -6.48e-09]
            if value(self.membrane_fixed_charge) == -140:
                # these params assume chi=-140 mM and lithium mem_conc range 50-80 mM & cobalt mem_conc range 80-110 mM
                self.diffusion_params_cobalt_cobalt = [-6.05e-06, 1.22e-08, 3.81e-09]
            return self.D_cobalt_cobalt[x, z] == -(
                (self.diffusion_params_cobalt_cobalt[0] * units.m**2 / units.h)
                + (
                    (
                        self.diffusion_params_cobalt_cobalt[1]
                        * units.m**5
                        / units.mol
                        / units.h
                    )
                    * (self.membrane_conc_mol_lithium[x, z])
                )
                + (
                    (
                        self.diffusion_params_cobalt_cobalt[2]
                        * units.m**5
                        / units.mol
                        / units.h
                    )
                    * (self.membrane_conc_mol_cobalt[x, z])
                )
            )

        self.D_cobalt_cobalt_calculation = Constraint(
            self.x_bar, self.z_bar, rule=_D_cobalt_cobalt_calculation
        )

        def _convection_coefficient_lithium_calculation(self, x, z):
            if value(self.membrane_fixed_charge) == 0:
                # these params assume chi=0 mM and lithium & cobalt mem_conc range 50-80 mM
                self.convection_params_lithium = [1, 0, 0]
            if value(self.membrane_fixed_charge) == -140:
                # these params assume chi=-140 mM and lithium mem_conc range 50-80 mM & cobalt mem_conc range 80-110 mM
                self.convection_params_lithium = [0.260, 0.00139, 0.00314]
            return self.convection_coefficient_lithium[x, z] == (
                self.convection_params_lithium[0]
                + (
                    self.convection_params_lithium[1]
                    * units.m**3
                    / units.mol
                    * self.membrane_conc_mol_lithium[x, z]
                )
                + (
                    self.convection_params_lithium[2]
                    * units.m**3
                    / units.mol
                    * self.membrane_conc_mol_cobalt[x, z]
                )
            )

        self.convection_coefficient_lithium_calculation = Constraint(
            self.x_bar, self.z_bar, rule=_convection_coefficient_lithium_calculation
        )

        def _convection_coefficient_cobalt_calculation(self, x, z):
            if value(self.membrane_fixed_charge) == 0:
                # these params assume chi=0 mM and lithium & cobalt mem_conc range 50-80 mM
                self.convection_params_cobalt = [1, 0, 0]
            if value(self.membrane_fixed_charge) == -140:
                # these params assume chi=-140 mM and lithium mem_conc range 50-80 mM & cobalt mem_conc range 80-110 mM
                self.convection_params_cobalt = [0.0860, 0.00198, 0.00448]
            return self.convection_coefficient_cobalt[x, z] == (
                self.convection_params_cobalt[0]
                + (
                    self.convection_params_cobalt[1]
                    * units.m**3
                    / units.mol
                    * self.membrane_conc_mol_lithium[x, z]
                )
                + (
                    self.convection_params_cobalt[2]
                    * units.m**3
                    / units.mol
                    * self.membrane_conc_mol_cobalt[x, z]
                )
            )

        self.convection_coefficient_cobalt_calculation = Constraint(
            self.x_bar, self.z_bar, rule=_convection_coefficient_cobalt_calculation
        )

        def _lithium_flux_membrane(self, x, z):
            if x == 0:
                return Constraint.Skip
            return self.mol_flux_lithium[x] == (
                (
                    self.convection_coefficient_lithium[x, z]
                    * self.membrane_conc_mol_lithium[x, z]
                    * self.volume_flux_water[x]
                )
                - (
                    self.D_lithium_lithium[x, z]
                    / self.membrane_thickness
                    * self.d_membrane_conc_mol_lithium_dz[x, z]
                )
                - (
                    self.D_lithium_cobalt[x, z]
                    / self.membrane_thickness
                    * self.d_membrane_conc_mol_cobalt_dz[x, z]
                )
            )

        self.lithium_flux_membrane = Constraint(
            self.x_bar, self.z_bar, rule=_lithium_flux_membrane
        )

        def _cobalt_flux_membrane(self, x, z):
            if x == 0:
                return Constraint.Skip
            return self.mol_flux_cobalt[x] == (
                (
                    self.convection_coefficient_cobalt[x, z]
                    * self.membrane_conc_mol_cobalt[x, z]
                    * self.volume_flux_water[x]
                )
                - (
                    self.D_cobalt_lithium[x, z]
                    / self.membrane_thickness
                    * self.d_membrane_conc_mol_lithium_dz[x, z]
                )
                - (
                    self.D_cobalt_cobalt[x, z]
                    / self.membrane_thickness
                    * self.d_membrane_conc_mol_cobalt_dz[x, z]
                )
            )

        self.cobalt_flux_membrane = Constraint(
            self.x_bar, self.z_bar, rule=_cobalt_flux_membrane
        )

        def _chloride_flux_membrane(self, x):
            if x == 0:
                return Constraint.Skip
            return 0 == (
                (self.config.property_package.charge["Li"] * self.mol_flux_lithium[x])
                + (self.config.property_package.charge["Co"] * self.mol_flux_cobalt[x])
                + (
                    self.config.property_package.charge["Cl"]
                    * self.mol_flux_chloride[x]
                )
            )

        self.chloride_flux_membrane = Constraint(
            self.x_bar, rule=_chloride_flux_membrane
        )

        # other physical constraints
        def _osmotic_pressure_calculation(self, x):
            return self.osmotic_pressure[x] == units.convert(
                (
                    Constants.gas_constant  # J / mol / K
                    * self.temperature
                    * (
                        (
                            self.config.property_package.num_solutes["Li"]
                            * self.config.property_package.sigma["Li"]
                            * (
                                self.retentate_conc_mol_comp[0, x, "Li"]
                                - self.permeate_conc_mol_comp[0, x, "Li"]
                            )
                        )
                        + (
                            self.config.property_package.num_solutes["Co"]
                            * self.config.property_package.sigma["Co"]
                            * (
                                self.retentate_conc_mol_comp[0, x, "Co"]
                                - self.permeate_conc_mol_comp[0, x, "Co"]
                            )
                        )
                        + (
                            self.config.property_package.num_solutes["Cl"]
                            * self.config.property_package.sigma["Cl"]
                            * (
                                self.retentate_conc_mol_comp[0, x, "Cl"]
                                - self.permeate_conc_mol_comp[0, x, "Cl"]
                            )
                        )
                    )
                ),
                to_units=units.bar,
            )

        self.osmotic_pressure_calculation = Constraint(
            self.x_bar, rule=_osmotic_pressure_calculation
        )

        def _electroneutrality_retentate(self, x):
            return 0 == (
                self.config.property_package.charge["Li"]
                * self.retentate_conc_mol_comp[0, x, "Li"]
                + self.config.property_package.charge["Co"]
                * self.retentate_conc_mol_comp[0, x, "Co"]
                + self.config.property_package.charge["Cl"]
                * self.retentate_conc_mol_comp[0, x, "Cl"]
            )

        self.electroneutrality_retentate = Constraint(
            self.x_bar, rule=_electroneutrality_retentate
        )

        def _electroneutrality_membrane(self, x, z):
            if x == 0:
                return Constraint.Skip
            return 0 == (
                self.config.property_package.charge["Li"]
                * self.membrane_conc_mol_lithium[x, z]
                + self.config.property_package.charge["Co"]
                * self.membrane_conc_mol_cobalt[x, z]
                + self.config.property_package.charge["Cl"]
                * self.membrane_conc_mol_chloride[x, z]
                + self.membrane_fixed_charge
            )

        self.electroneutrality_membrane = Constraint(
            self.x_bar, self.z_bar, rule=_electroneutrality_membrane
        )

        def _electroneutrality_permeate(self, x):
            return 0 == (
                self.config.property_package.charge["Li"]
                * self.permeate_conc_mol_comp[0, x, "Li"]
                + self.config.property_package.charge["Co"]
                * self.permeate_conc_mol_comp[0, x, "Co"]
                + self.config.property_package.charge["Cl"]
                * self.permeate_conc_mol_comp[0, x, "Cl"]
            )

        self.electroneutrality_permeate = Constraint(
            self.x_bar, rule=_electroneutrality_permeate
        )

        # boundary conditions
        # TODO: investigate numerics of substituting out H_Cl
        def _retentate_membrane_interface_lithium(self, x):
            if x == 0:
                return Constraint.Skip
            return (
                self.config.property_package.partition_coefficient["Li"]
                * self.config.property_package.partition_coefficient["Cl"]
                * self.retentate_conc_mol_comp[0, x, "Li"]
                * self.retentate_conc_mol_comp[0, x, "Cl"]
            ) == (
                self.membrane_conc_mol_lithium[x, 0]
                * self.membrane_conc_mol_chloride[x, 0]
            )

        self.retentate_membrane_interface_lithium = Constraint(
            self.x_bar, rule=_retentate_membrane_interface_lithium
        )

        def _retentate_membrane_interface_cobalt(self, x):
            if x == 0:
                return Constraint.Skip
            return (
                self.config.property_package.partition_coefficient["Co"]
                * self.config.property_package.partition_coefficient["Cl"]
                ** self.config.property_package.charge["Co"]
                * self.retentate_conc_mol_comp[0, x, "Co"]
                * self.retentate_conc_mol_comp[0, x, "Cl"]
                ** self.config.property_package.charge["Co"]
            ) == (
                self.membrane_conc_mol_cobalt[x, 0]
                * self.membrane_conc_mol_chloride[x, 0]
                ** self.config.property_package.charge["Co"]
            )

        self.retentate_membrane_interface_cobalt = Constraint(
            self.x_bar, rule=_retentate_membrane_interface_cobalt
        )

        def _membrane_permeate_interface_lithium(self, x):
            if x == 0:
                return Constraint.Skip
            return (
                self.config.property_package.partition_coefficient["Li"]
                * self.config.property_package.partition_coefficient["Cl"]
                * self.permeate_conc_mol_comp[0, x, "Li"]
                * self.permeate_conc_mol_comp[0, x, "Cl"]
            ) == (
                self.membrane_conc_mol_lithium[x, 1]
                * self.membrane_conc_mol_chloride[x, 1]
            )

        self.membrane_permeate_interface_lithium = Constraint(
            self.x_bar, rule=_membrane_permeate_interface_lithium
        )

        def _membrane_permeate_interface_cobalt(self, x):
            if x == 0:
                return Constraint.Skip
            return (
                self.config.property_package.partition_coefficient["Co"]
                * self.config.property_package.partition_coefficient["Cl"]
                ** self.config.property_package.charge["Co"]
                * self.permeate_conc_mol_comp[0, x, "Co"]
                * self.permeate_conc_mol_comp[0, x, "Cl"]
                ** self.config.property_package.charge["Co"]
            ) == (
                self.membrane_conc_mol_cobalt[x, 1]
                * self.membrane_conc_mol_chloride[x, 1]
                ** self.config.property_package.charge["Co"]
            )

        self.membrane_permeate_interface_cobalt = Constraint(
            self.x_bar, rule=_membrane_permeate_interface_cobalt
        )

    def discretize_model(self):
        discretizer = TransformationFactory("dae.finite_difference")
        discretizer.apply_to(
            self, wrt=self.x_bar, nfe=self.config.NFEx, scheme="BACKWARD"
        )
        discretizer.apply_to(
            self, wrt=self.z_bar, nfe=self.config.NFEz, scheme="BACKWARD"
        )

    def fix_initial_values(self):
        """
        Fix initial values for the two salt diafiltration unit model.
        """
        for x in self.x_bar:
            # chloride concentration gradient in retentate variable is created by default but
            # is not needed in model; fix to reduce number of variables
            self.d_retentate_conc_mol_comp_dx[0, x, "Cl"].fix(
                value(self.numerical_zero_tolerance)
            )
            # associated discretization equation not needed in model
            if x != 0:
                self.d_retentate_conc_mol_comp_dx_disc_eq[0, x, "Cl"].deactivate()

        # initial conditions
        self.retentate_flow_volume[0, 0].fix(
            self.feed_flow_volume[0] + self.diafiltrate_flow_volume[0]
        )
        self.retentate_conc_mol_comp[0, 0, "Li"].fix(
            (
                self.feed_flow_volume[0] * self.feed_conc_mol_comp[0, "Li"]
                + self.diafiltrate_flow_volume[0]
                * self.diafiltrate_conc_mol_comp[0, "Li"]
            )
            / (self.feed_flow_volume[0] + self.diafiltrate_flow_volume[0])
        )
        self.retentate_conc_mol_comp[0, 0, "Co"].fix(
            (
                self.feed_flow_volume[0] * self.feed_conc_mol_comp[0, "Co"]
                + self.diafiltrate_flow_volume[0]
                * self.diafiltrate_conc_mol_comp[0, "Co"]
            )
            / (self.feed_flow_volume[0] + self.diafiltrate_flow_volume[0])
        )

        # set "zero" initial values to a sufficiently small value
        self.permeate_flow_volume[0, 0].fix(value(self.numerical_zero_tolerance))
        self.permeate_conc_mol_comp[0, 0, "Li"].fix(
            value(self.numerical_zero_tolerance)
        )
        self.permeate_conc_mol_comp[0, 0, "Co"].fix(
            value(self.numerical_zero_tolerance)
        )
        for z in self.z_bar:
            self.membrane_conc_mol_lithium[0, z].fix(
                value(self.numerical_zero_tolerance)
            )
            self.membrane_conc_mol_cobalt[0, z].fix(
                value(self.numerical_zero_tolerance)
            )
            self.membrane_conc_mol_chloride[0, z].fix(
                value(self.numerical_zero_tolerance)
            )
        self.d_retentate_conc_mol_comp_dx[0, 0, "Li"].fix(
            value(self.numerical_zero_tolerance)
        )
        self.d_retentate_conc_mol_comp_dx[0, 0, "Co"].fix(
            value(self.numerical_zero_tolerance)
        )
        self.d_retentate_flow_volume_dx[0, 0].fix(value(self.numerical_zero_tolerance))

        # set fluxes at x=0 to numerically 0 (expected to be 0)
        self.volume_flux_water[0].fix(value(self.numerical_zero_tolerance))
        self.mol_flux_lithium[0].fix(value(self.numerical_zero_tolerance))
        self.mol_flux_cobalt[0].fix(value(self.numerical_zero_tolerance))
        self.mol_flux_chloride[0].fix(value(self.numerical_zero_tolerance))

    def add_scaling_factors(self):
        """
        Assigns scaling factors to certain variables and constraints to
        improve solver performance.
        """
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Add scaling factors for poorly scaled variables
        for x in self.x_bar:
            if x != 0:
                self.scaling_factor[self.retentate_flow_volume[0, x]] = 1e-2
                self.scaling_factor[self.permeate_flow_volume[0, x]] = 1e-2

                self.scaling_factor[self.volume_flux_water[x]] = 1e-1
                self.scaling_factor[self.mol_flux_lithium[x]] = 1e-1
                self.scaling_factor[self.mol_flux_cobalt[x]] = 1e-1
                self.scaling_factor[self.mol_flux_chloride[x]] = 1e-2

            for z in self.z_bar:
                self.scaling_factor[self.D_lithium_lithium[x, z]] = 1e6
                self.scaling_factor[self.D_lithium_cobalt[x, z]] = 1e7
                self.scaling_factor[self.D_cobalt_lithium[x, z]] = 1e7
                self.scaling_factor[self.D_cobalt_cobalt[x, z]] = 1e6

                self.scaling_factor[self.membrane_conc_mol_lithium[x, z]] = 1e-1
                self.scaling_factor[self.membrane_conc_mol_cobalt[x, z]] = 1e-1
                self.scaling_factor[self.membrane_conc_mol_chloride[x, z]] = 1e-2

        # Add scaling factors for poorly scaled constraints
        constraint_autoscale_large_jac(self)

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
            self.retentate_flow_volume[:, self.x_bar.last()]
        )
        self.retentate_outlet.add(self._retentate_flow_volume_ref, "flow_vol")
        self._retentate_conc_mol_comp_ref = Reference(
            self.retentate_conc_mol_comp[:, self.x_bar.last(), :]
        )
        self.retentate_outlet.add(self._retentate_conc_mol_comp_ref, "conc_mol_comp")

        self.permeate_outlet = Port(doc="Permeate Outlet Port")
        self._permeate_flow_volume_ref = Reference(
            self.permeate_flow_volume[:, self.x_bar.last()]
        )
        self.permeate_outlet.add(self._permeate_flow_volume_ref, "flow_vol")
        self._permeate_conc_mol_comp_ref = Reference(
            self.permeate_conc_mol_comp[:, self.x_bar.last(), :]
        )
        self.permeate_outlet.add(self._permeate_conc_mol_comp_ref, "conc_mol_comp")
