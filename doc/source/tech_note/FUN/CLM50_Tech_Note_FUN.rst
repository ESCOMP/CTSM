.. _rst_FUN:

Fixation and Uptake of Nitrogen (FUN)
=======================================

Introduction
-----------------

The Fixation and Uptake of Nitrogen model is based on work by :ref:`Fisher et al. (2010)<Fisheretal2010>`, :ref:`Brzostek et al. (2014)<Brzosteketal2014>`, and :ref:`Shi et al. (2016)<Shietal2016>`.  The concept of FUN is that in most cases, Nitrogen uptake requires the expenditure of energy in the form of carbon, and further, that there are numerous potential sources of Nitrogen in the environment which a plant may exchange for carbon. The ratio of carbon expended to Nitrogen acquired is referred to here as the cost, or exchange rate,  of N acquisition (:math:`E_{nacq}`, gC/gN)). There are eight pathways for N uptake:

1. Fixation by symbiotic bacteria in root nodules (for N fixing plants) (:math:`_{fix}`)
2. Retranslocation of N from senescing tissues (:math:`_{ret}`)
3. Active uptake of NH4 by arbuscular mycorrhizal plants (:math:`_{active,nh4}`)
4. Active uptake of NH4 by ectomycorrhizal plants (:math:`_{active,nh4}`)
5. Active uptake of NO3 by arbuscular mycorrhizal plants (:math:`_{active,no3}`)
6. Active uptake of NO3 by ectomycorrhizal plants (:math:`_{active,no3}`)
7. Nonmycorrhizal uptake of NH4 (:math:`_{nonmyc,no3}`)
8. Nonmycorrhizal uptake of NO3 (:math:`_{nonmyc,nh4}`)

The notation suffix for each pathway is given in parentheses here. At each timestep, each of these pathways is associated with a cost term (:math:`N_{cost,x}`), a payment in carbon (:math:`C_{nuptake,x}`), and an influx of Nitrogen (:math:`N_{uptake,x}`) where :math:`x` is one of the eight uptake streams listed above.

For each PFT, we define a fraction of the total C acquisition that can be used for N fixation (:math:`f_{fixers}`), which is broadly equivalent to the fraction of a given PFT that is capable of fixing Nitrogen, and thus represents an upper limit on the amount to which fixation can be increased in low n conditions.  For each PFT, the cost calculation is conducted twice. Once where fixation is possible and once where it is not. (:math:`f_{fixers}`)

For all of the active uptake pathways, whose cost depends on varying concentrations of N through the soil profile, the costs and fluxes are also determined by soil layer :math:`j`.

Boundary conditions of FUN
--------------------------------------------------------

Available Carbon
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The carbon available for FUN, :math:`C_{avail}` (gC m\ :sup:`-2`) is the total canopy  photosynthetic uptake (GPP), minus the maintenance respiration fluxes (:math:`m_r`) and multiplied by the time step in seconds (:math:`\delta t`). Thus, the remainder of this chapter considers fluxes per timestep, and integrates these fluxes as they are calculated.

 .. math::

   C_{avail} = (GPP - m_r) \delta t

Growth respiration is thus only calculated on the part of the carbon uptake that remains after expenditure of C by the FUN module.

Available Soil Nitrogen
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Cost of Nitrogen Fixation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The cost of fixation is derived from :ref:`Houlton et al. (2008)<Houltonetal2008>`.
 .. math::

   N_{cost,fix} = -s_{fix}/(1.25 e^{a_{fix} + b_{fix} . t_{soil}  (1 - 0.5 t_{soil}/ c_{fix}) })

Herein, :math:`a_{fix}`, :math:`b_{fix}` and :math:`c_{fix}` are all parameters of the temperature response function of fixation reported by Houlton et al. (2008) (:math:`exp[a+bT_s(1-0.5T_s/c)`).   t_{soil} is the soil temperature in C. The values of these parameters are fitted to empirical data as a=-3.62 :math:`\pm` 0.52, b=0.27:math:`\pm` 0.04 and c=25.15 :math:`\pm` 0.66. 1.25 converts from the temperature response function to a 0-1 limitation factor (as specifically employed by Houlton et al.).  This function is a 'rate' of uptake for a given temperature. Here we assimilated the rate of fixation into the cost term by assuming that the rate is analagous to a conductance for N, and inverting the term to produce a cost/resistance analagoue. We then multiply this temperature term by the minimum cost at optimal temperature (:math:`s_{fix}`) to give a temperature limited cost in terms of C to N ratios.

Cost of Active Uptake
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The cost of N uptake from soil, for each layer :math:`j`, is controlled by two uptake parameters that pertain respectively to the relationship between soil N content and N uptake, and root C density and N uptake.

For non-mycorrhizal uptake:

 .. math::

   N_{cost,nonmyc,j} = \frac{k_{n,nonmyc}}{N_{smin,j}} + \frac{k_{c,nonmyc}}{c_{root,j}}

and for active uptake:

 .. math::

   N_{cost,active,j} = \frac{k_{n,active}}{N_{smin,j}} + \frac{k_{c,active}}{c_{root,j}}

where :math:`k_{n,active}` varies according to whether we are considering ecto or arbuscular mycorrhizal uptake.

 .. math::
   :label: 18.2

   k_{n,active}  =
   \left\{\begin{array}{lr}
   k_{n,Eactive}&  e = 1\\
   k_{n,Aactive}&  e = 0
   \end{array}\right\}

where m=1 pertains to the fraction of the PFT that is ecotmycorrhizal, as opposed to arbuscular mycorrhizal.

Resolving N cost across simultaneous uptake streams
--------------------------------------------------------
The total cost of N uptake is calculated based on the assumption that carbon is partitioned to each stream in proportion to the inverse of the cost of uptake. So, more expensive pathways receive less carbon. Earlier versions of FUN :ref:`(Fisher et al., 2010)<Fisheretal2010>`) utilized a scheme whereby plants only took up N from the cheapest pathway. :ref:`Brzostek et al. (2014)<Brzosteketal2014>` introduced a scheme for the simultaneous uptake from different pathways. Here we calcualate a 'conductance' to N uptake (analagous to the inverse of the cost function conceptualized as a resistance term) :math:`N_{conductance}` ( gN/gC) as:

 .. math::

   N_{conductance,f}=  \sum{(1/N_{cost,x})}

From this, we then calculate the fraction of the carbon allocated to each pathway as

 .. math::

   C_{frac,x} = \frac{1/N_{cost,x}}{N_{conductance}}

These fractions are used later, to calculate the carbon expended on different uptake pathways.  Next, the N acquired from each uptake stream per unit C spent (:math:`N_{exch,x}`, gN/gC)  is determined as

 .. math::

   N_{exch,x} = \frac{C_{frac,x}}{N_{cost,x}}

We then determine the total amount of N uptake per unit C spent (:math:`N_{exch,tot}`, gN/gC) as the sum of all the uptake streams.

 .. math::
   N_{exch,tot} = \sum{N_{exch,x}}

and thus the subsequent overall N cost is

 .. math::
   N_{cost,tot} = 1/{N_{exch,tot}}

 Retranslocation is determined via a different set of mechanisms, once the :math:`N_{cost,tot}` is known.

Nitrogen Retranslocation
--------------------------------------------------------
The retranslocation uses an iterative algorithm to remove Nitrogen from each piece of falling litter.  There are two pathways for this, 'free' uptake which removes the labile N pool, and 'paid-for' uptake which uses C to extract N from increasingly more recalcitrant pools.

At each timestep, the pool of carbon in falling leaves (:math:`C_{fallingleaf}`, g m\ :sup:`-2`) is generated from the quantity of litterfall on that day (see Phenology chapter for details). The amount of N in the litter pool (:math:`N_{fallingleaf}`, g m\ :sup:`-2`) is calculated as the total leaf N multiplied by the fraction of the leaf pool passed to litter that timestep.

 .. math::

  N_{fallingleaf} = N_{leaf}.C_{fallingleaf}/C_{leaf}

The carbon available at the beginning of the iterative retranslocation calculation is equal to the :math:`C_{avail}` input into FUN.

 .. math::

  C_{avail,retrans,0} = C_{avail}

Free Retranslocation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Some part of the leaf Nitrogen pool is removed without the need for an C expenditure.  This 'free' N uptake amount, (:math:`N_{retrans,free}`, gN m\ :sup:`-2`) is calculated as

 .. math::

  N_{retrans,free}  = max(N_{fallingleaf} -  (C_{fallingleaf}/CN_{litter,min} ),0.0)

where :math:`CN_{litter,min}` is the minimum C:N ratio of the falling litter (currently set to 1.5 x the target C:N ratio).

The new :math:`N_{fallingleaf}` (gN m\ :sup:`-2`) is then determined as

 .. math::

  N_{fallingleaf} = N_{fallingleaf} - N_{retrans,free}

and the new litter C:N ratio as

 .. math::

  CN_{fallingleaf}=C_{fallingleaf}/N_{fallingleaf}

Paid-for Retranslocation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The remaining calculations conduct an iterative calculation to determine the degree to which N retranslocation from leaves is paid for as C:N ratios and thus cost increase as N is extracted.  The iteration continues until either

1. The cost of retranslocation (:math:`cost_{retrans}` increases beyond the cost of acquiring N from alternative pathways (:math:`N_{cost,tot}`).
2. :math:`CN_{fallingleaf}` rises to a maximum level, after which no more extraction is possible (representing unavoidable N loss) or
3. There is no more carbon left to pay for extraction.

First we calculate the cost of extraction (:math:`cost_{retrans}`, gC/gN) for the current leaf C:N ratio as

 .. math::

  cost_{retrans}= k_{retrans} / (1/CN_{fallingleaf})^{1.3}

where :math:`k_{retrans}`  is a parameter controlling the overall cost of resorption, which also increases exponentially as the C:N ratio increases

Next, we calculate the amount of C needed to be spent to increase the falling leaf C:N ratio by 1.0 in this iteration :math:`i` (:math:`C_{retrans_spent,i}`,  gC m\ :sup:`-2`) as:
 .. math::

  C_{retrans,spent,i}   = cost_{retrans}.(N_{fallingleaf} - C_{fallingleaf}/
                          (CN_{fallingleaf} + 1.0))

(wherein the retranslocation cost is assumed to not change over the increment of 1.0 in C:N ratio).   Next, we calculate whether this is larger than the remaining C available to spend.

 .. math::

  C_{retrans,spent,i} = min(C_{retrans,spent,i}, C_{avail,retrans,i})

The amount of N retranslocated from the leaf in this iteration (:math:`N_{retrans_paid,i}`,  gN m\ :sup:`-2`) is calculated, checking that it does not fall below zero:

 .. math::

  N_{retrans,paid,i} = min(N_{fallingleaf},C_{retrans,spent,i} / cost_{retrans})

The next step calculates the growth C which is accounted for by this amount of N extraction in this iteration (:math:`C_{retrans,accounted,i}`).  This is calculated using the current plant C:N ratio, and also for the additional C which will need to be spent on growth respiration to build this amount of new tissue.

 .. math::

  C_{retrans,accounted,i} = N_{retrans,paid,i} . CN_{plant} . (1.0 + gr_{frac})

Then the falling leaf N is updated:

 .. math::

  N_{fallingleaf}    = N_{fallingleaf} - N_{ret,i}

and the :math:`CN_{fallingleaf}` and cost_{retrans} are updated. The amount of available carbon that is either unspent on N acquisition nor accounted for by N uptake is updated:

 .. math::

  C_{avail,retrans,i+1}  = C_{avail,retrans,i} - C_{retrans,spent,i} - C_{retrans,accounted,i}

Outputs of Retranslocation algorithm.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The final output of the retranslocation calculation are the retranslocated N (:math:`N_{retrans}`,  gN m\ :sup:`-2`), C spent on retranslocation (:math:`C_{retrans_paid}`,  gC m\ :sup:`-2`), and C accounted for by retranslocation (:math:`C_{retrans_accounted}`,  gC m\ :sup:`-2`).

For paid-for uptake, we accumulate the total carbon spent on retranslocation (:math:`C_{spent_retrans}`),

 .. math::

  C_{retrans,spent} = \sum{C_{retrans,i}}

The total N acquired from retranslocation is

 .. math::

  N_{retrans} = N_{retrans,paid}+N_{retrans,free}

where N acquired by paid-for retranslocation is

 .. math::

  N_{retrans,paid} = \sum{N_{retrans,paid,i}}

The total carbon accounted for by retranslocation is the sum of the C accounted for by paid-for N uptake (:math:`N_{retrans_paid}`) and by free N uptake (:math:`N_{retrans_free}`).

 .. math::

  C_{retrans,accounted} = \sum{C_{retrans,accounted,i}}+N_{retrans,free}.CN_{plant} . (1.0 + gr_{frac})

The total available carbon in FUN to spend on fixation and active uptake (:math:`C_{tospend}`,  gC m\ :sup:`-2`) is calculated as the carbon available minus that account for by retranslocation:

 .. math::

  C_{tospend} = C_{avail} - C_{retrans,accounted}

Carbon expenditure on fixation and active uptake.
--------------------------------------------------------

At each model timestep, the overall cost of N uptake is calculated (see below) in terms of C:N ratios. The available carbon (:math:`C_{avail}`, g m\ :sup:`-2` s\ :sup:`-1`) is then allocated to two alternative outcomes, payment for N uptake, or conservation for growth. For each carbon conserved for growth, a corresponding quantity of N must be made available.  In the case where the plant target C:N ratio is fixed, the partitioning between carbon for growth (:math:`C_{growth}`) and carbon for N uptake  (:math:`C_{nuptake}`) is calculated by solving a system of simultaneous equations. First, the carbon available must equal the carbon spent on N uptake plus that saved for growth.

 .. math::

   C_{growth}+C_{nuptake}=C_{avail}

Second, the nitrogen acquired from expenditure of N (left hand side of term below) must equal the N that is required to match the growth carbon (right hand side of term below).

 .. math::

   C_{nuptake}/N_{cost} =C_{growth}/CN_{target}

The solution to these two equated terms can be used to estimate the ideal :math:`C_{nuptake}` as follows,

 .. math::
   C_{nuptake} =C_{tospend}/ ( (1.0+f_{gr}*(CN_{target} / N_{cost}) + 1) .

and the other C and N fluxes can be determined following the logic above.

Modifications to allow variation in C:N ratios
--------------------------------------------------------
The original FUN model as developed by :ref:`Fisher et al. (2010)<Fisheretal2010>` and :ref:`Brzostek et al. (2014)<Brzosteketal2014>` assumes a fixed plant tissue C:N ratio. This means that in the case where N is especially limiting, all excess carbon will be utilized in an attempt to take up more Nitrogen. It has been repeatedly observed, however, that in these circumstances in real life, plants have some flexibility in the C:N stoichiometry of their tissues, and therefore, this assumption may not be realistic. However, the degree to which the C:N ratio varies with N availability is poorly documented, and existing global nitrogen models use a variety of heuristic methods by which to incorporate changing C:N ratios (Zaehle and Friend 2010; Ghimire et al. 2016). This algorithm exists as a placeholder to allow variable C:N ratios to occur, and to allow exploration of how much the parameters controlling their flexibility has on model outcomes. Incorporation of emerging understanding of the controls on tissue stoichiometry should ultimately replace this scheme.

Thus, in CLM5, we introduce the capacity for tissue C:N ratios to be prognostic, rather than static. Overall N and C availability (:math:`N_{uptake}` and :math:`C_{growth}`) and hence tissue C:N ratios, are both determined by FUN.  Allocation to individual tissues is discussed in the allocation chapter

Here we introduce an algorithm which adjusts the C expenditure on uptake to allow varying tissue C:N ratios. Increasing C spent on uptake will directly reduce the C:N ratio, and reducing C spent on uptake (retaining more for tissue growth) will increase it. C spent on uptake is impacted by both the N cost in the environment, and the existing tissue C:N ratio of the plant. The output of this algorithm is :math:`\gamma_{FUN}`, the fraction of the ideal :math:`C_{nuptake}` calculated from the FUN equation above

 .. math::
   C_{nuptake} = C_{nuptake}.\gamma_{FUN}

Response of C expenditure to Nitrogen uptake cost
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The environmental cost of Nitrogen (:math:`N_{cost,tot}`) is used to determine :math:`\gamma_{FUN}`.

 .. math::
   \gamma_{FUN} = max(0.0,1.0 - (N_{cost,tot}-a_{cnflex})/b_{cnflex})

where :math:`a_{cnflex}` and :math:`b_{cnflex}` are parameters fitted to give flexible C:N ranges over the operating range of N costs of the model. Calibration of these parameters should be subject to future testing in idealized experimental settings; they are here intended as a placeholder to allow some flexible stoichiometry, in the absence of adequate understanding of this process.  Here :math:`a_{cnflex}` operates as the :math:`N_{cost,tot}` above which there is a modification in the C expenditure (to allow higher C:N ratios), and :math:`b_{cnflex}` is the scalar which determines how much the C expenditure is modified for a given discrepancy between :math:`a_{cnflex}` and the actual cost of uptake.

Response of C expenditure to plant C:N ratios
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We first calculate a :math:`\delta_{CN}`, which is the difference between the target C:N (:math:`target_{CN}`) a model parameter, and the existing C:N ratio (:math:`CN_{plant}`)

 .. math::

  CN_{plant} = \frac{C_{leaf} + C_{leaf,storage}}{N_{leaf} + N_{leaf,storage})}

and
 .. math::
   \delta_{CN} = CN_{plant} - target_{CN}

We then increase :math:`\gamma_{FUN}` to  account for situations where (even if N is expensive) plant C:N ratios have increased too far from the target.  Where  :math:`\delta_{CN}` is negative, we reduce C spent on N uptake and retain more C for growth

 .. math::

   \gamma_{FUN}  =
   \left\{\begin{array}{lr}
   \gamma_{FUN}+ 0.5.(delta_{CN}/c_{flexcn})& delta_{CN} > 0\\
   \gamma_{FUN}+(1-\gamma_{FUN}).min(1,\delta_{CN}/c_{flexcn}) &  delta_{CN} < 0
   \end{array}\right\}

We then restrict the degree to which C expenditure can be reduced (to prevent unrealistically high C:N ratios) as

 .. math::
   \gamma_{FUN} = max(min(1.0,\gamma_{FUN}),0.5)

Calculation of N uptake streams from active uptake and fixation
----------------------------------------------------------------

Once the final :math:`C_{nuptake}` is known, the fluxes of C to the individual pools can be derived as

 .. math::

   C_{nuptake,x}  = C_{frac,x}.C_{nuptake}

 .. math::

   N_{uptake,x}  = \frac{C_{nuptake}}{N_{cost}}

Following this, we determine whether the extraction estimates exceed the pool size for each source of N.  Where :math:`N_{active,no3} + N_{nonmyc,no3} > N_{avail,no3}`, we calculate the unmet uptake, :math:`N_{unmet,no3}`

 .. math::

   N_{unmet,no3}  = N_{active,no3} + N_{nonmyc,no3} - N_{avail,no3}

then modify both fluxes to account

 .. math::

   N_{active,no3} = N_{active,no3} +  N_{unmet,no3}.\frac{N_{active,no3}}{N_{active,no3}+N_{nonmyc,no3}}

 .. math::

   N_{nonmyc,no3} = N_{nonmyc,no3} +  N_{unmet,no3}.\frac{N_{nonmyc,no3}}{N_{active,no3}+N_{nonmyc,no3}}

and similarly, for NH4, where :math:`N_{active,nh4} + N_{nonmyc,nh4} > N_{avail,nh4}`, we calculate the unmet uptake, :math:`N_{unmet,no3}`

 .. math::

   N_{unmet,nh4}  = N_{active,nh4} + N_{nonmyc,nh4} - N_{avail,nh4}

then modify both fluxes to account

 .. math::

   N_{active,nh4} = N_{active,nh4} +  N_{unmet,nh4}.\frac{N_{active,nh4}}{N_{active,nh4}+N_{nonmyc,nh4}}

 .. math::

   N_{nonmyc,nh4} = N_{nonmyc,nh4} +  N_{unmet,nh4}.\frac{N_{nonmyc,nh4}}{N_{active,nh4}+N_{nonmyc,nh4}}

and then update the C spent to account for hte new lower N acquisition in that layer/pool.

 .. math::

   C_{active,nh4} = N_{active,nh4}.N_{cost,active,nh4}\\
   C_{active,no3} = N_{active,no3}.N_{cost,active,no3}\\
   C_{nonmyc,no3} = N_{nonmyc,no3}.N_{cost,nonmyc,no3}\\
   C_{nonmyc,no3} = N_{nonmyc,no3}.N_{cost,nonmyc,no3}\\

Following this, we determine how much carbon is accounted for for each soil layer.

 .. math::

   C_{accounted,x,j}  =  C_{spent,j,x} - (N_{acquired,j,x}.CN_{plant}.(1.0+ gr_{frac}))

