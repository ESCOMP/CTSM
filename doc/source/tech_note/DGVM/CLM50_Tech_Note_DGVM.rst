.. _rst_Dynamic Global Vegetation and FATES:

Dynamic Global Vegetation and FATES
===================================

What has changed
^^^^^^^^^^^^^^^^^^^^

- Deprecation of the dynamic global vegetation model (DGVM): The CLM5.0 model contains the legacy 'CNDV' code, which runs the CLM biogeochemistry model in combination with the LPJ-derived dynamics vegetation model introduced in CLM3. While this capacity has not technically been removed from the model, the DGVM has not been tested in the development of CLM5 and is no longer scientifically supported.

- Introduction of FATES: The Functionally Assembled Terrestrial Ecosystem Simulator (FATES) is the actively developed DGVM for the CLM5.

FATES
^^^^^^^^^^^^^^^^^^^^

FATES is the "Functionally Assembled Terrestrial Ecosystem Simulator". It is an external module which can run within a given "Host Land Model" (HLM) like CLM.

FATES was derived from the CLM Ecosystem Demography model (CLM(ED)), which was documented in:

Fisher, R. A., Muszala, S., Verteinstein, M., Lawrence, P., Xu, C., McDowell, N. G., Knox, R. G., Koven, C., Holm, J., Rogers, B. M., Spessa, A., Lawrence, D., and Bonan, G.: Taking off the training wheels: the properties of a dynamic vegetation model without climate envelopes, CLM4.5(ED), Geosci. Model Dev., 8, 3593-3619, https://doi.org/10.5194/gmd-8-3593-2015, 2015.

The Ecosystem Demography ('ED'), concept within FATES is derived from the work of :ref:`Moorcroft et al. (2001)<mc_2001>` and is a cohort model of vegetation competition and co-existence, allowing a representation of the biosphere which accounts for the division of the land surface into successional stages, and for competition for light between height structured cohorts of representative trees of various plant functional types.

The implementation of the Ecosystem Demography concept within FATES links the surface flux and canopy physiology concepts in CLM with numerous additional developments necessary to accommodate the new model. These include a version of the SPITFIRE (Spread and InTensity of Fire) model of :ref:`Thonicke et al. (2010)<thonickeetal2010>`, and an adoption of the concept of `Perfect Plasticity Approximation` approach of :ref:`Purves et al. 2008<purves2008>`, :ref:`Lichstein et al. 2011<lichstein2011>` and :ref:`Weng et al. 2014<weng2014>`, in accounting for the spatial arrangement of crowns. Novel algorithms accounting for the fragmentation of coarse woody debris into chemical litter streams, for the physiological optimization of canopy thickness, for the accumulation of seeds in the seed bank, for multi-layer multi-PFT radiation transfer, for drought-deciduous and cold-deciduous phenology, for carbon storage allocation, and for tree mortality under carbon stress, are also included.

Further reading
^^^^^^^^^^^^^^^^^^^^

For more information about FATES, including a Users Guide and Technical Note, please see the `FATES documentation`_.

.. _FATES documentation: https://fates-users-guide.readthedocs.io/en/latest/index.html
