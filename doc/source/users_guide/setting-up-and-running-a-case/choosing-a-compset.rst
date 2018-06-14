.. _choosing-a-compset:

====================
 Choosing a compset
====================

When setting up a new case one of the first choices to make is which "component set" (or compset) to use. 
The compset refers to which component models are used as well as specific settings for them. 
We label the different types of compsets with a different letter of the alphabet from "A" (for all data model) to "X" (for all dead model). 
The compsets of interest when working with CLM are the "I" compsets (which contain CLM with a data atmosphere model and a stub ocean, and stub sea-ice models), "E" and "F" compsets (which contain CLM with the active atmosphere model (CAM), prescribed sea-ice model, and a data ocean model), and "B" compsets which have all active components. 
Below we go into details on the "I" compsets which emphasize CLM as the only active model, and just mention the two other categories.

When working with CLM you usually want to start with a relevant "I" compset before moving to the more complex cases that involve other active model components. 
The "I" compsets can exercise CLM in a way that is similar to the coupled modes, but with much lower computational cost and faster turnaround times.

Compsets coupled to data atmosphere and stub ocean/sea-ice ("I" compsets)
-------------------------------------------------------------------------

`Supported CLM Configurations <CLM-URL>`_ are listed in `Table 1-1 <CLM-1.1-Choosing-a-compset-using-CLM#table-1-1-scientifically-supported-i-compsets>`_ for the Scientifically Supported compsets (have been scientifically validated with long simulations) and in `Table 1-2 <CLM-1.1-Choosing-a-compset-using-CLM#table-1-2-functionally-supported-i-compsets>`_ for the Functionally Supported compsets (we've only checked that they function).


**Scientifically Supported I Compsets:**

+------------+--------------+-----------+-----------------+----------------+
| Short Name | Description  | Atm.      | Compset Alias   | Period         |
|            |              | Forcing   | Name            |                |
+============+==============+===========+=================+================+
| +|version|SP   | Satellite    | CRUNCEP   | 1850CRUCLM45    | 1850           |
|            | phenology    |           |                 |                |
|            | with new     |           |                 |                |
|            | biogeophys   |           |                 |                |
+------------+--------------+-----------+-----------------+----------------+
| +|version|SP   | New          | CRUNCEP   | I1850Clm50BgcCropCru| 1850           |
|            | biogeophys   |           +-----------------+----------------+
|            | + CENTURY-   |           | I20TRCRUCLM45BGC| 20th Century   |
|            | like         |           |                 |                |
|            | vertically   |           |                 |                |
|            | resolved     |           |                 |                |
|            | soil         |           |                 |                |
|            | BGC + CH4    |           |                 |                |
|            | emissions,   |           |                 |                |
|            | nitrogen     |           |                 |                |
|            | updates      |           |                 |                |
+------------+--------------+-----------+-----------------+----------------+
| +|version|CN   | New          | CRUNCEP   | I1850CRUCLM45CN | 1850           |
|            | biogeophys   |           |                 |                |
|            | + CN soil    |           |                 |                |
|            | BGC, updates |           |                 |                |
+------------+--------------+-----------+-----------------+----------------+
| CLM4SP     | As in        | Qian      | I1850           | 1850           |
|            | CCSM4/CESM1  |           +-----------------+----------------+
|            | release      |           | I               | 2000           |
|            |              |           +-----------------+----------------+
|            |              |           | I20TR           | 20th Century   |
+------------+--------------+-----------+-----------------+----------------+
| CLM4CN     | As in        | Qian      | I1850CN         | 1850           |
|            | CCSM4/CESM1  |           +-----------------+----------------+
|            | release      |           | ICN             | 2000           |
|            |              |           +-----------------+----------------+
|            |              |           | I20TRCN         | 20th Century   |
|            |              |           +-----------------+----------------+
|            |              |           | IRCP26CN        | RCP 2.6 to 2100|
|            |              |           +-----------------+----------------+
|            |              |           | IRCP45CN        | RCP 4.5 to 2100|
|            |              |           +-----------------+----------------+
|            |              |           | IRCP60CN        | RCP 6.0 to 2100|
|            |              |           +-----------------+----------------+
|            |              |           | IRCP85CN        | RCP 8.5 to 2100|
+------------+--------------+-----------+-----------------+----------------+

**Functionally Supported I Compsets:**

+------------+-------------------+-------------------+-----------------+--------+
| Short Name | Description       | Atm.              | Compset Alias   | Period |
|            |                   | Forcing           | Name            |        |
+============+===================+===================+=================+========+
| +|version|BGC- | ICRUCLM45BGCCROP  | New biogeophys +  | CRUNCEP         | 2000   |
| CROP       |                   | CENTURY-like      |                 |        |
|            |                   | vertically        |                 |        |
|            |                   | resolved soil     |                 |        |
|            |                   | BGC + CH4         |                 |        |
|            |                   | emissions,        |                 |        |
|            |                   | nitrogen updates  |                 |        |
|            |                   | with prognostic   |                 |        |
|            |                   | CROP              |                 |        |
|            |                   |                   |                 |        |
+------------+-------------------+-------------------+-----------------+--------+
| +|version|BGC- | I1850Clm50BgcCropCruDV| New biogeophys    | CRUNCEP         | 1850   |
| DV         |                   | + CENTURY-like    |                 |        |
|            |                   | vertically        |                 |        |
|            |                   | resolved soil     |                 |        |
|            |                   | BGC + CH4         |                 |        |
|            |                   | emissions,        |                 |        |
|            |                   | nitrogen updates  |                 |        |
|            |                   | with DV           |                 |        |
|            |                   |                   |                 |        |
|            |                   |                   |                 |        |
|            |                   |                   |                 |        |
+------------+-------------------+-------------------+-----------------+--------+
| +|version|SP-  | ICLM45VIC         | Satellite         | Qian            | 2000   |
| VIC        |                   | phenology with new|                 |        |
|            |                   | biogeophys with   |                 |        |
|            |                   | VIC hydrology     |                 |        |
+------------+-------------------+-------------------+-----------------+--------+
|CLM4CN-CROP | ICNCROP           | As in CCSM4/CESM1 | Qian            | 2000   |
|            |                   | release           |                 |        |
+------------+-------------------+-------------------+-----------------+--------+
|CLM4CN-DV   | ICNDV             | As in CCSM4/CESM1 | Qian            | 1850   |
|            |                   | release           |                 |        |
+------------+-------------------+-------------------+-----------------+--------+

Here is the entire list of compsets available. 
Note that using the "-user_compset" option even more combinations are possible. 
In the list below we give the alias name and then the long name which describes each component in parenthesis. 
Alias (Long-name with time-period and each component)

1. ``I`` (2000_DATM%QIA_CLM40%SP_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I1850`` (1850_DATM%QIA_CLM40%SP_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I1850CLM45`` (1850_DATM%QIA_CLM45%SP_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I1850CLM45BGC`` (1850_DATM%QIA_CLM45%BGC_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I1850CLM45CN`` (1850_DATM%QIA_CLM45%CN_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I1850CLM45CNF`` (1850_DATM%QIA_CLM45%CN_SICE_SOCN_RTM%FLOOD_SGLC_SWAV)

#. ``I1850CN`` (1850_DATM%QIA_CLM40%CN_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I1850CRU`` (1850_DATM%CRU_CLM40%SP_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I1850CRUCLM45`` (1850_DATM%CRU_CLM45%SP_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I1850Clm50BgcCropCru`` (1850_DATM%CRU_CLM45%BGC_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I1850Clm50BgcCropCruDV`` (1850_DATM%CRU_CLM45%BGCDV_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I1850CRUCLM45CN`` (1850_DATM%CRU_CLM45%CN_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I1850CRUCN`` (1850_DATM%CRU_CLM40%CN_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I1850SPINUPCLM45BGC`` (1850_DATM%S1850_CLM45%BGC_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I1850Clm50BgcSpinup`` (1850_DATM%S1850_CLM40%CN_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I1PT`` (2000_DATM%1PT_CLM40%SP_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I1PTCLM45`` (2000_DATM%1PT_CLM45%SP_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I20TR`` (20TR_DATM%QIA_CLM40%SP_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I20TRCLM45`` (20TR_DATM%QIA_CLM45%SP_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I20TRCLM45CN`` (20TR_DATM%QIA_CLM45%CN_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I20TRCN`` (20TR_DATM%QIA_CLM40%CN_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I20TRCRU`` (20TR_DATM%CRU_CLM40%SP_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I20TRCRUCLM45`` (20TR_DATM%CRU_CLM45%SP_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I20TRCRUCLM45BGC`` (20TR_DATM%CRU_CLM45%BGC_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I20TRCRUCLM45CN`` (20TR_DATM%CRU_CLM45%CN_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I20TRCRUCN`` (20TR_DATM%CRU_CLM40%CN_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I4804`` (4804_DATM%QIA_CLM40%SP_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I4804CLM45`` (4804_DATM%QIA_CLM45%SP_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I4804CLM45CN`` (4804_DATM%QIA_CLM45%CN_SICE_SOCN_RTM_SGLC_SWAV)

#. ``I4804CN`` (4804_DATM%QIA_CLM40%CN_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ICLM45`` (2000_DATM%QIA_CLM45%SP_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ICLM45BGC`` (2000_DATM%QIA_CLM45%BGC_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ICLM45BGCCROP`` (2000_DATM%QIA_CLM45%BGC-CROP_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ICLM45BGCDV`` (2000_DATM%QIA_CLM45%BGCDV_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ICLM45BGCDVCROP`` (2000_DATM%QIA_CLM45%BGCDV-CROP_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ICLM45BGCNoVS`` (2000_DATM%QIA_CLM45%NoVS_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ICLM45CN`` (2000_DATM%QIA_CLM45%CN_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ICLM45CNCROP`` (2000_DATM%QIA_CLM45%CN-CROP_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ICLM45CNDV`` (2000_DATM%QIA_CLM45%CNDV_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ICLM45CNTEST`` (2003_DATM%QIA_CLM45%CN_SICE_SOCN_RTM_SGLC_SWAV_TEST)

#. ``ICLM45CRUBGC`` (2000_DATM%CRU_CLM45%BGC_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ICLM45GLCMEC`` (2000_DATM%QIA_CLM45%CN_SICE_SOCN_RTM_CISM1_SWAV_TEST)

#. ``ICLM45SNCRFRC`` (2000_DATM%QIA_CLM45%SP-SNCR_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ICLM45USUMB`` (2000_DATM%1PT_CLM45%SP_SICE_SOCN_RTM_SGLC_SWAV_CLMUSRDAT%1x1_US-UMB)

#. ``ICLM45VIC`` (2000_DATM%QIA_CLM45%SP-VIC_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ICLM45alaskaCN`` (2000_DATM%QIA_CLM45%CN_SICE_SOCN_RTM_SGLC_SWAV_CLMUSRDAT%13x12pt_f19_alaskaUSA)

#. ``ICN`` (2000_DATM%QIA_CLM40%CN_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ICNCROP`` (2000_DATM%QIA_CLM40%CN-CROP_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ICNDV`` (2000_DATM%QIA_CLM40%CNDV_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ICNDVCROP`` (2000_DATM%QIA_CLM40%CNDV-CROP_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ICNTEST`` (2003_DATM%QIA_CLM40%CN_SICE_SOCN_RTM_SGLC_SWAV_TEST)

#. ``ICRU`` (2000_DATM%CRU_CLM40%SP_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ICRUCLM45`` (2000_DATM%CRU_CLM45_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ICRUCLM45BGC`` (2000_DATM%CRU_CLM45%BGC_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ICRUCLM45BGCCROP`` (2000_DATM%CRU_CLM45%BGC-CROP_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ICRUCLM45BGCTEST`` (2003_DATM%CRU_CLM45%BGC_SICE_SOCN_RTM_SGLC_SWAV_TEST)

#. ``ICRUCLM45CN`` (2000_DATM%CRU_CLM45%CN_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ICRUCN`` (2000_DATM%CRU_CLM40%CN_SICE_SOCN_RTM_SGLC_SWAV)

#. ``IG`` (2000_DATM%QIA_CLM40%SP_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IG1850`` (1850_DATM%QIA_CLM40%SP_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IG1850CLM45`` (1850_DATM%QIA_CLM45%SP_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IG1850CLM45CN`` (1850_DATM%QIA_CLM45%CN_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IG1850CN`` (1850_DATM%QIA_CLM40%CN_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IG20TR`` (20TR_DATM%QIA_CLM40%SP_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IG20TRCLM45`` (20TR_DATM%QIA_CLM45%SP_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IG20TRCLM45CN`` (20TR_DATM%QIA_CLM45%CN_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IG20TRCN`` (20TR_DATM%QIA_CLM40%CN_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IG4804`` (4804_DATM%QIA_CLM40%SP_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IG4804CLM45`` (4804_DATM%QIA_CLM45%SP_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IG4804CLM45CN`` (4804_DATM%QIA_CLM45%CN_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IG4804CN`` (4804_DATM%QIA_CLM40%CN_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IGCLM45`` (2000_DATM%QIA_CLM45%SP_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IGCLM45CN`` (2000_DATM%QIA_CLM45%CN_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IGCN`` (2000_DATM%QIA_CLM40%CN_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IGLCMEC`` (2000_DATM%QIA_CLM40%CN_SICE_SOCN_RTM_CISM1_SWAV_TEST)

#. ``IGRCP26CLM45CN`` (RCP2_DATM%QIA_CLM45%CN_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IGRCP26CN`` (RCP2_DATM%QIA_CLM40%CN_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IGRCP45CLM45CN`` (RCP4_DATM%QIA_CLM45%CN_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IGRCP45CN`` (RCP4_DATM%QIA_CLM40%CN_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IGRCP60CLM45CN`` (RCP6_DATM%QIA_CLM45%CN_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IGRCP60CN`` (RCP6_DATM%QIA_CLM40%CN_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IGRCP85CLM45CN`` (RCP8_DATM%QIA_CLM45%CN_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IGRCP85CN`` (RCP8_DATM%QIA_CLM40%CN_SICE_SOCN_RTM_CISM1_SWAV)

#. ``IRCP26CLM45CN`` (RCP2_DATM%QIA_CLM45%CN_SICE_SOCN_RTM_SGLC_SWAV)

#. ``IRCP26CN`` (RCP2_DATM%QIA_CLM40%CN_SICE_SOCN_RTM_SGLC_SWAV)

#. ``IRCP45CLM45CN`` (RCP4_DATM%QIA_CLM45%CN_SICE_SOCN_RTM_SGLC_SWAV)

#. ``IRCP45CN`` (RCP4_DATM%QIA_CLM40%CN_SICE_SOCN_RTM_SGLC_SWAV)

#. ``IRCP60CLM45CN`` (RCP6_DATM%QIA_CLM45%CN_SICE_SOCN_RTM_SGLC_SWAV)

#. ``IRCP60CN`` (RCP6_DATM%QIA_CLM40%CN_SICE_SOCN_RTM_SGLC_SWAV)

#. ``IRCP85CLM45CN`` (RCP8_DATM%QIA_CLM45%CN_SICE_SOCN_RTM_SGLC_SWAV)

#. ``IRCP85CN`` (RCP8_DATM%QIA_CLM40%CN_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ISNCRFRC`` (2000_DATM%QIA_CLM40%SP-SNCR_SICE_SOCN_RTM_SGLC_SWAV)

#. ``ITEST`` (2003_DATM%QIA_CLM40%SP_SICE_SOCN_RTM_SGLC_SWAV_TEST)

#. ``ITESTCLM45`` (2003_DATM%QIA_CLM45%SP_SICE_SOCN_RTM_SGLC_SWAV_TEST)

#. ``IUSUMB`` (2000_DATM%1PT_CLM40%SP_SICE_SOCN_RTM_SGLC_SWAV_CLMUSRDAT%1x1_US-UMB)

#. ``IalaskaCN`` (2000_DATM%QIA_CLM40%CN_SICE_SOCN_RTM_SGLC_SWAV_CLMUSRDAT%13x12pt_f19_alaskaUSA)

Compsets coupled to active atmosphere with data ocean
-----------------------------------------------------
CAM compsets are compsets that start with "E" or "F" in the name. They are described more fully in the scripts documentation or the CAM documentation. "E" compsets have a slab ocean model while "F" compsets have a data ocean model.

Fully coupled compsets with fully active ocean, sea-ice, and atmosphere
-----------------------------------------------------------------------
Fully coupled compsets are compsets that start with "B" in the name. They are described more fully in the scripts documentation.

Conclusion to choosing a compset
--------------------------------
We've introduced the basic type of compsets that use CLM and given some further details for the "standalone CLM" (or "I" compsets). 
The `config_compsets.xml <CLM-URL>`_ lists all of the compsets and gives a full description of each of them. 
In the next section we look into customizing the setup time options for compsets using CLM.
