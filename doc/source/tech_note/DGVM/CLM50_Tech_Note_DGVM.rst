.. _rst_Dynamic Global Vegetation Model:

Dynamic Global Vegetation
===================================

What has changed
^^^^^^^^^^^^^^^^^^^^

- Deprecation of the dynamic global vegetation model (DGVM): The CLM5.0 model contains the legacy 'CNDV' code, which runs the CLM biogeochemistry model in combination with the LPJ-derived dynamics vegetation model introduced in CLM3. While this capacity has not technically been removed from the model, the DGVM has not been tested in the development of CLM5 and is no longer scientifically supported. 

- Introduction of FATES: The Functionally Assembled Terrestrial Ecosystem Simulator (FATES) is the actively developed DGVM for the CLM5. 


.. _rst_FATES:

Technical Documentation for FATES
===================================

FATES is the "Functionally Assembled Terrestrial Ecosystem Simulator". It is an external module which can run within a given "Host Land Model" (HLM). Currently (November 2017) implementations are supported in both the Community Land Model(CLM) and in the land model of the E3SM Dept. of Energy Earth System Model. 

FATES was derived from the CLM Ecosystem Demography model (CLM(ED)), which was documented in:

Fisher, R. A., Muszala, S., Verteinstein, M., Lawrence, P., Xu, C., McDowell, N. G., Knox, R. G., Koven, C., Holm, J., Rogers, B. M., Spessa, A., Lawrence, D., and Bonan, G.: Taking off the training wheels: the properties of a dynamic vegetation model without climate envelopes, CLM4.5(ED), Geosci. Model Dev., 8, 3593-3619, https://doi.org/10.5194/gmd-8-3593-2015, 2015.

and this technical note was first published as an appendix to that paper. 

https://pdfs.semanticscholar.org/396c/b9f172cb681421ed78325a2237bfb428eece.pdf

Introduction
^^^^^^^^^^^^^^^^^^^

The Ecosystem Demography ('ED'), concept within FATES is derived from the work of :ref:`Moorcroft et al. (2001)<mc_2001>`

and is a cohort model of vegetation competition and co-existence, allowing a representation of the biosphere which accounts for the division of the land surface into successional stages, and for competition for light between height structured cohorts of representative trees of various plant functional types. 

The implementation of the Ecosystem Demography
concept within FATES links the surface flux and canopy physiology concepts in the CLM/E3SM
with numerous additional developments necessary to accommodate the new
model also documented here. These include a version of the SPITFIRE
(Spread and InTensity of Fire) model of :ref:`Thonicke et al. (2010)<thonickeetal2010>`, and an adoption of the concept of
`Perfect Plasticity Approximation` approach of
:ref:`Purves et al. 2008<purves2008>`, :ref:`Lichstein et al. 2011<lichstein2011>` and :ref:`Weng et al. 2014<weng2014>`, in accounting
for the spatial arrangement of crowns. Novel algorithms accounting for
the fragmentation of coarse woody debris into chemical litter streams,
for the physiological optimisation of canopy thickness, for the
accumulation of seeds in the seed bank, for multi-layer multi-PFT
radiation transfer, for drought-deciduous and cold-deciduous phenology,
for carbon storage allocation, and for tree mortality under carbon
stress, are also included and presented here.

Numerous other implementations of the
Ecosystem Demography concept exist (See :ref:`Fisher et al. (2018)<Fisheretal2018>` for a review of these) Therefore, to avoid confusion between the
concept of 'Ecosystem Demography' and the implementation of this concept
in different models, the CLM(ED) implementation described by :ref:`Fisher et al. (2015)<Fisheretal2015>` will hereafter be called 'FATES' (the Functionally Assembled Terrestrial Ecosystem Simulator).

The representation of ecosystem heterogeneity in FATES
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The terrestrial surface of the Earth is heterogeneous for many reasons, driven
by variations in climate, edaphic history, ecological variability,
geological forcing and human interventions. Land surface models
represent this variability first by introducing a grid structure to the
land surface, allowing different atmospheric forcings to operate in each
grid cell, and subsequently by representing 'sub-grid' variability in
the surface properties. In the CLM, the land surface is divided into
numerous 'landunits' corresponding to the underlying condition of the
surface (e.g. soils, ice, lakes, bare ground) and then 'columns'
referring to elements of the surface that share below ground resources
(water & nutrients). Within the soil landunit, for example, there are
separate columns for crops, and for natural vegetation, as these are
assumed to use separate resource pools. The FATES model at present
only operates on the naturally vegetated column. The soil column is
sub-divided into numerous tiles, that correspond to statistical
fractions of the potentially vegetated land area. In the CLM 4.5 (and
all previous versions of the model), sub-grid tiling operates on the
basis of plant functional types (PFTs). That is, each piece of land is
assumed to be occupied by only one plant functional type, with multiple
PFT-specific tiles sharing a common soil water and nutrient pool. This
PFT-based tiling structure is the standard method used by most land
surface models deployed in climate prediction.

The introduction of the Ecosystem Demography concept introduces
significant alterations to the representation of the land surface in the
CLM. In FATES, the tiling structure represents the disturbance
history of the ecosystem. Thus, some fraction of the land surface is
characterized as 'recently disturbed', some fraction has escaped
disturbance for a long time, and other areas will have intermediate
disturbances. Thus the ED concept essentially discretizes the trajectory
of succession from disturbed ground to 'mature' ecosystems. Within
FATES, each "disturbance history class" is referred to as a ‘patch’.
The word "patch"  has many possible interpretations, so it is important
to note that: **there is no spatial location associated with the concept
of a 'patch' . It refers to a fraction of the potential vegetated area
consisting of all parts of the ecosystem with similar disturbance
history.**

The 'patch' organizational structure in CLM thus replaces the previous
'PFT' structure in the organization heirarchy. The original hierarchical
land surface organizational structure of CLM as described in
:ref:`Oleson et al. 2013<olesonetal2013>` may be depicted as:

.. math::

   \mathbf{gridcell} \left\{
   \begin{array}{cc} 
   \mathbf{landunit} &   \\ 
   \mathbf{landunit} &\left\{ 
   \begin{array}{ll} 
   \mathbf{column}&\\
   \mathbf{column}&\left\{ 
   \begin{array}{ll} 
   \mathbf{pft}&\\
   \mathbf{pft}&\\
   \mathbf{pft}&\\
   \end{array}\right.\\ 
   \mathbf{column}&\\
   \end{array}\right.\\ 
   \mathbf{landunit} &   \\
   \end{array}\right.

and the new structure is altered to the following:

.. math::

   \mathbf{gridcell} \left\{
   \begin{array}{cc} 
   \mathbf{landunit} &   \\ 
   \mathbf{landunit} &\left\{ 
   \begin{array}{ll} 
   \mathbf{column}&\\
   \mathbf{column}&\left\{ 
   \begin{array}{ll} 
   \mathbf{patch}&\\
   \mathbf{patch}&\\
   \mathbf{patch}&\\
   \end{array}\right.\\ 
   \mathbf{column}&\\
   \end{array}\right.\\ 
   \mathbf{landunit} &   \\
   \end{array}\right.

Thus, each gridcell becomes a matrix of 'patches' that are
conceptualized by their 'age since disturbance' in years. This is the
equivalent of grouping together all those areas of a gridcell that are
'canopy gaps', into a single entity, and all those areas that are
'mature forest' into a single entity.

Cohortized representation of tree populations
---------------------------------------------

Each common-disturbance-history patch is a notional ecosystem that might
in reality contain numerous individual plants which vary in their
physiological attributes, in height and in spatial position. One way of
addressing this heterogeneity is to simulate a forest of specific
individuals, and to monitor their behavior through time. This is the
approach taken by "gap" and individual-based models 
(:ref:`Smith et al. 2001<smith2001>`, :ref:`Sato et al. 2007<sato2007>`, :ref:`Uriarte et al. 2009<uriarte2009>`, :ref:`Fyllas et al. 2014 <fyllas2014>`). The
depiction of individuals typically implies that the outcome of the model
is stochastic. This is because we lack the necessary detailed knowledge
to simulate the individual plant's fates. Thus gap models imply both
stochastic locations and mortality of plants. Thus, (with a genuinely
random seed) each model outcome is different, and an ensemble of model
runs is required to generate an average representative solution. Because
the random death of large individual trees can cause significant
deviations from the mean trajectory for a small plot (a typical
simulated plot size is 30m x 30 m) the number of runs required to
minimize these deviations is large and computationally expensive. For
this reason, models that resolve individual trees typically use a
physiological timestep of one day or longer (e.g.  :ref:`Smith et al. 2001<smith2001>`, :ref:`Xiaidong et al. 2005 <xiaodong2005>`, :ref:`Sato et al. 2007<sato2007>`

The approach introduced by the Ecosystem Demography model
:ref:`Moorcroft et al. 2001<mc_2001>` is to group the hypothetical population
of plants into "cohorts". In the notional ecosystem, after the
land-surface is divided into common-disturbance-history patches, the
population in each patch is divided first into plant functional types
(the standard approach to representing plant diversity in large scale
vegetation models), and then each plant type is represented as numerous
height classes. Importantly, **for each PFT/height class bin, we model
*one* representative individual plant, which tracks the average
properties of this `cohort` of individual plants.** Thus, each
common-disturbance-history patch is typically occupied by a set of
cohorts of different plant functional types, and different height
classes within those plant functional types. Each cohort is associated
with a number of identical trees, :math:`n_{coh}` (where :math:`{coh}`
denotes the identification or index number for a given cohort)..

The complete hierarchy of elements in FATES is therefore now
described as follows:

.. math::

   \mathbf{gridcell}\left\{
   \begin{array}{cc} 
   \mathbf{landunit} &   \\ 
   \mathbf{landunit} &\left\{ 
   \begin{array}{ll} 
   \mathbf{column}&\\
   \mathbf{column}&\left\{ 
   \begin{array}{ll} 
   \mathbf{patch}&\\
   \mathbf{patch}&\left\{ 
   \begin{array}{ll} 
   \mathbf{cohort}&\\
   \mathbf{cohort}&\\
   \mathbf{cohort}&\\
   \end{array}\right.\\ 
   \mathbf{patch}&\\
   \end{array}\right.\\ 
   \mathbf{column}&\\
   \end{array}\right.\\ 
   \mathbf{landunit} &   \\
   \end{array}\right.

Discretization of cohorts and patches
-------------------------------------

Newly disturbed land and newly recruited seedlings can in theory be
generated at each new model timestep as the result of germination and
disturbance processes. If the new patches and cohorts established at
*every* timestep were tracked by the model structure, the computational
load would of course be extremely high (and thus equivalent to an
individual-based approach). A signature feature of the ED model is the
system by which `functionally equivalent` patches and cohorts are fused
into single model entities to save memory and computational time.

This functionality requires that criteria are established for the
meaning of `functional equivalence`, which are by necessity slightly
subjective, as they represent ways of abstracting reality into a more
tractable mathematical representation. As an example of this, for
height-structured cohorts, we calculate the relativized differences in
height (:math:`h_{coh}`, m) between two cohorts of the same pft,
:math:`p` and :math:`q` as

.. math:: d_{hite,p,q} = \frac{\mathrm{abs}.(h_{p-}h_{q})}{\frac{1}{2}(h_{p}+h_{q})}

If :math:`d_{hite,p,q}` is smaller than some threshold :math:`t_{ch}`,
and they are of the same plant functional type, the two cohorts are
considered equivalent and merged to form a third cohort :math:`r`, with
the properties of cohort :math:`p` and :math:`q` averaged such that they
conserve mass. The model parameter :math:`t_{ch}` can be adjusted to
adjust the trade-off between simulation accuracy and computational load.
There is no theoretical optimal value for this threshold but it may be
altered to have finer or coarser model resolutions as needed.

Similarly, for common-disturbance-history patches, we again assign
a threshold criteria, which is then compared to the difference between
patches :math:`m` and :math:`n`, and if the difference is less than some
threshold value (:math:`t_{p}`) then patches are merged together,
otherwise they are kept separate. However, in contrast with
height-structured cohorts, where the meaning of the difference criteria
is relatively clear, how the landscape should be divided into
common-disturbance-history units is less clear. Several alternative
criteria are possible, including Leaf Area Index, total biomass and
total stem basal area.

In this implementation of FATES we assess the amount of
above-ground biomass in each PFT/plant diameter bin. Biomass is first
grouped into fixed diameter bins for each PFT (:math:`ft`) and a
significant difference in any bin will cause patches to remain
separated. This means that if two patches have similar total biomass,
but differ in the distribution of that biomass between diameter classes
or plant types, they remain as separate entities. Thus

.. math:: B_{profile,m,dc,ft} = \sum_{d_{c,min}}^{d_{c,max}} (B_{ag,coh}n_{coh})

:math:`B_{profile,m,dc,ft}` is the binned above-ground biomass profile
for patch :math:`m`,\ :math:`d_{c}` is the diameter class.
:math:`d_{c,min}` and :math:`d_{c,max}` are the lower and upper
boundaries for the :math:`d_{c}` diameter class. :math:`B_{ag,coh}` and
:math:`n_{coh}` depict the biomass (KgC m\ :math:`^{-2}`) and the number
of individuals of each cohort respectively. A difference matrix between
patches :math:`m` and :math:`n` is thus calculated as

.. math:: d_{biomass,mn,dc,ft} = \frac{\rm{abs}\it(B_{profile,m,hc,ft}-B_{profile,n,hc,ft})}{\frac{1}{2}(B_{profile,m,hc,ft}+B_{profile,n,hc,ft})}

If all the values of :math:`d_{biomass,mn,hc,ft}` are smaller than the
threshold, :math:`t_{p}`, then the patches :math:`m` and :math:`n` are
fused together to form a new patch :math:`o`.

To increase computational efficiency and to simplify the coding
structure of the model, the maximum number of patches is capped at
:math:`P_{no,max}`. To force the fusion of patches down to this number,
the simulation begins with a relatively sensitive discretization of
patches (:math:`t_{p}` = 0.2) but if the patch number exceeds the
maximum, the fusion routine is repeated iteratively until the two most
similar patches reach their fusion threshold. This approach maintains an
even discretization along the biomass gradient, in contrast to, for
example, simply fusing the oldest or youngest patches together.

The area of the new patch (:math:`A_{patch,o}`, m\ :math:`^{2}`)
is the sum of the area of the two existing patches,

.. math:: A_{patch,o}  = A_{patch,n}  + A_{patch,m}

and the cohorts ‘belonging’ to patches :math:`m` and :math:`n` now
co-occupy patch :math:`o`. The state properties of :math:`m` and
:math:`n` (litter, seed pools, etc. ) are also averaged in accordance
with mass conservation .

Linked Lists: the general code structure of FATES
---------------------------------------------------

The number of patches in each natural vegetation column and the
number of cohorts in any given patch are variable through time because
they are re-calculated for each daily timestep of the model. The more
complex an ecosystem, the larger the number of patches and cohorts. For
a slowly growing ecosystem, where maximum cohort size achieved between
disturbance intervals is low, the number of cohorts is also low. For
fast-growing ecosystems where many plant types are viable and maximum
heights are large, more cohorts are required to represent the ecosystem
with adequate complexity.

In terms of variable structure, the creation of an array whose size
could accommodate every possible cohort would mean defining the maximum
potential number of cohorts for every potential patch, which would
result in very large amounts of wasted allocated memory, on account of
the heterogeneity in the number of cohorts between complex and simple
ecosystems (n.b. this does still happen for some variables at restart
timesteps). To resolve this, the cohort structure in FATES model
does not use an array system for internal calculations. Instead it uses
a system of *linked lists* where each cohort structure is linked to the
cohorts larger than and smaller than itself using a system of pointers.
The shortest cohort in each patch has a ‘shorter’ pointer that points to
the *null* value, and the tallest cohort has a ‘taller’ pointer that
points to the null value. 

Instead of iterating along a vector indexed by :math:`coh`, the code
structures typically begin at the tallest cohort in a given patch, and
iterate until a null pointer is encountered. 

Using this structure, it is therefore possible to have an unbounded upper limit on cohort number, and also to easily alter the ordering of  cohorts if, for example, a cohort of one functional type begins to  grow faster than a competitor of another functional type, and the cohort list can easily be re-ordered by altering the pointer structure. Each cohort has `pointers` indicating to which patch and gridcell it belongs. The patch system is analogous to the cohort system, except that patches are ordered in terms of their relative age, with pointers to older and younger patches where cp\ :math:`_1` is the oldest:


Indices used in FATES
-----------------------

Some of the indices used in FATES are similar to those used in the
standard CLM4.5 model; column (:math:`c`), land unit(\ :math:`l`), grid
cell(\ :math:`g`) and soil layer (:math:`j`). On account of the
additional complexity of the new representation of plant function,
several additional indices are introduced that describe the
discritization of plant type, fuel type, litter type, plant height,
canopy identity, leaf vertical structure and fuel moisture
characteristics. To provide a reference with which to interpret the
equations that follow, they are listed here.

.. raw:: latex

   \bigskip

.. raw:: latex

   \captionof{table}{Table of subscripts used in this document  }

+------------------+-----------------------+
| Parameter Symbol | Parameter Name        |
+==================+=======================+
| *ft*             | Plant Functional Type |
+------------------+-----------------------+
| *fc*             | Fuel Class            |
+------------------+-----------------------+
| *lsc*            | Litter Size Class     |
+------------------+-----------------------+
| *coh*            | Cohort Index          |
+------------------+-----------------------+
| *patch*          | Patch Index           |
+------------------+-----------------------+
| *Cl*             | Canopy Layer          |
+------------------+-----------------------+
| *z*              | Leaf Layer            |
+------------------+-----------------------+
| *mc*             | Moisture Class        |
+------------------+-----------------------+

.. raw:: latex

   \bigskip 

Cohort State Variables
----------------------

The unit of allometry in the ED model is the cohort. Each cohort
represents a group of plants with similar functional types and heights
that occupy portions of column with similar disturbance histories. The
state variables of each cohort therefore consist of several pieces of
information that fully describe the growth status of the plant and its
position in the ecosystem structure, and from which the model can be
restarted. The state variables of a cohort are as follows:

.. raw:: latex

   \bigskip

.. raw:: latex

   \captionof{table}{State Variables of  `cohort' sructure}

+-----------------+-----------------+-----------------+-----------------+
| Quantity        | Variable name   | Units           | Notes           |
+=================+=================+=================+=================+
| Plant           | :math:`{\it{ft} | integer         |                 |
| Functional Type | _{coh}}`        |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| Number of       | :math:`n_{coh}` | n per           |                 |
| Individuals     |                 | 10000m\ :math:` |                 |
|                 |                 | ^{-2}`          |                 |
+-----------------+-----------------+-----------------+-----------------+
| Height          | :math:`h_{coh}` | m               |                 |
+-----------------+-----------------+-----------------+-----------------+
| Diameter        | :math:`\it{dbh_ | cm              |                 |
|                 | {coh}}`         |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| Structural      | :math:`{b_{stru | KgC             | Stem wood       |
| Biomass         | c,coh}}`        | plant\ :math:`^ | (above and      |
|                 |                 | {-1}`           | below ground)   |
+-----------------+-----------------+-----------------+-----------------+
| Alive Biomass   | :math:`{b_{aliv | KgC             | Leaf, fine root |
|                 | e,coh}}`        | plant\ :math:`^ | and sapwood     |
|                 |                 | {-1}`           |                 |
+-----------------+-----------------+-----------------+-----------------+
| Stored Biomass  | :math:`{b_{stor | KgC             | Labile carbon   |
|                 | e,coh}}`        | plant\ :math:`^ | reserve         |
|                 |                 | {-1}`           |                 |
+-----------------+-----------------+-----------------+-----------------+
| Leaf memory     | :math:`{l_{memo | KgC             | Leaf mass when  |
|                 | ry,coh}}`       | plant\ :math:`^ | leaves are      |
|                 |                 | {-1}`           | dropped         |
+-----------------+-----------------+-----------------+-----------------+
| Canopy Layer    | :math:`{C_{l,co | integer         | 1 = top layer   |
|                 | h}}`            |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| Phenological    | :math:`{S_{phen | integer         | 1=leaves off.   |
| Status          | ,coh}}`         |                 | 2=leaves on     |
+-----------------+-----------------+-----------------+-----------------+
| Canopy trimming | :math:`C_{trim, | fraction        | 1.0=max leaf    |
|                 | coh}`           |                 | area            |
+-----------------+-----------------+-----------------+-----------------+
| Patch Index     | :math:`{p_{coh} | integer         | To which patch  |
|                 | }`              |                 | does this       |
|                 |                 |                 | cohort belong?  |
+-----------------+-----------------+-----------------+-----------------+

Patch State Variables
---------------------

A patch, as discuss earlier, is a fraction of the landscape which
contains ecosystems with similar structure and disturbance history. A
patch has no spatial location. The state variables, which are
‘ecosystem’ rather than ‘tree’ scale properties, from which the model
can be restarted, are as follows

.. raw:: latex

   \bigskip

.. raw:: latex

   \captionof{table}{State variables of `patch' structure}

+-------------+-------------+-------------+-------------+-------------+
| Quantity    | Variable    | Units       | Indexed By  |             |
|             | name        |             |             |             |
+=============+=============+=============+=============+=============+
| Area        | :math:`\it{ | m\ :math:`^ | -           |             |
|             | A_{patch}}` | {2}`        |             |             |
+-------------+-------------+-------------+-------------+-------------+
| Age         | :math:`age_ | years       | -           |             |
|             | {patch}`    |             |             |             |
+-------------+-------------+-------------+-------------+-------------+
| Seed        | :math:`seed_| KgC         | :math:`ft`  |             |
|             | {patch}`    | m\ :math:`^ |             |             |
|             |             | {-2}`       |             |             |
+-------------+-------------+-------------+-------------+-------------+
| Leaf Litter | :math:`l_{l | KgC         | :math:`ft`  |             |
|             | itter,patch | m\ :math:`^ |             |             |
|             | }`          | {-2}`       |             |             |
+-------------+-------------+-------------+-------------+-------------+
| Root Litter | :math:`r_{l | KgC         | :math:`ft`  |             |
|             | itter,patch | m\ :math:`^ |             |             |
|             | }`          | {-2}`       |             |             |
+-------------+-------------+-------------+-------------+-------------+
| AG Coarse   | :math:`     | KgC         | Size Class  |             |
| Woody       | {CWD}_{A    | m\ :math:`^ | (lsc)       |             |
| Debris      | G,patch}`   | {-2}`       |             |             |
+-------------+-------------+-------------+-------------+-------------+
| BG Coarse   | :math:`     | KgC         | Size Class  |             |
| Woody       | {CWD}_{B    | m\ :math:`^ | (lsc)       |             |
| Debris      | G,patch}`   | {-2}`       |             |             |
+-------------+-------------+-------------+-------------+-------------+
| Canopy      | :math:`S_{c | -           | Canopy      |             |
| Spread      | ,patch}`    |             | Layer       |             |
+-------------+-------------+-------------+-------------+-------------+
| Column      | :math:`{l_{ | integer     | -           |             |
| Index       | patch}}`    |             |             |             |
+-------------+-------------+-------------+-------------+-------------+


Model Structure
---------------

Code concerned with the Ecosystem Demography model interfaces with the
CLM model in four ways: i) During initialization, ii) During the
calculation of surface processes (albedo, radiation absorption, canopy
fluxes) each model time step (typically half-hourly), iii) During the
main invokation of the ED model code at the end of each day. Daily
cohort-level NPP is used to grow plants and alter the cohort structures,
disturbance processes (fire and mortality) operate to alter the patch
structures, and all fragmenting carbon pool dynamics are calculated. iv)
during restart reading and writing. The net assimilation (NPP) fluxes
attributed to each cohort are accumulated throughout each daily cycle
and passed into the ED code as the major driver of vegetation dynamics.

Initialization of vegetation from bare ground
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the model is restarted from a bare ground state (as opposed to a
pre-existing vegetation state), the state variables above are
initialized as follows. First, the number of plants per PFT is allocated
according to the initial seeding density (:math:`S_{init}`, individuals
per m\ :math:`^{2}`) and the area of the patch :math:`A_{patch}`, which
in the first timestep is the same as the area of the notional ecosystem
:math:`A_{tot}`. The model has no meaningful spatial dimension, but we
assign a notional area such that the values of ‘:math:`n_{coh}`’ can be
attributed. The default value of :math:`A_{tot}` is one hectare (10,000
m\ :math:`^{2}`), but the model will behave identically irrepective of
the value of this parameter.

.. math:: n_{coh,0} = S_{init}A_{patch}

Each cohort is initialized at the minimum canopy height
:math:`h_{min,ft}`, which is specified as a parameter for each plant
functional type and denotes the smallest size of plant which is tracked
by the model. Smaller plants are not considered, and their emergence
from the recruitment processes is unresolved and therefore implicitly
parameterized in the seedling establishment model.. The diameter of each
cohort is then specified using the log-linear allometry between stem
diameter and canopy height

.. math:: \mathit{dbh}_{coh} = 10^{\frac{\log_{10}(h_{coh}) - c_{allom}}{m_{allom}}  }

where the slope of the log-log relationship, :math:`m_{allom}` is 0.64
and the intercept :math:`c_{allom}` is 0.37. The structural biomass
associated with a plant of this diameter and height is given (as a
function of wood density, :math:`\rho`, g cm\ :math:`^{-3}`)

.. math:: b_{struc,coh} =c_{str}h_{coh}^{e_{str,hite}} dbh_{coh}^{e_{str,dbh}} \rho_{ft}^{e_{str,dens}}

taken from the original ED1.0 allometry
:ref:`Moorcroft et al. 2001<mc_2001>` (values of the allometric constants in
Table `[table:allom] <#table:allom>`__. The maximum amount of leaf
biomass associated with this diameter of tree is calculated according to
the following allometry

.. math:: b_{max,leaf,coh} =c_{leaf}\it{dbh}_{coh}^{e_{leaf,dbh}} \rho_{ft}^{e_{leaf,dens}}

from this quantity, we calculate the active/fine root biomass
:math:`b_{root,coh}` as

.. math:: b_{root,coh} =  b_{max,leaf,coh}\cdot f_{frla}

where :math:`f_{frla}` is the fraction of fine root biomass to leaf
biomass, assigned per PFT

.. raw:: latex

   \captionof{table}{Parameters needed for model initialization.}

+-----------------+-----------------+-----------------+-----------------+
| Parameter       | Parameter Name  | Units           | Default Value   |
| Symbol          |                 |                 |                 |
+=================+=================+=================+=================+
| :math:`h_{min}` | Minimum plant   | m               | 1.5             |
|                 | height          |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`S_{init}`| Initial         | Individuals     |                 |
|                 | Planting        | m\ :math:`^{-2}`|                 |
|                 | density         |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`A_{tot}` | Model area      | m\ :math:`^{2}` | 10,000          |
+-----------------+-----------------+-----------------+-----------------+

[table:init]

Allocation of biomass
^^^^^^^^^^^^^^^^^^^^^

Total live biomass :math:`b_{alive}` is the state variable of the model
that describes the sum of the three live biomass pools leaf
:math:`b_{leaf}`, root :math:`b_{root}` and sapwood :math:`b_{sw}` (all
in kGC individual\ :math:`^{-1}`). The quantities are constrained by the
following

.. math:: b_{alive} = b_{leaf} + b_{root} + b_{sw}

Sapwood volume is a function of tree height and leaf biomass

.. math:: b_{sw} = b_{leaf}\cdot h_{coh}\cdot f_{swh}

where :math:`f_{swh}` is the ratio of sapwood mass (kgC) to leaf mass
per unit tree height (m). Also, root mass is a function of leaf mass

.. math:: b_{root} = b_{leaf}\cdot f_{swh}

Thus

.. math:: b_{alive} = b_{leaf} + b_{leaf}\cdot f_{frla} + b_{leaf}\cdot h_{coh}\cdot f_{swh}

Rearranging gives the fraction of biomass in the leaf pool
:math:`f_{leaf}` as

.. math:: f_{leaf} = \frac{1}{1+h_{coh}\cdot f_{swh}+f_{frla} }

Thus, we can determine the leaf fraction from the height at the tissue
ratios, and the phenological status of the cohort :math:`S_{phen,coh}`.

.. math:: b_{leaf} = b_{alive} \cdot l _{frac}

To divide the live biomass pool at restart, or whenever it is
recalculated, into its consituent parts, we first

.. math::

   b_{leaf} = \left\{ \begin{array}{ll}
   b_{alive} \cdot l _{frac}&\textrm{for } S_{phen,coh} = 1\\
   &\\
   0&\textrm{for } S_{phen,coh} = 0\\
   \end{array} \right.

Because sometimes the leaves are dropped, using leaf biomass as a
predictor of root and sapwood would produce zero live biomass in the
winter. To account for this, we add the LAI memory variable
:math:`l_{memory}` to the live biomass pool to account for the need to
maintain root biomass when leaf biomass is zero. Thus, to calculated the
root biomass, we use

.. math:: b_{root} = (b_{alive}+l_{memory})\cdot l_{frac} \cdot f_{frla}

To calculated the sapwood biomass, we use

.. math:: b_{sw} = (b_{alive}+l_{memory})\cdot l_{frac} \cdot f_{swh} \cdot h_{coh}

.. raw:: latex

   \captionof{table}{Allometric Constants}

+-----------------+-----------------+-----------------+-----------------+
| Parameter       | Parameter Name  | Units           | Default Value   |
| Symbol          |                 |                 |                 |
+=================+=================+=================+=================+
| :math:`c_{allom | Allometry       |                 | 0.37            |
| }`              | intercept       |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`m_{allom | Allometry slope |                 | 0.64            |
| }`              |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`c_{str}` | Structural      |                 | 0.06896         |
|                 | biomass         |                 |                 |
|                 | multiplier      |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`e_{str,d | Structural      |                 | 1.94            |
| bh}`            | Biomass dbh     |                 |                 |
|                 | exponent        |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`e_{str,h | Structural      |                 | 0.572           |
| ite}`           | Biomass height  |                 |                 |
|                 | exponent        |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`e_{str,d | Structural      |                 | 0.931           |
| ens}`           | Biomass density |                 |                 |
|                 | exponent        |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`c_{leaf}`| Leaf biomass    |                 | 0.0419          |
|                 | multiplier      |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`e_{leaf, | Leaf biomass    |                 | 1.56            |
| dbh}`           | dbh exponent    |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`e_{leaf, | Leaf biomass    |                 | 0.55            |
| dens}`          | density         |                 |                 |
|                 | exponent        |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`f_{swh}` | Ratio of        | m\ :math:`^{-1}`|                 |
|                 | sapwood mass to |                 |                 |
|                 | height          |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`f_{frla}`| Ratio of fine   | -               | 1.0             |
|                 | root mass to    |                 |                 |
|                 | leaf mass       |                 |                 |
+-----------------+-----------------+-----------------+-----------------+

[table:allom]

Canopy Structure and the Perfect Plasticity Approximation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

During initialization and every subsequent daily ED timestep, the canopy
structure model is called to determine how the leaf area of the
different cohorts is arranged relative to the incoming radiation, which
will then be used to drive the radiation and photosynthesis
calculations. This task requires that some assumptions are made about 1)
the shape and depth of the canopy within which the plant leaves are
arranged and 2) how the leaves of different cohorts are arranged
relative to each other. This set of assumptions are critical to model
performance in ED-like cohort based models, since they determine how
light resources are partitioned between competing plants of varying
heights, which has a very significant impact on how vegetation
distribution emerges from competition
:ref:`Fisher et al. 2010<Fisheretal2010>`.

The standard ED1.0 model makes a simple 'flat disk' assumption, that the
leaf area of each cohort is spread in an homogenous layer at one exact
height across entire the ground area represented by each patch. FATES has diverged from this representation due to (at least) two problematic emergent properties that we identified as generating unrealistic behaviours espetially for large-area patches.

1. Over-estimation of light competition . The vertical stacking of
cohorts which have all their leaf area at the same nominal height means
that when one cohort is only very slightly taller than it’s competitor,
it is completely shaded by it. This means that any small advantage in
terms of height growth translates into a large advantage in terms of
light competition, even at the seedling stage. This property of the
model artificially exaggerates the process of light competition. In
reality, trees do not compete for light until their canopies begin to
overlap and canopy closure is approached.

2. Unrealistic over-crowding. The 'flat-disk' assumption has no
consideration of the spatial extent of tree crowns. Therefore it has no
control on the packing density of plants in the model. Given a mismatch
between production and mortality, entirely unrealistic tree densities
are thus possible for some combinations of recruitment, growth and
mortality rates.

To account for the filling of space in three dimensions using the
one-dimensional representation of the canopy employed by CLM, we
implement a new scheme derived from that of
:ref:`Purves et al. 2008<purves2008>`. Their argument follows the development
of an individual-based variant of the SORTIE model, called SHELL, which
allows the location of individual plant crowns to be highly flexible in
space. Ultimately, the solutions of this model possess an emergent
property whereby the crowns of the plants simply fill all of the
available space in the canopy before forming a distinct understorey.

Purves et al. developed a model that uses this feature, called the
‘perfect plasticity approximation’, which assumes the plants are able to
perfectly fill all of the available canopy space. That is, at canopy
closure, all of the available horizontal space is filled, with
negligible gaps, owing to lateral tree growth and the ability of tree
canopies to grow into the available gaps (this is of course, an
over-simplified but potential useful ecosystem property). The ‘perfect
plasticity approximation’ (PPA) implies that the community of trees is
subdivided into discrete canopy layers, and by extension, each cohort
represented by FATES model is assigned a canopy layer status flag,
:math:`C_L`. In this version, we set the maximum number of canopy layers
at 2 for simplicity, although is possible to have a larger number of
layers in theory. :math:`C_{L,coh}` = 1 means that all the trees of
cohort :math:`coh` are in the upper canopy (overstory), and
:math:`C_{L,coh}` = 2 means that all the trees of cohort :math:`coh` are
in the understorey.

In this model, all the trees in the canopy experience full light on
their uppermost leaf layer, and all trees in the understorey experience
the same light (full sunlight attenuated by the average LAI of the upper
canopy) on their uppermost leaves, as described in the radiation
transfer section (more nuanced versions of this approach may be
investigated in future model versions). The canopy is assumed to be
cylindrical, the lower layers of which experience self-shading by the
upper layers.

To determine whether a second canopy layer is required, the model needs
to know the spatial extent of tree crowns. Crown area,
:math:`A_{crown}`, m\ :math:`^{2}`, is defined as

.. math:: A_{crown,coh}  = \pi (dbh_{coh} S_{c,patch,Cl})^{1.56}

where :math:`A_{crown}` is the crown area of a single tree canopy
(m:math:`^{-2}`) and :math:`S_{c,patch,Cl}` is the ‘canopy spread’
parameter (m cm^-1) of this canopy layer, which is assigned as a
function of canopy space filling, discussed below. In contrast to
:ref:`Purves et al. 2008<purves2008>` , we use an exponent, identical to that
for leaf biomass, of 1.56, not 2.0, such that tree leaf area index does
not change as a function of diameter.

To determine whether the canopy is closed, we calculate the total canopy
area as:

.. math:: A_{canopy} = \sum_{coh=1}^{nc,patch}{A_{crown,coh}.n_{coh}}

where :math:`nc_{patch}` is the number of cohorts in a given patch. If
the area of all crowns :math:`A_{canopy}` (m:math:`^{-2}`) is larger
than the total ground area of a patch (:math:`A_{patch}`), then some
fraction of each cohort is demoted to the understorey.

Under these circumstances, the `extra` crown area :math:`A_{loss}`
(i.e., :math:`A_{canopy}` - :math:`A_p`) is moved into the understorey.
For each cohort already in the canopy, we determine a fraction of trees
that are moved from the canopy (:math:`L_c`) to the understorey.
:math:`L_c` is calculated as :ref:`Fisher et al. 2010<Fisheretal2010>`

.. math:: L_{c}= \frac{A_{loss,patch} w_{coh}}{\sum_{coh=1}^{nc,patch}{w_{coh}}} ,

where :math:`w_{coh}` is a weighting of each cohort determined by basal
diameter :math:`dbh` (cm) and the competitive exclusion coefficient
:math:`C_{e}`

.. math:: w_{coh}=dbh_{coh}C_{e}.

The higher the value of :math:`C_e` the greater the impact of tree
diameter on the probability of a given tree obtaining a position in the
canopy layer. That is, for high :math:`C_e` values, competition is
highly deterministic. The smaller the value of :math:`C_e`, the greater
the influence of random factors on the competitive exclusion process,
and the higher the probability that slower growing trees will get into
the canopy. Appropriate values of :math:`C_e` are poorly constrained but
alter the outcome of competitive processes.

The process by which trees are moved between canopy layers is complex
because 1) the crown area predicted for a cohort to lose may be larger
than the total crown area of the cohort, which requires iterative
solutions, and 2) on some occasions (e.g. after fire), the canopy may
open up and require ‘promotion’ of cohorts from the understorey, and 3)
canopy area may change due to the variations of canopy spread values (
:math:`S_{c,patch,Cl}`, see the section below for details) when
fractions of cohorts are demoted or promoted. Further details can be
found in the code references in the footnote.

Horizontal Canopy Spread
-------------------------

:ref:`Purves et al. 2008<purves2008>` estimated the ratio between canopy and
stem diameter :math:`c_{p}` as 0.1 m cm\ :math:`^{-1}` for canopy trees
in North American forests, but this estimate was made on trees in closed
canopies, whose shape is subject to space competition from other
individuals. Sapling trees have no constraints in their horizontal
spatial structure, and as such, are more likely to display their leaves
to full sunlight. Also, prior to canopy closure, light interception by
leaves on the sides of the canopy is also higher than it would be in a
closed canopy forest. If the ‘canopy spread’ parameter is constant for
all trees, then we simulate high levels of self-shading for plants in
unclosed canopies, which is arguably unrealistic and can lower the
productivity of trees in areas of unclosed canopy (e.g. low productivity
areas of boreal or semi-arid regions where LAI and canopy cover might
naturally be low). We here interpret the degree of canopy spread,
:math:`S_{c}` as a function of how much tree crowns interfere with each
other in space, or the total canopy area :math:`A_{canopy}`. However
:math:`A_{canopy}` itself is a function of :math:`S_{c}`, leading to a
circularity. :math:`S_{c}` is thus solved iteratively through time.

Each daily model step, :math:`A_{canopy}` and the fraction of the
gridcell occupied by tree canopies in each canopy layer
(:math:`A_{f,Cl}` = :math:`A_{canopy,Cl}`/:math:`A_{patch}`) is
calculated based on :math:`S_{c}` from the previous timestep. If
:math:`A_{f}` is greater than a threshold value :math:`A_{t}`,
:math:`S_{c}` is increased by a small increment :math:`i`. The threshold
:math:`A_{t}` is, hypothetically, the canopy fraction at which light
competition begins to impact on tree growth. This is less than 1.0 owing
to the non-perfect spatial spacing of tree canopies. If :math:`A_{f,Cl}`
is greater than :math:`A_{t}`, then :math:`S_{c}` is reduced by an
increment :math:`i`, to reduce the spatial extent of the canopy, thus.

.. math::

   S_{c,patch,Cl,t+1} = \left\{ \begin{array}{ll}
   S_{c,patch,Cl,t} + i& \textrm{for $A_{f,Cl} < A_{t}$}\\
   &\\
   S_{c,patch,Cl,t} - i& \textrm{for $A_{f,Cl} > A_{t}$}\\
   \end{array} \right.

The values of :math:`S_{c}` are bounded to upper and lower limits. The
lower limit corresponds to the observed canopy spread parameter for
canopy trees :math:`S_{c,min}` and the upper limit corresponds to the
largest canopy extent :math:`S_{c,max}`

.. math::

   S_{c,patch,Cl} = \left\{ \begin{array}{ll}
   S_{c,min}& \textrm{for } S_{c,patch,Cl}< S_{c,\rm{min}}\\
   &\\
   S_{c,max}& \textrm{for } S_{c,patch,Cl} > S_{c,\rm{max}}\\
   \end{array} \right.

This iterative scheme requires two additional parameters (:math:`i` and
:math:`A_{t}`). :math:`i` affects the speed with which canopy spread
(and hence leaf are index) increase as canopy closure is neared.
However, the model is relatively insensitive to the choice of either
:math:`i` or :math:`A_{t}`.

Definition of Leaf and Stem Area Profile
----------------------------------------

Within each patch, the model defines and tracks cohorts of multiple
plant functional types that exist either in the canopy or understorey.
Light on the top leaf surface of each cohort in the canopy is the same,
and the rate of decay through the canopy is also the same for each PFT.
Therefore, we accumulate all the cohorts of a given PFT together for the
sake of the radiation and photosynthesis calculations (to avoid separate
calculations for every cohort).

Therefore, the leaf area index for each patch is defined as a
three-dimensional array :math:`\mathit{lai}_{Cl,ft,z}` where :math:`C_l`
is the canopy layer, :math:`ft` is the functional type and :math:`z` is
the leaf layer within each canopy. This three-dimensional structure is
the basis of the radiation and photosynthetic models. In addition to a
leaf area profile matrix, we also define, for each patch, the area which
is covered by leaves at each layer as :math:`\mathit{carea}_{Cl,ft,z}`.

Each plant cohort is already defined as a member of a single canopy
layer and functional type. This means that to generate the
:math:`x_{Cl,ft,z}` matrix, it only remains to divide the leaf area of
each cohort into leaf layers. First, we determine how many leaf layers
are occupied by a single cohort, by calculating the ‘tree LAI’ as the
total leaf area of each cohort divided by its crown area (both in
m\ :math:`^{2}`)

.. math:: \mathit{tree}_{lai,coh} = \frac{b_{leaf,coh}\cdot\mathrm{sla}_{ft}}{A_{crown,coh}}

where :math:`\mathrm{sla}_{ft}` is the specific leaf area in
m\ :math:`^{2}` KgC\ :math:`^{-1}` and :math:`b_{leaf}` is in kGC per
plant.

Stem area index (SAI) is ratio of the total area of all woody stems on a
plant to the area of ground covered by the plant. During winter in
deciduous areas, the extra absorption by woody stems can have a
significant impact on the surface energy budget. However, in previous
`big leaf` versions of the CLM, computing the circumstances under which
stem area was visible in the absence of leaves was difficult and the
algorithm was largely heuristic as a result. Given the multi-layer
canopy introduced for FATES, we can determine the leaves in the higher
canopy layers will likely shade stem area in the lower layers when
leaves are on, and therefore stem area index can be calculated as a
function of woody biomass directly.

Literature on stem area index is particularly poor, as it’s estimation
is complex and not particularly amenable to the use of, for example,
assumptions of random distribution in space that are typically used to
calculate leaf area from light interception.
:ref:`Kucharik et al. 1998<kucharik1998>` estimated that SAI visible from an
LAI2000 sensor was around 0.5 m^2 m^-2. Low et al. 2001
estimate that the wood area index for Ponderosa Pine forest is
0.27-0.33. The existing CLM(CN) algorithm sets the minimum SAI at 0.25
to match MODIS observations, but then allows SAI to rise as a function
of the LAI lost, meaning than in some places, predicted SAI can reach
value of 8 or more. Clearly, greater scientific input on this quantity
is badly needed. Here we determine that SAI is a linear function of
woody biomass, to at very least provide a mechanistic link between the
existence of wood and radiation absorbed by it. The non-linearity
between how much woody area exists and how much radiation is absorbed is
provided by the radiation absorption algorithm. Specifically, the SAI of
an individual cohort (:math:`\mathrm{tree}_{sai,coh}`, m\ :math:`^{2}`
m\ :math:`^{-2}`) is calculated as follows,

.. math:: \mathrm{tree}_{sai,coh} = k_{sai}\cdot b_{struc,coh} ,

where :math:`k_{sai}` is the coefficient linking structural biomass to
SAI. The number of occupied leaf layers for cohort :math:`coh`
(:math:`n_{z,coh}`) is then equal to the rounded up integer value of the
tree SAI (:math:`{tree}_{sai,coh}`) and LAI (:math:`{tree}_{lai,coh}`)
divided by the layer thickness (i.e., the resolution of the canopy layer
model, in units of vegetation index (:math:`lai`\ +\ :math:`sai`) with a
default value of 1.0, :math:`\delta _{vai}` ),

.. math:: n_{z,coh} = {\frac{\mathrm{tree}_{lai,coh}+\mathrm{tree}_{sai,coh}}{\delta_{vai}}}.

The fraction of each layer that is leaf (as opposed to stem) can then be
calculated as

.. math:: f_{leaf,coh} = \frac{\mathrm{tree}_{lai,coh}}{\mathrm{tree}_{sai,coh}+\mathrm{tree}_{lai,coh}}.

Finally, the leaf area in each leaf layer pertaining to this cohort is
thus

.. math::

   \mathit{lai}_{z,coh}  = \left\{ \begin{array}{ll}
    \delta_{vai} \cdot f_{leaf,coh} \frac{A_{canopy,coh}}{A_{canopy,patch}}& \textrm{for $i=1,..., i=n_{z,coh}-1$}\\
   &\\
    \delta_{vai} \cdot f_{leaf,coh} \frac{A_{canopy,coh}}{A_{canopy,patch}}\cdot r_{vai}& \textrm{for $i=n_{z,coh}$}\\
   \end{array} \right.

and the stem area index is

.. math::

   \mathit{sai}_{z,coh}  = \left\{ \begin{array}{ll}
    \delta_{vai} \cdot (1-f_{leaf,coh})\frac{A_{canopy,coh}}{A_{canopy,patch}}& \textrm{for $i=1,..., i=n_{z,coh}-1$}\\
   &\\
    \delta_{vai} \cdot (1-f_{leaf,coh}) \frac{A_{canopy,coh}}{A_{canopy,patch}}\cdot r_{vai}& \textrm{for $i=n_{z,coh}$}\\
   \end{array} \right.

where :math:`r_{vai}` is the remainder of the canopy that is below the
last full leaf layer

.. math:: r_{vai} =(\mathrm{tree}_{lai,coh} + \mathrm{tree}_{sai,coh}) - (\delta _{vai} \cdot (n_{z,coh} -1)).

:math:`A_{canopy,patch}` is the total canopy area occupied by plants in
a given patch (m:math:`^{2}`) and is calculated as follows,

.. math:: A_{canopy,patch} = \textrm{min}\left( \sum_{coh=1}^{coh = ncoh}A_{canopy,coh}, A_{patch}  \right).

The canopy is conceived as a cylinder, although this assumption could be
altered given sufficient evidence that canopy shape was an important
determinant of competitive outcomes, and the area of ground covered by
each leaf layer is the same through the cohort canopy. With the
calculated SAI and LAI, we are able to calculate the complete canopy
profile. Specifically, the relative canopy area for the cohort
:math:`{coh}` is calculated as

.. math:: \mathit{area}_{1:nz,coh}  =  \frac{A_{crown,coh}}{A_{canopy,patch}}.

The total occupied canopy area for each canopy layer (:math:`Cl`), plant
functional type (:math:`ft`) and leaf layer (:math:`z`) bin is thus

.. math::
  
   \mathit{c}_{area,Cl,ft,z} = \sum_{coh=1}^{coh=ncoh} area_{1:nz,coh} 

where :math:`ft_{coh}=ft`  and  :math:`Cl_{coh} = Cl.`

All of these quantities are summed across cohorts to give the complete
leaf and stem area profiles,

.. math::

   \mathit{lai} _{Cl,ft,z} = \sum_{coh=1}^{coh=ncoh} \mathit{lai}_{z,coh}  

.. math::

   \mathit{sai}_{Cl,ft,z} = \sum_{coh=1}^{coh=ncoh} \mathit{sai}_{z,coh}  
   

Burial of leaf area by snow
---------------------------

The calculations above all pertain to the total leaf and stem area
indices which charecterize the vegetation structure. In addition, the
model must know when the vegetation is covered by snow, and by how much,
so that the albedo and energy balance calculations can be adjusted
accordingly. Therefore, we calculated a ‘total’ and ‘exposed’
:math:`lai` and :math:`sai` profile using a representation of the bottom
and top canopy heights, and the depth of the average snow pack. For each
leaf layer :math:`z` of each cohort, we calculate an ‘exposed fraction
:math:`f_{exp,z}` via consideration of the top and bottom heights of
that layer :math:`h_{top,z}` and :math:`h_{bot,z}` (m),

.. math::

   \begin{array}{ll}
   h_{top,z} = h_{coh} - h_{coh}\cdot f_{crown,ft}\cdot\frac{z}{n_{z,coh}}& \\
   &\\
   h_{bot,z} = h_{coh} - h_{coh}\cdot f_{crown,ft}\cdot\frac{z+1}{n_{z,coh}}&\\
   \end{array}

where :math:`f_{crown,ft}` is the plant functional type (:math:`ft`)
specific fraction of the cohort height that is occupied by the crown.
Specifically, the ‘exposed fraction :math:`f_{exp,z}` is calculated as
follows,

.. math::

   f_{exp,z}\left\{ \begin{array}{ll}
   = 1.0 &  h_{bot,z}> d_{snow}\\
   &\\
   = \frac{d_{snow} -h_{bot,z}}{h_{top,z}-h_{bot,z}}  & h_{top,z}> d_{snow}, h_{bot,z}< d_{snow}\\
   &\\
   = 0.0 & h_{top,z}< d_{snow}\\
   \end{array} \right.

The resulting exposed (:math:`elai, esai`) and total
(:math:`tlai, tsai`) leaf and stem area indicies are calculated as

.. math::

   \begin{array}{ll}
   \mathit{elai} _{Cl,ft,z} &= \mathit{lai} _{Cl,ft,z} \cdot f_{exp,z}\\
   \mathit{esai} _{Cl,ft,z} &= \mathit{sai} _{Cl,ft,z} \cdot f_{exp,z}\\
   \mathit{tlai} _{Cl,ft,z} &= \mathit{lai} _{Cl,ft,z}\\
   \mathit{tsai} _{Cl,ft,z} &= \mathit{sai} _{Cl,ft,z} \
   \end{array} ,

and are used in the radiation interception and photosynthesis algorithms
described later.

+-------------+-------------+-------------+-------------+-------------+
| Parameter   | Parameter   | Units       | Notes       | Indexed by  |
| Symbol      | Name        |             |             |             |
+=============+=============+=============+=============+=============+
| :math:`     | Thickness   | m\ :math:`^ |             |             |
| \delta_     | of single   | {-2}`\ m\ : |             |             |
| {vai}`      | canopy      | math:`^{-2}`|             |             |
|             | layer       |             |             |             |
+-------------+-------------+-------------+-------------+-------------+
| :math:`C_e` | Competitive | none        |             |             |
|             | Exclusion   |             |             |             |
|             | Parameter   |             |             |             |
+-------------+-------------+-------------+-------------+-------------+
| :math:`c_{p | Minimum     | m\ :math:`^ |             |             |
| ,min}`      | canopy      | {2}`        |             |             |
|             | spread      | cm\ :math:` |             |             |
|             |             | ^{-1}`      |             |             |
+-------------+-------------+-------------+-------------+-------------+
| :math:`c_{p | Competitive | m\ :math:`^ |             |             |
| ,max}`      | Exclusion   | {2}`        |             |             |
|             | Parameter   | cm\ :math:` |             |             |
|             |             | ^{-1}`      |             |             |
+-------------+-------------+-------------+-------------+-------------+
| :math:`i`   | Incremental | m\ :math:`^ |             |             |
|             | change in   | {2}`        |             |             |
|             | :math:`c_p` | cm\ :math:` |             |             |
|             |             | ^{-1}`      |             |             |
|             |             | y\ :math:`^ |             |             |
|             |             | {-1}`       |             |             |
+-------------+-------------+-------------+-------------+-------------+
| :math:`A_t` | Threshold   | none        |             |             |
|             | canopy      |             |             |             |
|             | closure     |             |             |             |
+-------------+-------------+-------------+-------------+-------------+
| :math:`f_{c | Crown       | none        |             | :math:`ft`  |
| rown,ft}`   | fraction    |             |             |             |
+-------------+-------------+-------------+-------------+-------------+
| :math:`k_{s | Stem area   | m^2 KgC^-1  |             |             |
| ai}`        | per unit    |             |             |             |
|             | woody       |             |             |             |
|             | biomass     |             |             |             |
+-------------+-------------+-------------+-------------+-------------+


Radiation Transfer
^^^^^^^^^^^^^^^^^^^

Fundamental Radiation Transfer Theory
-------------------------------------

The first interaction of the land surface with the properties of
vegetation concerns the partitioning of energy into that which is
absorbed by vegetation, reflected back into the atmosphere, and absorbed
by the ground surface. Older versions of the CLM have utilized a
"two-stream" approximation
:ref:`Sellers 1985<sellers1985>`, :ref:`Sellers et al. 1986<sellers1996>` that provided an
empirical solution for the radiation partitioning of a multi-layer
canopy for two streams, of diffuse and direct light. However,
implementation of the Ecosystem Demography model requires a) the
adoption of an explicit multiple layer canopy b) the implementation of a
multiple plant type canopy and c) the distinction of canopy and
under-storey layers, in-between which the radiation streams are fully
mixed. The radiation mixing between canopy layers is necessary as the
position of different plants in the under-storey is not defined
spatially or relative to the canopy trees above. In this new scheme, we
thus implemented a one-dimensional scheme that traces the absorption,
transmittance and reflectance of each canopy layer and the soil,
iterating the upwards and downwards passes of radiation through the
canopy until a pre-defined accuracy tolerance is reached. This approach
is based on the work of :ref:`Norman 1979<norman1979>`.

Here we describe the basic theory of the radiation transfer model for
the case of a single homogenous canopy, and in the next section we
discuss how this is applied to the multi layer multi PFT canopy in the
FATES implementation. The code considers the fractions of a single
unit of incoming direct and a single unit of incoming diffuse light,
that are absorbed at each layer of the canopy for a given solar angle
(:math:`\alpha_{s}`, radians). Direct radiation is extinguished through
the canopy according to the coefficient :math:`k_{dir}` that is
calculated from the incoming solar angle and the dimensionless leaf
angle distribution parameter (:math:`\chi`) as

.. math:: k_{dir} = g_{dir} / \sin(\alpha_s)\\

where

.. math:: g_{dir} = \phi_1 + \phi_2 \cdot \sin(\alpha_s)\\

and

.. math::

   \begin{array} {l}
   \phi_1 = 0.5 - 0.633\chi_{l} - 0.33\chi_l ^2\\
   \phi_2 =0.877 (1 - 2\phi_1)\\

   \end{array}

The leaf angle distribution is a descriptor of how leaf surfaces are
arranged in space. Values approaching 1.0 indicate that (on average) the
majority of leaves are horizontally arranged with respect to the ground.
Values approaching -1.0 indicate that leaves are mostly vertically
arranged, and a value of 0.0 denotes a canopy where leaf angle is random
(a ‘spherical’ distribution).

According to Beer’s Law, the fraction of light that is transferred
through a single layer of vegetation (leaves or stems) of thickness
:math:`\delta_{vai}`, without being intercepted by any surface, is

.. math:: \mathit{tr}_{dir} = e^{-k_{dir}  \delta_{vai}}

and the incident direct radiation transmitted to each layer of the
canopy (:math:`dir_{tr,z}`) is thus calculated from the cumulative leaf
area ( :math:`L_{above}` ) shading each layer (:math:`z`):

.. math:: \mathit{dir}_{tr,z} = e^{-k_{dir}  L_{above,z}}

The fraction of the leaves :math:`f_{sun}` that are exposed to direct
light is also calculated from the decay coefficient :math:`k_{dir}`.

.. math::

   \begin{array}{l}
   f_{sun,z} = e^{-k_{dir}  L_{above,z}}\\
    \rm{and} 
   \\ f_{shade,z} = 1-f_{sun,z}
   \end{array}

where :math:`f_{shade,z}` is the fraction of leaves that are shaded
from direct radiation and only receive diffuse light.

Diffuse radiation, by definition, enters the canopy from a spectrum of
potential incident directions, therefore the un-intercepted transfer
(:math:`tr_{dif}`) through a leaf layer of thickness :math:`\delta_l` is
calculated as the mean of the transfer rate from each of 9 different
incident light directions (:math:`\alpha_{s}`) between 0 and 180 degrees
to the horizontal.

.. math:: \mathit{tr}_{dif} = \frac{1}{9} \sum\limits_{\alpha_s=5\pi/180}^{\alpha_s=85\pi/180} e^{-k_{dir,l} \delta_{vai}} \\ \\

.. math:: tr_{dif}= \frac{1}{9} \pi \sum_{\alpha s=0}^{ \pi / 2}  \frac{e^{-gdir} \alpha_s}{\delta_{vai} \cdot \rm{sin}(\alpha_s) \rm{sin}(\alpha_s) \rm{cos}(\alpha_s)}

The fraction (1-:math:`tr_{dif}`) of the diffuse radiation is
intercepted by leaves as it passes through each leaf layer. Of this,
some fraction is reflected by the leaf surfaces and some is transmitted
through. The fractions of diffuse radiation reflected from
(:math:`\mathit{refl}_{dif}`) and transmitted though
(:math:`\mathit{tran}_{dif}`) each layer of leaves are thus,
respectively

.. math::

   \begin{array}{l}
   \mathit{refl_{dif}} = (1 - tr_{dif})  \rho_{l,ft}\\
   \mathit{tran}_{dif} = (1 - tr_{dif})  \tau_{l,ft} + tr_{dif}
   \end{array}

where :math:`\rho_{l,ft}` and :math:`\tau_{l,ft}` are the fractions of
incident light reflected and transmitted by individual leaf surfaces.

Once we know the fractions of light that are transmitted and reflected
by each leaf layer, we begin the process of distributing light through
the canopy. Starting with the first leaf layer (:math:`z`\ =1), where
the incident downwards diffuse radiation (:math:`\mathit{dif}_{down}`)
is 1.0, we work downwards for :math:`n_z` layers, calculating the
radiation in the next layer down (:math:`z+1`) as:

.. math:: \mathit{dif}_{down,z+1} = \frac{\mathit{dif}_{down,z} \mathit{tran}_{dif} }    {1 - \mathit{r}_{z+1}  \mathit{refl}_{dif}}

Here, :math:`\mathit{dif}_{down,z} \mathit{tran}_{dif}` calculates the
fraction of incoming energy transmitted downwards onto layer
:math:`z+1`. This flux is then increased by the additional radiation
:math:`r_z` that is reflected upwards from further down in the canopy to
layer :math:`z`, and then is reflected back downwards according to the
reflected fraction :math:`\mathit{refl_{dif}}`. The more radiation in
:math:`\mathit{r}_{z+1}  \mathit{refl}_{dif}`, the smaller the
denominator and the larger the downwards flux. :math:`r` is also
calculated sequentially, starting this time at the soil surface layer
(where :math:`z = n_z+1`)

.. math:: r_{nz+1} = alb_s

where :math:`alb_s` is the soil albedo characteristic. The upwards
reflected fraction :math:`r_z` for each leaf layer, moving upwards, is
then :ref:`Norman 1979<norman1979>`

.. math:: r_z  = \frac{r_{z+1}  \times \mathit{tran}_{dif}  ^{2} }{ (1 - r_{z+1}  \mathit{refl_{dif}}) + \mathit{refl_{dif}}}.

The corresponding upwards diffuse radiation flux is therefore the
fraction of downwards radiation that is incident on a particular layer,
multiplied by the fraction that is reflected from all the lower layers:

.. math:: \mathit{dif}_{up,z} = r_z \mathit{dif}_{down,z+1}

Now we have initial conditions for the upwards and downwards diffuse
fluxes, these must be modified to account for the fact that, on
interception with leaves, direct radiation is transformed into diffuse
radiation. In addition, the initial solutions to the upwards and
downwards radiation only allow a single ‘bounce’ of radiation through
the canopy, so some radiation which might be intercepted by leaves
higher up is potentially lost. Therefore, the solution to this model is
iterative. The iterative solution has upwards and a downwards components
that calculate the upwards and downwards fluxes of total radiation at
each leaf layer (:math:`rad_{dn, z}` and :math:`rad_{up, z}`) . The
downwards component begins at the top canopy layer (:math:`z=1`). Here
we define the incoming solar diffuse and direct radiation
(:math:`\it{solar}_{dir}` and :math:`\it{solar}_{dir}` respectively).

.. math::

   \begin{array}{l}
    \mathit{dif}_{dn,1} =  \it{solar}_{dif} \\
   \mathit{rad}_{dn, z+1} = \mathit{dif}_{dn,z} \cdot  \mathit{tran}_{dif}  +\mathit{dif}_{up,z+1}   \cdot  \mathit{refl}_{dif}   + \mathit{solar}_{dir}  \cdot  dir_{tr,z}  (1- tr_{dir})  \tau_l.
   \end{array}

The first term of the right-hand side deals with the diffuse radiation
transmitted downwards, the second with the diffuse radiation travelling
upwards, and the third with the direct radiation incoming at each layer
(:math:`dir_{tr,z}`) that is intercepted by leaves
(:math:`1-  tr_{dir}`) and then transmitted through through the leaf
matrix as diffuse radiation (:math:`\tau_l`). At the bottom of the
canopy, the light reflected off the soil surface is calculated as

.. math:: rad _{up, nz} =  \rm{\it{dif}}_{down,z}  \cdot  salb_{dif} +\it{solar}_{dir} \cdot dir_{tr,z} salb_{dir}.

The upwards propagation of the reflected radiation is then

.. math:: rad_{up, z} = \mathit{dif}_{up,z+1} \cdot  \mathit{tran}_{dif}  +\mathit{dif}_{dn,z}   \cdot  \mathit{refl}_{dif}   + \it{solar}_{dir}  \cdot  dir_{tr,z}  (1- tr_{dir})  \rho_l.

Here the first two terms deal with the diffuse downwards and upwards
fluxes, as before, and the third deals direct beam light that is
intercepted by leaves and reflected upwards. These upwards and downwards
fluxes are computed for multiple iterations, and at each iteration,
:math:`rad_{up, z}` and :math:`rad_{down, z}` are compared to their
values in the previous iteration. The iteration scheme stops once the
differences between iterations for all layers is below a predefined
tolerance factor, (set here at :math:`10^{-4}`). Subsequently, the
fractions of absorbed direct (:math:`abs_{dir,z}`) and diffuse
(:math:`abs_{dif,z}`) radiation for each leaf layer then

.. math:: abs_{dir,z} = \it{solar}_{dir}   \cdot dir_{tr,z} \cdot (1- tr_{dir}) \cdot (1 - \rho_l-\tau_l)

.. math:: abs_{dif,z} = (\mathit{dif}_{dn,z} +  \mathit{dif}_{up,z+1} ) \cdot (1 - tr_{dif}) \cdot (1 - \rho_l-\tau_l).

and, the radiation energy absorbed by the soil for the diffuse and
direct streams is is calculated as

.. math:: \it{abs}_{soil} = \mathit{dif}_{down,nz+1} \cdot (1 -  salb_{dif}) +\it{solar}_{dir}   \cdot dir_{tr,nz+1} \cdot (1-  salb_{dir}).

Canopy level albedo is denoted as the upwards flux from the top leaf
layer

.. math:: \it{alb}_{canopy}=  \frac{\mathit{dif}_{up,z+1}  }{  \it{solar}_{dir} + \it{solar}_{dif}}

and the division of absorbed energy into sunlit and shaded leaf
fractions, (required by the photosynthesis calculations), is

.. math:: abs_{sha,z} = abs_{dif,z} \cdot f_{sha}

.. math:: abs_{sun,z} =  abs_{dif,z} \cdot f_{sun}+ abs_{dir,z}

Resolution of radiation transfer theory within the FATES canopy structure
-------------------------------------------------------------------------

The radiation transfer theory above, was described with reference to a
single canopy of one plant functional type, for the sake of clarity of
explanation. The FATES model, however, calculates radiative and
photosynthetic fluxes for a more complex hierarchical structure within
each patch/time-since-disturbance class, as described in the leaf area
profile section. Firstly, we denote two or more canopy layers (denoted
:math:`C_l`). The concept of a ‘canopy layer’ refers to the idea that
plants are organized into discrete over and under-stories, as predicted
by the Perfect Plasticity Approximation
(:ref:`Purves et al. 2008<purves2008>`, :ref:`Fisher et al. 2010<Fisheretal2010>`). Within each canopy layer
there potentially exist multiple cohorts of different plant functional
types and heights. Within each canopy layer, :math:`C_l`, and functional
type, :math:`ft`, the model resolves numerous leaf layers :math:`z`,
and, for some processes, notably photosynthesis, each leaf layer is
split into a fraction of sun and shade leaves, :math:`f_{sun}` and
:math:`f_{sha}`, respectively.

The radiation scheme described in Section is solved explicitly for this
structure, for both the visible and near-infrared wavebands, according
to the following assumptions.

-  A *canopy layer* (:math:`C_{L}`) refers to either the over or understorey

-  A *leaf layer* (:math:`z`) refers to the discretization of the LAI
   within the canopy of a given plant functional type.

-  All PFTs in the same canopy layer have the same solar radiation
   incident on the top layer of the canopy

-  Light is transmitted through the canopy of each plant functional type independently

-  Between canopy layers, the light streams from different plant
   functional types are mixed, such that the (undefined) spatial
   location of plants in lower canopy layers does not impact the amount
   of light received.

-  Where understorey layers fill less area than the overstorey layers,
   radiation is directly transferred to the soil surface.

-  All these calculations pertain to a single patch, so we omit the
   `patch` subscript for simplicity in the following discussion.

Within this framework, the majority of the terms in the radiative
transfer scheme are calculated with indices of :math:`C_L`,
:math:`\it{ft}` and :math:`z`. In the following text, we revisit the
simplified version of the radiation model described above, and explain
how it is modified to account for the more complex canopy structure used
by FATES.

Firstly, the light penetration functions, :math:`k_{dir}` and
:math:`g_{dir}` are described as functions of :math:`\it{ft}`, because
the leaf angle distribution, :math:`\chi_l`, is a pft-specific
parameter. Thus, the diffuse irradiance transfer rate, :math:`tr_{dif}`
is also :math:`\it{ft}` specific because :math:`g_{dir}`, on which it
depends, is a function of :math:`\chi_l`.

The amount of direct light reaching each leaf layer is a function of the
leaves existing above the layer in question. If a leaf layer ‘:math:`z`’
is in the top canopy layer (the over-storey), it is only shaded by
leaves of the same PFT so :math:`k_{dir}` is unchanged from equation. If
there is more than one canopy layer (:math:`C_{l,max}>1`), then the
amount of direct light reaching the top leaf surfaces of the
second/lower layer is the weighted average of the light attenuated by
all the parallel tree canopies in the canopy layer above, thus.

.. math:: dir_{tr,Cl,:,1} =\sum_{ft=1}^{npft}{(dir_{tr,Cl,ft,z_{max}} \cdot c_{area,Cl-1,ft,z_{max}})}

where :math:`\it{pft}_{wt}` is the areal fraction of each canopy layer
occupied by each functional type and :math:`z_{max}` is the index of the
bottom canopy layer of each pft in each canopy layer (the subscripts
:math:`C_l` and :math:`ft` are implied but omitted from all
:math:`z_{max}` references to avoid additional complications)

Similarly, the sunlit fraction for a leaf layer ‘:math:`z`’ in the
second canopy layer (where :math:`C_l > 1`) is

.. math:: f_{sun,Cl,ft,z} = W_{sun,Cl} \cdot e^{k_{dir,ft,laic,z}}

where :math:`W_{sun,Cl}` is the weighted average sunlit fraction in the
bottom layer of a given canopy layer.

.. math:: W_{sun,Cl} = \sum_{ft=1}^{npft}{(f_{sun,Cl-1,ft,zmax} \cdot  c_{area,Cl-1,ft,zmax})}

Following through the sequence of equations for the simple single pft
and canopy layer approach above, the :math:`\mathit{refl}_{dif}` and
:math:`\mathit{tran}_{dif}` fluxes are also indexed by :math:`C_l`,
:math:`\it{ft}`, and :math:`z`. The diffuse radiation reflectance ratio
:math:`r_z` is also calculated in a manner that homogenizes fluxes
between canopy layers. For the canopy layer nearest the soil
(:math:`C_l` = :math:`C_{l,max}`). For the top canopy layer
(:math:`C_l`\ =1), a weighted average reflectance from the lower layers
is used as the baseline, in lieu of the soil albedo. Thus:

.. math:: r_{z,Cl,:,1} =  \sum_{ft=1}^{npft}{(r_{z,Cl-1,ft,1}   \it{pft}_{wt,Cl-1,ft,1})}

For the iterative flux resolution, the upwards and downwards fluxes are
also averaged between canopy layers, thus where :math:`C_l>1`

.. math:: rad_{dn, Cl,ft,1} = \sum_{ft=1}^{npft}{(rad_{dn, Cl-1,ft,zmax} \cdot  \it{pft}_{wt,Cl-1,ft,zmax})}

and where :math:`C_l` =1, and :math:`C_{l,max}>1`

.. math:: rad_{up,Cl,ft,zmax} = \sum_{ft=1}^{npft}{(rad_{up, Cl+1,ft,1} \cdot  \it{pft}_{wt,Cl+1,ft,1})}

The remaining terms in the radiation calculations are all also indexed
by :math:`C_l`, :math:`ft` and :math:`z` so that the fraction of
absorbed radiation outputs are termed :math:`abs_{dir,Cl,ft,z}` and
:math:`abs_{dif,Cl,ft,z}`. The sunlit and shaded absorption rates are
therefore

.. math:: abs_{sha,Cl,ft,z} = abs_{dif,Cl,ft,z}\cdot f_{sha,Cl,ft,z}

and

.. math:: abs_{sun,Cl,ft,z} =  abs_{dif,Cl,ft,z} \cdot f_{sun,Cl,ft,z}+ abs_{dir,Cl,ft,z}

The albedo of the mixed pft canopy is calculated as the weighted average
of the upwards radiation from the top leaf layer of each pft where
:math:`C_l`\ =1:

.. math:: \it{alb}_{canopy}=  \sum_{ft=1}^{npft}{\frac{\mathit{dif}_{up,1,ft,1}    \it{pft}_{wt,1,ft,1}} {\it{solar}_{dir} + \it{solar}_{dif}}}

The radiation absorbed by the soil after passing through through
under-storey vegetation is:

.. math:: \it{abs}_{soil}=  \sum_{ft=1}^{npft}{ \it{pft}_{wt,1,ft,1}( \mathit{dif}_{down,nz+1} (1 -  salb_{dif}) +\it{solar}_{dir}   dir_{tr,nz+1}  (1-  salb_{dir}))}

to which is added the diffuse flux coming directly from the upper
canopy and hitting no understorey vegetation.

.. math:: \it{abs}_{soil}=  \it{abs}_{soil}+dif_{dn,2,1}  (1-  \sum_{ft=1}^{npft}{\it{pft}_{wt,1,ft,1}})  (1 -  salb_{dif})

and the direct flux coming directly from the upper canopy and hitting
no understorey vegetation.

.. math:: \it{abs}_{soil}=  \it{abs}_{soil}+\it{solar}_{dir} dir_{tr,2,1}(1-  \sum_{ft=1}^{npft}{\it{pft}_{wt,1,ft,1}})  (1 -  salb_{dir})

These changes to the radiation code are designed to be structurally
flexible, and the scheme may be collapsed down to only include on canopy
layer, functional type and pft for testing if necessary.

.. raw:: latex

   \captionof{table}{Parameters needed for radiation transfer model. }

+-----------------+-----------------+-----------------+-----------------+
| Parameter       | Parameter Name  | Units           | indexed by      |
| Symbol          |                 |                 |                 |
+=================+=================+=================+=================+
| :math:`\chi`    | Leaf angle      | none            | *ft*            |
|                 | distribution    |                 |                 |
|                 | parameter       |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`\rho_l`  | Fraction of     | none            | *ft*            |
|                 | light reflected |                 |                 |
|                 | by leaf surface |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`\tau_l`  | Fraction of     | none            | *ft*            |
|                 | light           |                 |                 |
|                 | transmitted by  |                 |                 |
|                 | leaf surface    |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`alb_s`   | Fraction of     | none            | direct vs       |
|                 | light reflected |                 | diffuse         |
|                 | by soil         |                 |                 |
+-----------------+-----------------+-----------------+-----------------+

.. raw:: latex

   \bigskip 

Photosynthesis
^^^^^^^^^^^^^^^^^^^^

Fundamental photosynthetic physiology theory
--------------------------------------------

In this section we describe the physiological basis of the
photosynthesis model before describing its application to the FATES
canopy structure. This description in this section is largely repeated
from the Oleson et al. CLM4.5 technical note but included here for
comparison with its implementation in FATES. Photosynthesis in C3
plants is based on the model of :ref:`Farquhar 1980<Farquharetal1980>` as
modified by :ref:`Collatz et al. (1991)<Collatzetal1991>`. Photosynthetic assimilation
in C4 plants is based on the model of :ref:`Collatz et al. (1991)<Collatzetal1991>`.
In both models, leaf photosynthesis, :math:`\textrm{gpp}`
(:math:`\mu`\ mol CO\ :math:`_2` m\ :math:`^{-2}` s\ :math:`^{-1}`) is
calculated as the minimum of three potentially limiting fluxes,
described below:

.. math:: \textrm{gpp} = \rm{min}(w_{j}, w_{c},w_{p}).

The RuBP carboxylase (Rubisco) limited rate of carboxylation
:math:`w_{c}` (:math:`\mu`\ mol CO\ :math:`_{2}` m\ :math:`^{-2}`
s\ :math:`^{-1}`) is determined as

.. math::

   w_{c}=  \left\{ \begin{array}{ll}
   \frac{V_{c,max}(c_{i} - \Gamma_*)}{ci+K_{c}(1+o_{i}/K_{o})} & \textrm{for $C_{3}$ plants}\\
   &\\
   V_{c,max}& \textrm{for $C_{4}$ plants}\\
   \end{array} \right.
   c_{i}-\Gamma_*\ge 0

where :math:`c_{i}` is the internal leaf CO\ :math:`_{2}` partial
pressure (Pa) and :math:`o_i (0.209P_{atm}`) is the O\ :math:`_{2}`
partial pressure (Pa). :math:`K_{c}` and :math:`K_{o}` are the
Michaelis-Menten constants (Pa) for CO\ :math:`_{2}` and
O\ :math:`_{2}`. These vary with vegetation temperature :math:`T_v`
(:math:`^{o}`\ C) according to an Arrhenious function described in
:ref:`Oleson et al. 2013<olesonetal2013>`. :math:`V_{c,max}` is the leaf layer
photosynthetic capacity (:math:`\mu` mol CO\ :math:`_2` m\ :math:`^{-2}`
s\ :math:`^{-1}`).

The maximum rate of carboxylation allowed by the capacity to regenerate
RuBP (i.e., the light-limited rate) :math:`w_{j}` (:math:`\mu`\ mol
CO\ :math:`_2` m\ :math:`^{-2}` s\ :math:`^{-1}`) is

.. math::

   w_j=  \left\{ \begin{array}{ll}
   \frac{J(c_i - \Gamma_*)}{4ci+8\Gamma_*} & \textrm{for C$_3$ plants}\\
   &\\
   4.6\phi\alpha & \textrm{for C$_4$ plants}\\
   \end{array} \right.
   c_i-\Gamma_*\ge 0

To find :math:`J`, the electron transport rate (:math:`\mu` mol
CO\ :math:`_2` m\ :math:`^{-2}` s\ :math:`^{-1}`), we solve the
following quadratic term and take its smaller root,

.. math:: \Theta_{psII}J^{2}-(I_{psII} +J_{max})J+I_{psII}J_{max} =0

where :math:`J_{max}` is the maximum potential rate of electron
transport (:math:`\mu`\ mol m\ :math:`_{-2}` s\ :math:`^{-1}`),
:math:`I_{PSII}` is the is the light utilized in electron transport by
photosystem II (:math:`\mu`\ mol m\ :math:`_{-2}` s\ :math:`^{-1}`) and
:math:`\Theta_{PSII}` is is curvature parameter. :math:`I_{PSII}` is
determined as

.. math:: I_{PSII} =0.5 \Phi_{PSII}(4.6\phi)

where :math:`\phi` is the absorbed photosynthetically active radiation
(Wm:math:`^{-2}`) for either sunlit or shaded leaves (:math:`abs_{sun}`
and :math:`abs_{sha}`). :math:`\phi` is converted to photosynthetic
photon flux assuming 4.6 :math:`\mu`\ mol photons per joule. Parameter
values are :math:`\Phi_{PSII}` = 0.7 for C3 and :math:`\Phi_{PSII}` =
0.85 for C4 plants.

The export limited rate of carboxylation for C3 plants and the PEP
carboxylase limited rate of carboxylation for C4 plants :math:`w_e`
(also in :math:`\mu`\ mol CO\ :math:`_2` m\ :math:`^{-2}`
s\ :math:`^{-1}`) is

.. math::

   w_e=  \left\{ \begin{array}{ll}
   3 T_{p,0} & \textrm{for $C_3$ plants}\\
   &\\
   k_{p} \frac{c_i}{P_{atm}}& \textrm{for $C_4$ plants}.\\
   \end{array} \right.

:math:`T_{p}` is the triose-phosphate limited rate of photosynthesis,
which is equal to :math:`0.167 V_{c,max0}`. :math:`k_{p}` is the initial
slope of C4 CO\ :math:`_{2}` response curve. The Michaelis-Menten
constants :math:`K_{c}` and :math:`K_{o}` are modeled as follows,

.. math:: K_{c} = K_{c,25}(a_{kc})^{\frac{T_v-25}{10}},

.. math:: K_{o} = K_{o,25}(a_{ko})^{\frac{T_v-25}{10}},

where :math:`K_{c,25}` = 30.0 and :math:`K_{o,25}` = 30000.0 are values
(Pa) at 25 :math:`^{o}`\ C, and :math:`a_{kc}` = 2.1 and :math:`a_{ko}`
=1.2 are the relative changes in :math:`K_{c,25}` and :math:`K_{o,25}`
respectively, for a 10\ :math:`^{o}`\ C change in temperature. The
CO\ :math:`_{2}` compensation point :math:`\Gamma_{*}` (Pa) is

.. math:: \Gamma_* = \frac{1}{2} \frac{K_c}{K_o}0.21o_i

where the term 0.21 represents the ratio of maximum rates of oxygenation
to carboxylation, which is virtually constant with temperature
:ref:`Farquhar, 1980<Farquharetal1980>`.

Resolution of the photosynthesis theory within the FATES canopy structure.
--------------------------------------------------------------------------

The photosynthesis scheme is modified from the CLM4.5 model to give
estimates of photosynthesis, respiration and stomatal conductance for a
three dimenstional matrix indexed by canopy level (:math:`C_l`), plant
functional type (:math:`ft`) and leaf layer (:math:`z`). We conduct the
photosynthesis calculations at each layer for both sunlit and shaded
leaves. Thus, the model also generates estimates of :math:`w_{c},w_{j}`
and :math:`w_{e}` indexed in the same three dimensional matrix. In this
implementation, some properties (stomatal conductance parameters,
top-of-canopy photosynthetic capacity) vary with plant functional type,
and some vary with both functional type and canopy depth (absorbed
photosynthetically active radiation, nitrogen-based variation in
photosynthetic properties). The remaining drivers of photosynthesis
(:math:`P_{atm}`, :math:`K_c`, :math:`o_i`, :math:`K_o`, temperature,
atmospheric CO\ :math:`_2`) remain the same throughout the canopy. The
rate of gross photosynthesis (:math:`gpp_{Cl,ft,z}`)is the smoothed
minimum of the three potentially limiting processes (carboxylation,
electron transport, export limitation), but calculated independently for
each leaf layer:

.. math:: \textrm{gpp}_{Cl,ft,z} = \rm{min}(w_{c,Cl,ft,z},w_{j,Cl,ft,z},w_{e,Cl,ft,z}).

For :math:`w_{c,Cl,ft,z},`, we use

.. math::

   w_{c,Cl,ft,z}=  \left\{ \begin{array}{ll}
   \frac{V_{c,max,Cl,ft,z}(c_{i,Cl,ft,z}- \Gamma_*)}{c_{i,Cl,ft,z}+K_c(1+o_i/K_o)} & \textrm{for $C_3$ plants}\\
   &\\
   V_{c,max,Cl,ft,z}& \textrm{for $C_4$ plants}\\
   \end{array} \right.
   c_{i,Cl,ft,z}-\Gamma_*\ge 0

where :math:`V_{c,max}` now varies with PFT, canopy depth and layer
(see below). Internal leaf :math:`CO_{2}` (:math:`c_{i,Cl,ft,z})` is
tracked seperately for each leaf layer. For the light limited rate
:math:`w_j`, we use

.. math::

   w_j=  \left\{ \begin{array}{ll}
   \frac{J(c_i - \Gamma_*)4.6\phi\alpha}{4ci+8\Gamma_*} & \textrm{for C$_3$ plants}\\
   &\\
   4.6\phi\alpha & \textrm{for C$_4$ plants}\\
   \end{array} \right.

where :math:`J` is calculated as above but based on the absorbed
photosynthetically active radiation( :math:`\phi_{Cl,ft,z}`) for either
sunlit or shaded leaves in Wm\ :math:`^{-2}`. Specifically,

.. math::

   \phi_{Cl,ft,z}=  \left\{ \begin{array}{ll}
   abs_{sun,Cl,ft,z}& \textrm{for sunlit leaves}\\
   &\\
   abs_{sha,Cl,ft,z}& \textrm{for shaded leaves}\\
   \end{array} \right.

The export limited rate of carboxylation for C3 plants and the PEP
carboxylase limited rate of carboxylation for C4 plants :math:`w_c`
(also in :math:`\mu`\ mol CO\ :math:`_2` m\ :math:`^{-2}`
s\ :math:`^{-1}`) is calculated in a similar fashion,

.. math::

   w_{e,Cl,ft,z}=  \left\{ \begin{array}{ll}
   0.5V_{c,max,Cl,ft,z} & \textrm{for $C_3$ plants}\\
   &\\
   4000 V_{c,max,Cl,ft,z} \frac{c_{i,Cl,ft,z}}{P_{atm}}& \textrm{for $C_4$ plants}.\\
   \end{array} \right.

Variation in plant physiology with canopy depth
-----------------------------------------------

Both :math:`V_{c,max}` and :math:`J_{max}` vary with vertical depth in
the canopy on account of the well-documented reduction in canopy
nitrogen through the leaf profile, see :ref:`Bonan et al. 2012<bonanetal2012>` for
details). Thus, both :math:`V_{c,max}` and :math:`J_{max}` are indexed
by by :math:`C_l`, :math:`ft` and :math:`z` according to the nitrogen
decay coefficient :math:`K_n` and the amount of vegetation area shading
each leaf layer :math:`V_{above}`,

.. math::

   \begin{array}{ll}
   V_{c,max,Cl,ft,z} & = V_{c,max0,ft} e^{-K_{n,ft}V_{above,Cl,ft,z}},\\
   J_{max,Cl,ft,z} & = J_{max0,ft} e^{-K_{n,ft}V_{above,Cl,ft,z}},\\
   \end{array}

where :math:`V_{c,max,0}` and :math:`J_{max,0}` are the top-of-canopy
photosynthetic rates. :math:`V_{above}` is the sum of exposed leaf area
index (:math:`\textrm{elai}_{Cl,ft,z}`) and the exposed stem area index
(:math:`\textrm{esai}_{Cl,ft,z}`)( m\ :math:`^{2}` m\ :math:`^{-2}` ).
Namely,

.. math:: V_{Cl,ft,z} = \textrm{elai}_{Cl,ft,z} + \textrm{esai}_{Cl,ft,z}.

The vegetation index shading a particular leaf layer in the top canopy
layer is equal to

.. math::

   \begin{array}{ll}
   V_{above,Cl,ft,z}= \sum_{1}^{z} V_{Cl,ft,z} & \textrm{for $Cl= 1$. }
   \end{array}

For lower canopy layers, the weighted average vegetation index of the
canopy layer above (:math:`V_{canopy}`) is added to this within-canopy
shading. Thus,

.. math::

   \begin{array}{ll}
   V_{above,Cl,ft,z}=  \sum_{1}^{z}  V_{Cl,ft,z} + V_{canopy,Cl-1} & \textrm{for $Cl >1$, }\\
   \end{array}

where :math:`V_{canopy}` is calculated as

.. math:: V_{canopy,Cl} =  \sum_{ft=1}^{\emph{npft}} {\sum_{z=1}^{nz(ft)} (V_{Cl,ft,z} \cdot  \it{pft}_{wt,Cl,ft,1}).}

:math:`K_{n}` is the coefficient of nitrogen decay with canopy depth.
The value of this parameter is taken from the work of
:ref:`Lloyd et al. 2010<Lloydetal2010>` who determined, from 204 vertical profiles
of leaf traits, that the decay rate of N through canopies of tropical
rainforests was a function of the :math:`V_{cmax}` at the top of the
canopy. They obtain the following term to predict :math:`K_{n}`,

.. math:: K_{n,ft} = e^{0.00963 V_{c,max0,ft} - 2.43},

where :math:`V_{cmax}` is again in :math:`\mu`\ mol CO\ :math:`_2`
m\ :math:`^{-2}` s\ :math:`^{-1}`.

Water Stress on gas exchange
----------------------------

The top of canopy leaf photosynthetic capacity, :math:`V_{c,max0}`, is
also adjusted for the availability of water to plants as

.. math:: V_{c,max0,25} = V_{c,max0,25}  \beta_{sw},

where the adjusting factor :math:`\beta_{sw}` ranges from one when the
soil is wet to zero when the soil is dry. It depends on the soil water
potential of each soil layer, the root distribution of the plant
functional type, and a plant-dependent response to soil water stress,

.. math:: \beta_{sw} = \sum_{j=1}^{nj}w_{j}r_{j},

where :math:`w_{j}` is a plant wilting factor for layer :math:`j` and
:math:`r_{j}` is the fraction of roots in layer :math:`j`.The plant
wilting factor :math:`w_{j}` is

.. math::

   w_{j}=  \left\{ \begin{array}{ll}
   \frac{\psi_c-\psi_{j}}{\psi_c - \psi_o} (\frac{\theta_{sat,j} - \theta_{ice,j}}{\theta_{sat,j}})& \textrm{for $T_i >$-2C}\\
   &\\
   0 & \textrm{for $T_{j} \ge$-2C}\\
   \end{array} \right.

where :math:`\psi_{i}` is the soil water matric potential (mm) and
:math:`\psi_{c}` and :math:`\psi_{o}` are the soil water potential (mm)
when stomata are fully closed or fully open, respectively. The term in
brackets scales :math:`w_{i}` the ratio of the effective porosity (after
accounting for the ice fraction) relative to the total porosity.
:math:`w_{i}` = 0 when the temperature of the soil layer (:math:`T_{i}`
) is below some threshold (-2:math:`^{o}`\ C) or when there is no liquid
water in the soil layer (:math:`\theta_{liq,i} \le 0`). For more details
on the calculation of soil matric potential, see the CLM4.5 technical
note.

Variation of water stress and water uptake within tiles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The remaining drivers of the photosynthesis model remain constant
(atmospheric CO\ :math:`_2` and O\ :math:`^2` and canopy temperature)
throughout the canopy, except for the water stress index
:math:`\beta_{sw}`. :math:`\beta_{sw}` must be indexed by :math:`ft`,
because plants of differing functional types have the capacity to have
varying root depth, and thus access different soil moisture profile and
experience differing stress functions. Thus, the water stress function
applied to gas exchange calculation is now calculated as

.. math:: \beta_{sw,ft} = \sum_{j=1}^{nj}w_{j,ft} r_{j,ft},

where :math:`w_{j}` is the water stress at each soil layer :math:`j`
and :math:`r_{j,ft}` is the root fraction of each PFT’s root mass in
layer :math:`j`. Note that this alteration of the :math:`\beta_{sw}`
parameter also necessitates recalculation of the vertical water
extraction profiles. In the original model, the fraction of extraction
from each layer (:math:`r_{e,j,patch}`) is the product of a single root
distribution, because each patch only has one plant functional type. In
FATES, we need to calculate a new weighted patch effective rooting
depth profile :math:`r_{e,j,patch}` as the weighted average of the
functional-type level stress functions and their relative contributions
to canopy conductance. Thus for each layer :math:`j`, the extraction
fraction is summed over all PFTs as

.. math:: r_{e,j,patch} =  \sum_{ft=1}^{ft=npft} \frac{w_{j,ft}}{\sum_{j=1}^{=nj} w_{j,ft} }\frac{G_{s,ft}}{G_{s,canopy}},

where :math:`nj` is the number of soil layers, :math:`G_{s,canopy}`\ is
the total canopy (see section 9 for details) and :math:`G_{s,ft}` is the
canopy conductance for plant functional type :math:`ft`,

.. math:: G_{s,ft}= \sum_{1}w_{ncoh,ft} {gs_{can,coh} n_{coh} }.

Aggregation of assimilated carbon into cohorts
----------------------------------------------

The derivation of photosynthetic rates per leaf layer, as above, give us
the estimated rate of assimilation for a unit area of leaf at a given
point in the canopy in :math:`\mu`\ mol CO\ :math:`_2` m\ :math:`^{-2}`
s\ :math:`_{-1}`. To allow the integration of these rates into fluxes
per individual tree, or cohort of trees (gCO:math:`_2`
tree\ :math:`^{-1}` s\ :math:`^{-1}`), they must be multiplied by the
amount of leaf area placed in each layer by each cohort. Each cohort is
described by a single functional type, :math:`ft` and canopy layer
:math:`C_l` flag, so the problem is constrained to integrating these
fluxes through the vertical profile (:math:`z`).

We fist make a weighted average of photosynthesis rates from sun
(:math:`\textrm{gpp}_{sun}`, :math:`\mu`\ mol CO\ :math:`_2`
m\ :math:`^{-2}` s\ :math:`^{-1}`) and shade leaves (
:math:`\textrm{gpp}_{shade}`, :math:`\mu`\ mol CO\ :math:`_2`
m\ :math:`^{-2}` s\ :math:`^{-1}`) as

.. math:: \textrm{gpp}_{Cl,ft,z} =\textrm{gpp}_{sun,Cl,ft,z} f_{sun,Cl,ft,z}+ \textrm{gpp}_{sha,Cl,ft,z}(1-f_{sun,Cl,ft,z}).

The assimilation per leaf layer is then accumulated across all the leaf
layers in a given cohort (*coh*) to give the cohort-specific gross
primary productivity (:math:`\mathit{GPP}_{coh}`),

.. math:: \textit{GPP}_{coh} = 12\times 10^{-9}\sum_{z=1}^{nz(coh)}gpp_{Cl,ft,z} A_{crown,coh} \textrm{elai}_{Cl,ft,z}

The :math:`\textrm{elai}_{l,Cl,ft,z}` is the exposed leaf area which is
present in each leaf layer in m\ :math:`^{2}` m\ :math:`^{-2}`. (For all
the leaf layers that are completely occupied by a cohort, this is the
same as the leaf fraction of :math:`\delta_{vai}`). The fluxes are
converted from :math:`\mu`\ mol into mol and then multiplied by 12 (the
molecular weight of carbon) to give units for GPP\ :math:`_{coh}` of KgC
cohort\ :math:`^{-1}` s\ :math:`^{-1}`. These are integrated for each
timestep to give KgC cohort\ :math:`^{-1}` day\ :math:`^{-1}`

.. raw:: latex

   \captionof{table}{Parameters needed for photosynthesis model.}

+-----------------+-----------------+-----------------+-----------------+
| Parameter       | Parameter Name  | Units           | indexed by      |
| Symbol          |                 |                 |                 |
+=================+=================+=================+=================+
| :math:`V_{c,max | Maximum         | :math:`\mu` mol | *ft*            |
| 0}`             | carboxylation   | CO :math:`_2`   |                 |
|                 | capacity        | m :math:`^{-2}` |                 |
|                 |                 | s :math:`^{-1}` |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`r_b`     | Base Rate of    | gC              |                 |
|                 | Respiration     | gN\ :math:`^{-1 |                 |
|                 |                 | } s^{-1}`)      |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`q_{10}`  | Temp. Response  |                 |                 |
|                 | of stem and     |                 |                 |
|                 | root            |                 |                 |
|                 | respiration     |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`R_{cn,le | CN ratio of     | gC/gN           | *ft*            |
| af,ft}`         | leaf matter     |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`R_{cn,ro | CN ratio of     | gC/gN           | *ft*            |
| ot,ft}`         | root matter     |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`f_{gr}`  | Growth          | none            |                 |
|                 | Respiration     |                 |                 |
|                 | Fraction        |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`\psi_c`  | Water content   | Pa              | *ft*            |
|                 | when stomata    |                 |                 |
|                 | close           |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`\psi_o`  | Water content   | Pa              | *ft*            |
|                 | above which     |                 |                 |
|                 | stomata are     |                 |                 |
|                 | open            |                 |                 |
+-----------------+-----------------+-----------------+-----------------+

.. raw:: latex

   \bigskip 

Plant respiration
^^^^^^^^^^^^^^^^^^

Plant respiration per individual :math:`R_{plant,coh}` (KgC individual
:math:`^{-1}` s\ :math:`^{-1}`) is the sum of two terms, growth and
maintenance respiration :math:`R_{g,coh}` and :math:`R_{m,coh}`

.. math:: R_{plant} = R_{g,coh}+ R_{m,coh}

Maintenance respiration is the sum of the respiration terms from four
different plant tissues, leaf, :math:`R_{m,leaf,coh}`, fine root
:math:`R_{m,froot,coh}`, coarse root :math:`R_{m,croot,coh}`\ and stem
:math:`R_{m,stem,coh}`, all also in (KgC individual :math:`^{-1}`
s\ :math:`^{-1}`) .

.. math:: R_{m,coh} = R_{m,leaf,coh}+ R_{m,froot,coh}+R_{m,croot,coh}+R_{m,stem,coh}

To calculate canopy leaf respiration, which varies through we canopy, we
first determine the top-of-canopy leaf respiration rate
(:math:`r_{m,leaf,ft,0}`, gC s\ :math:`^{-1}` m\ :math:`^{-2}`) is
calculated from a base rate of respiration per unit leaf nitrogen
derived from :ref:`Ryan et al. 1991<ryan1991>`. The base rate for leaf
respiration (:math:`r_{b}`) is 2.525 gC/gN s\ :math:`^{-1}`,

.. math:: r_{m,leaf,ft,0} = r_{b} N_{a,ft}(1.5^{(25-20)/10})

where :math:`r_b` is the base rate of metabolism (2.525 x
10\ :math:`^6` gC/gN s\ :math:`^{-1}`. This base rate is adjusted
assuming a Q\ :math:`_{10}` of 1.5 to scale from the baseline of 20C to
the CLM default base rate temperature of 25C. For use in the
calculations of net photosynthesis and stomatal conductance, leaf
respiration is converted from gC s\ :math:`^{-1}` m\ :math:`^{-2}`, into
:math:`\mu`\ mol CO\ :math:`_2` m\ :math:`^{-2}` s\ :math:`^{-1}`
(:math:`/12\cdot 10^{-6}`).

This top-of-canopy flux is scaled to account for variation in
:math:`N_a` through the vertical canopy, in the same manner as the
:math:`V_{c,max}` values are scaled using :math:`V_{above}`.

.. math:: r_{leaf,Cl,ft,z}  = r_{m,leaf,ft,0} e^{-K_{n,ft}V_{above,Cl,ft,z}}\beta_{ft}f(t)

Leaf respiration is also adjusted such that it is reduced by drought
stress, :math:`\beta_{ft}`, and canopy temperature, :math:`f(t_{veg})`.
For details of the temperature functions affecting leaf respiration see
the CLM4 technical note, Section 8, Equations 8.13 and 8.14. The
adjusted leaf level fluxes are scaled to individual-level (gC individual
:math:`^{-1}` s\ :math:`^{-1}`) in the same fashion as the
:math:`\rm{GPP}_{coh}` calculations

.. math:: \rm{R}_{m,leaf,coh} = 12\times 10^{-9}\sum_{z=1}^{nz(coh)}r_{leaf,Cl,ft,z} A_{crown} \textrm{elai}_{Cl,ft,z}

The stem and the coarse-root respiration terms are derived using the
same base rate of respiration per unit of tissue Nitrogen.

.. math:: R_{m,croot,coh} =  10^{-3}r_b t_c \beta_{ft} N_{\rm{livecroot,coh}}

.. math:: R_{m,stem,coh} =   10^{-3}r_b t_c \beta_{ft} N_{\rm{stem,coh}}

Here, :math:`t_c` is a temperature relationship based on a
:math:`q_{10}` value of 1.5, where :math:`t_v` is the vegetation
temperature. We use a base rate of 20 here as, again, this is the
baseline temperature used by :ref:`Ryan et al. 1991<ryan1991>`. The
10\ :math:`^{-3}` converts from gC invididual\ :math:`^{-1}`
s\ :math:`^{-1}` to KgC invididual\ :math:`^{-1}` s\ :math:`^{-1}`

.. math:: t_c=q_{10}^{(t_{v} - 20)/10}

The tissue N contents for live sapwood are derived from the leaf CN
ratios, and for fine roots from the root CN ratio as:

.. math:: N_{\rm{stem,coh}}  = \frac{B_{\rm{sapwood,coh}}}{ R_{cn,leaf,ft}}

and

.. math:: N_{\rm{livecroot,coh}}  = \frac{ B_{\rm{root,coh}}w_{frac,ft}}{R_{cn,root,ft}}

where :math:`B_{\rm{sapwood,coh}}` and :math:`B_{\rm{root,coh}}` are
the biomass pools of sapwood and live root biomass respectively (KgC
individual) and :math:`w_{frac,ft}` is the fraction of coarse root
tissue in the root pool (0.5 for woody plants, 0.0 for grasses and
crops). We assume here that stem CN ratio is the same as the leaf C:N
ratio, for simplicity. The final maintenance respiration term is derived
from the fine root respiration, which accounts for gradients of
temperature in the soil profile and thus calculated for each soil layer
:math:`j` as follows:

.. math:: R_{m,froot,j } = \frac{(1 - w_{frac,ft})B_{\rm{root,coh}}b_r\beta_{ft}}{10^3R_{cn,leaf,ft}}   \sum_{j=1}^{nj}t_{c,soi,j} r_{i,ft,j}

:math:`t_{c,soi}` is a function of soil temperature in layer :math:`j`
that has the same form as that for stem respiration, but uses vertically
resolved soil temperature instead of canopy temperature. In the CLM4.5,
only coarse and not fine root respriation varies as a function of soil
depth, and we maintain this assumption here, although it may be altered
in later versions. The growth respiration, :math:`R_{g,coh}` is a fixed
fraction :math:`f_{gr}` of the carbon remaining after maintenance
respiration has occurred.

.. math:: R_{g,coh}=\textrm{max}(0,GPP_{g,coh} - \it R\rm_{m,coh})f_{gr}

.. raw:: latex

   \captionof{table}{Parameters needed for plant respiration model.  }

+-----------------+-----------------+-----------------+-----------------+
| Parameter       | Parameter Name  | Units           | indexed by      |
| Symbol          |                 |                 |                 |
+=================+=================+=================+=================+
| :math:`-K_{n,ft | Rate of         | none            | -               |
| }`              | reduction of N  |                 |                 |
|                 | through the     |                 |                 |
|                 | canopy          |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`r_b`     | Base Rate of    | gC              |                 |
|                 | Respiration     | gN\ :math:`^{-1 |                 |
|                 |                 | } s^{-1}`)      |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`q_{10}`  | Temp. Response  |                 |                 |
|                 | of stem and     |                 |                 |
|                 | root            |                 |                 |
|                 | respiration     |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`R_{cn,le | CN ratio of     | gC/gN           | *ft*            |
| af,ft}`         | leaf matter     |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`R_{cn,ro | CN ratio of     | gC/gN           | *ft*            |
| ot,ft}`         | root matter     |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`f_{gr}`  | Growth          | none            |                 |
|                 | Respiration     |                 |                 |
|                 | Fraction        |                 |                 |
+-----------------+-----------------+-----------------+-----------------+

.. raw:: latex

   \bigskip 

Stomatal Conductance
^^^^^^^^^^^^^^^^^^^^

Fundamental stomatal conductance theory
---------------------------------------

Stomatal conductance is unchanged in concept from the CLM4.5 approach.
Leaf stomatal resistance is calculated from the Ball-Berry conductance
model as described by :ref:`Collatz et al. (1991)<Collatzetal1991>` and implemented in
a global climate model by :ref:`Sellers et al. 1996<sellers1996>`. The model
relates stomatal conductance (i.e., the inverse of resistance) to net
leaf photosynthesis, scaled by the relative humidity at the leaf surface
and the CO\ :math:`_2` concentration at the leaf surface. The primary
difference between the CLM implementation and that used by
:ref:`Collatz et al. (1991)<Collatzetal1991>` and :ref:`Sellers et al. (1996)<sellers1996>` is
that they used net photosynthesis (i.e., leaf photosynthesis minus leaf
respiration) instead of gross photosynthesis. As implemented here,
stomatal conductance equals the minimum conductance (:math:`b`) when
gross photosynthesis (:math:`A`) is zero. Leaf stomatal resistance is

.. math:: \frac{1}{r_{s}} = m_{ft} \frac{A}{c_s}\frac{e_s}{e_i}P_{atm}+b_{ft} \beta_{sw}

where :math:`r_{s}` is leaf stomatal resistance (s m\ :math:`^2`
:math:`\mu`\ mol\ :math:`^{-1}`), :math:`b_{ft}` is a plant functional
type dependent parameter equivalent to :math:`g_{0}` in the Ball-Berry
model literature. This parameter is also scaled by the water stress
index :math:`\beta_{sw}`. Similarly, :math:`m_{ft}` is the slope of the
relationship between the assimilation, :math:`c_s` and humidty dependant
term and the stomatal conductance, and so is equivalent to the
:math:`g_{1}` term in the stomatal literature. :math:`A` is leaf
photosynthesis (:math:`\mu`\ mol CO\ :math:`_2` m\ :math:`^{-2}`
s\ :math:`^{-1}`), :math:`c_s` is the CO\ :math:`_2` partial pressure at
the leaf surface (Pa), :math:`e_s` is the vapor pressure at the leaf
surface (Pa), :math:`e_i` is the saturation vapor pressure (Pa) inside
the leaf at the vegetation temperature conductance (:math:`\mu`\ mol
m\ :math:`^{-2}` s\ :math:`^{-1}`) when :math:`A` = 0 . Typical values
are :math:`m_{ft}` = 9 for C\ :math:`_3` plants and :math:`m_{ft}` = 4
for C\ :math:`_4` plants (
:ref:`Collatz et al. 1991<Collatzetal1991>`, :ref:`Collatz, 1992<Collatzetal1992>`, :ref:`Sellers et al 1996<sellersetal1996>`).
:ref:`Sellers et al. 1996<sellers1996>` used :math:`b` = 10000 for C\ :math:`_3`
plants and :math:`b` = 40000 for C\ :math:`_4` plants. Here, :math:`b`
was chosen to give a maximum stomatal resistance of 20000 s
m\ :math:`^{-1}`. These terms are nevertheless plant strategy dependent,
and have been found to vary widely with plant type
:ref:`Medlyn et al. 2011<Medlynetal2011>`.

Resistance is converted from units of s m\ :math:`^2 \mu`
mol\ :math:`^{-1}` to s m\ :math:`^{-1}` as: 1 s m\ :math:`^{-1}` =
:math:`1\times 10^{-9}`\ R\ :math:`_{\rm{gas}} \theta_{\rm{atm}}P_{\rm{atm}}`
(:math:`\mu`\ mol\ :math:`^{-1}` m\ :math:`^{2}` s), where
R\ :math:`_{gas}` is the universal gas constant (J K\ :math:`^{-1}`
kmol\ :math:`^{-1}`) and :math:`\theta_{atm}` is the atmospheric
potential temperature (K).

Resolution of stomatal conductance theory in the FATES canopy structure
-----------------------------------------------------------------------

The stomatal conductance is calculated, as with photosynthesis, for each
canopy, PFT and leaf layer. The CLM code requires a single canopy
conductance estimate to be generated from the multi-layer multi-PFT
array. In previous iterations of the CLM, sun and shade-leaf specific
values have been reported and then averaged by their respective leaf
areas. In this version, the total canopy condutance
:math:`G_{s,canopy}`, is calculated as the sum of the cohort-level
conductance values.

.. math:: G_{s,canopy} =  \sum{ \frac{gs_{can,coh} n_{coh} }{A_{patch}}}

Cohort conductance is the sum of the inverse of the leaf resistances at
each canopy layer (:math:`r_{s,z}` ) multipled by the area of each
cohort.

.. math:: gs_{can,coh} =\sum_{z=1}^{z=nv,coh}{\frac{ A_{crown,coh}}{r_{s,cl,ft,z}+r_{b}}}


.. raw:: latex

   \captionof{table}{Parameters needed for stomatal conductance model.  }

+------------------+--------------------------+-------+------------+
| Parameter Symbol | Parameter Name           | Units | indexed by |
+==================+==========================+=======+============+
| :math:`b_{ft}`   | Slope of Ball-Berry term | none  | *ft*       |
+------------------+--------------------------+-------+------------+
| :math:`m_{ft}`   | Slope of Ball-Berry term | none  | *ft*       |
+------------------+--------------------------+-------+------------+


Allocation and Growth
^^^^^^^^^^^^^^^^^^^^^

Total assimilation carbon enters the ED model each day as a
cohort-specific Net Primary Productivity :math:`\mathit{NPP}_{coh}`,
which is calculated as

.. math:: \mathit{NPP}_{coh} = \mathit{GPP}_{coh} - R_{plant,coh}

This flux of carbon is allocated between the demands of tissue turnover,
of carbohydrate storage and of growth (increase in size of one or many
plant organs). Priority is explicitly given to maintenance respiration,
followed by tissue maintenance and storage, then allocation to live
biomass and then to the expansion of structural and live biomass pools.
All fluxes here are first converted into in KgC
individual\ :math:`^{-1}` year\ :math:`^{-1}` and ultimately integrated
using a timesteps of 1/365 years for each day.

Tissue maintenance demand
-------------------------

We calculate a ‘tissue maintenance’ flux. The magnitude of this flux is
such that the quantity of biomass in each pool will remain constant,
given background turnover rates. For roots, this maintenenace demand is
simply

.. math:: r_{md,coh}  = b_{root}\cdot\alpha_{root,ft}

Where :math:`\alpha_{root,ft}` is the root turnover rate in y^-1. Given
that, for deciduous trees, loss of leaves is assumed to happen only one
per growing season, the algorithm is dependent on phenological habit
(whether or not this PFT is evergreen), thus

.. math::

   l_{md,coh} = \left\{ \begin{array}{ll}
    b_{leaf}\cdot\alpha_{leaf,ft}&\textrm{for } P_{evergreen}= 1\\
   &\\
   0&\textrm{for }  P_{evergreen}= 0\\
   \end{array} \right.

Leaf litter resulting from deciduous senescence is handled in the
phenology section. The total quantity of maintenance demand
(:math:`t_{md,coh}`. KgC individual y\ :math:`^{-1}`) is therefore

.. math:: t_{md,coh}  = l_{md,coh} + r_{md,coh}

Allocation to storage and turnover
----------------------------------

The model must now determine whether the NPP input is sufficient to meet
the maintenance demand and keep tissue levels constant. To determine
this, we introduce the idea of ‘carbon balance’ :math:`C_{bal,coh}` (KgC
individual\ :math:`^{-1}`) where

.. math:: C_{bal,coh}= \mathit{NPP}_{coh} - t_{md,coh}\cdot f_{md,min,ft}

where :math:`f_{md,min,ft}` is the minimum fraction of the maintenance
demand that the plant must meet each timestep, which is indexed by *ft*
and represents a life-history-strategy decision concerning whether
leaves should remain on in the case of low carbon uptake (a risky
strategy) or not be replaced (a conservative strategy). Subsequently, we
determine a flux to the storage pool, where the flux into the pool, as a
fraction of :math:`C_{bal,coh}`, is proportional to the discrepancy
between the target pool size and the actual pool size :math:`f_{tstore}`
where

.. math:: f_{tstore} =  \mathrm{max}\left(0,\frac{b_{store}}{b_{leaf}\cdot S_{cushion}}\right)

The allocation to storage is a fourth power function of
:math:`f_{tstore}` to mimic the qualitative behaviour found for carbon
allocation in arabidopsis by :ref:`Smith et al. 2007<smith2007>`.

.. math::

   \frac{\delta b_{store}}{\delta t} = \left\{ \begin{array}{ll}
   C_{bal,coh} \cdot e^{-f_{tstore}^{4}} &\textrm{for }C_{bal,coh}>0\\
   &\\
   C_{bal,coh} &\textrm{for }C_{bal,coh}\leq0\\
   \end{array} \right.

If the carbon remaining after the storage and minimum turnover fluxes
have been met, the next priority is the remaining flux to leaves
:math:`t_{md}\cdot(1-f_{md,min})`. If the quantity of carbon left
:math:`(C_{bal,coh}-\frac{\delta b_{store}}{\delta t})` is insufficient
to supply this amount of carbon, then the store of alive carbon is
depleted (to represent those leaves that have fallen off and not been
replaced)

.. math::

   \frac{\delta b_{alive}}{\delta t} = \left\{ \begin{array}{ll}
   0 &\textrm{ for } (C_{bal,coh}-\frac{\delta b_{store}}{\delta t}) > t_{md}\cdot(1-f_{md,min})\\
   &\\
   t_{md}\cdot(1-f_{md,min}) - \left(C_{bal,coh}-\frac{\delta b_{store}}{\delta t}\right)&\textrm{ for } (C_{bal,coh}-\frac{\delta b_{store}}{\delta t}) \leq t_{md}\cdot(1-f_{md,min})\\
   \end{array} \right.

correspondingly, the carbon left over for growth (:math:`C_{growth}`:
(KgC individual\ :math:`^{-1}` year\ :math:`^{-1}`) is therefore

.. math::

   C_{growth} = \left\{ \begin{array}{ll}
   C_{bal,coh}-\frac{\delta b_{store}}{\delta t} &\textrm{ for } (C_{bal,coh}-\frac{\delta b_{store}}{\delta t}) > 0\\
   &\\
   0&\textrm{ for } (C_{bal,coh}-\frac{\delta b_{store}}{\delta t}) \leq 0\\
   \end{array} \right.

to allocate the remaining carbon (if there is any), we first ascertain
whether the live biomass pool is at its target, or whether is has been
depleted by previous low carbon timesteps. Thus

.. math::

   \begin{array}{lll}
   b_{alive,target}&= b_{leaf,target}  (1+ f_{frla}+f_{swh}h_{coh}) &\textrm{for } S_{phen,coh} = 2\\
   b_{alive,target}&= b_{leaf,target}  ( f_{frla}+f_{swh}h_{coh})&\textrm{for } S_{phen,coh} = 1\\
   \end{array}

where the target leaf biomass :math:`b_{leaf.target}` ((Kg C
individual\ :math:`^{-1}`)) is the allometric relationship between dbh
and leaf biomass, ameliorated by the leaf trimming fraction (see
‘control of leaf area’ below)

.. math:: b_{leaf.target} = c_{leaf}\cdot dbh_{coh}^{e_{leaf,dbh}} \rho_{ft} ^{e_{leaf,dens}}\cdot C_{trim,coh}

:math:`\rho_{ft}` is the wood density, in g cm\ :math:`^{3}`.

Allocation to Seeds
-------------------

The fraction remaining for growth (expansion of live and structural
tissues) :math:`f_{growth}` is 1 minus that allocated to seeds.

.. math:: f_{growth,coh} = 1 - f_{seed,coh}

Allocation to seeds only occurs if the alive biomass is not below its
target, and then is a predefined fixed fraction of the carbon remaining
for growth. Allocation to clonal reproduction (primarily for grasses)
occurs when :math:`\textrm{max}_{dbh}` is achieved.

.. math::

   f_{seed,coh} = \left\{ \begin{array}{ll}
   R_{frac,ft}&\textrm{ for } \textrm{max}_{dbh} < dbh_{coh} \\
   &\\
   \left( R_{frac,ft}+C_{frac,ft} \right) &\textrm{ for } \textrm{max}_{dbh} \geq dbh_{coh} \\
   \end{array} \right.

the total amount allocated to seed production (:math:`p_{seed,coh}` in
KgC individual :math:`^{-1}` y\ :math:`^{-1}`) is thus

.. math:: p_{seed,coh} = C_{growth}\cdot f_{seed,coh}

Allocation to growing pools
---------------------------

The carbon is then partitioned into carbon available to grow the
:math:`b_{alive}` and :math:`b_{struc}` pools. A fraction :math:`v_{a}`
is available to live biomass pools, and a fraction :math:`v_{s}` is
available to structural pools.

.. math:: \frac{\delta b_{alive}}{\delta t} = C_{growth}\cdot  f_{growth} v_{a}

.. math:: \frac{\delta b_{struc}}{\delta t} = C_{growth}\cdot  f_{growth} v_{s}

If the alive biomass is lower than its ideal target, all of the
available carbon is directed into that pool. Thus:

.. math::

   v_{a}= \left\{ \begin{array}{ll}
   \frac{1}{1+u}&\textrm{ for } b_{alive} \geq b_{alive,target} \\
   &\\
   1.0&\textrm{ for } b_{alive} <  b_{alive,target} \\
   \end{array} \right.

.. math::

   v_{s}= \left\{ \begin{array}{ll}
   \frac{u}{1+u}&\textrm{ for } b_{alive} \geq b_{alive,target} \\
   &\\
   0.0&\textrm{ for } b_{alive} <  b_{alive,target} \\
   \end{array} \right.

In this case, the division of carbon between the live and structural
pools :math:`u` is derived as the inverse of the sum of the rates of
change in live biomass with respect to structural:

.. math:: u = \frac{1}{\frac{\delta b_{leaf}}{ \delta b_{struc} } + \frac{\delta b_{root}}{ \delta b_{struc} } +\frac{\delta b_{sw}}{ \delta b_{struc} } }

To calculate all these differentials, we first start with
:math:`\delta b_{leaf}/\delta b_{struc}`, where

.. math:: \frac{\delta b_{leaf}}{ \delta b_{struc}}= \frac{\frac{\delta \mathrm{dbh}}{\delta b_{struc}}}   {\frac{\delta \mathrm{dbh} }{\delta b_{leaf}} }

The rates of change of dbh with respect to leaf and structural biomass
are the differentials of the allometric equations linking these terms to
each other. Hence,

.. math:: \frac{\delta \mathrm{dbh} }{\delta b_{leaf}}=\frac{1}{b_{trim,coh}}\cdot (e_{leaf,dbh}-1)\exp  {\big(c_{leaf} \mathrm{dbh}^{(e_{leaf,dbh})-1} \rho_{ft}^{e_{leaf,dens}} \big)}

and where :math:`\mathrm{dbh}_{coh} >   \mathrm{dbh}_{max}`

.. math:: \frac{\delta b_{struc}}{\delta \mathrm{dbh}}  = e_{str,dbh} \cdot c_{str}\cdot e_{str,hite} h_{coh}^{e_{str,dbh}-1}   \mathrm{dbh}_{coh}^{e_{str,dbh}} \rho_{ft}^{e_{str,dens}}

If :math:`\mathrm{dbh}_{coh} \leq   \mathrm{dbh}_{max}` then we must
also account for allocation for growing taller as:

.. math:: \frac{\delta b_{struc}}{\delta \mathrm{dbh}} =  \frac{\delta b_{struc}}{\delta \mathrm{dbh}} + \frac{\delta h}{\delta \mathrm{dbh}} \cdot  \frac{\delta b_{struc}}{\delta \mathrm{dbh} }

where

.. math:: \frac{\delta h}{\delta \mathrm{dbh}}= 1.4976  \mathrm{dbh}_{coh}^{m_{allom}-1}

.. math:: \frac{\delta  \mathrm{dbh} }{\delta b_{struc}} =\frac{1}{ \frac{\delta b_{struc}}{\delta \mathrm{dbh}} }

Once we have the :math:`\delta b_{leaf}/\delta b_{struc}`, we calculate
:math:`\delta b_{root}/\delta b_{struc}` as

.. math:: \frac{\delta  b_{root}}{\delta b_{struc}} =\frac{\delta b_{leaf}}{\delta b_{struc}}\cdot f_{frla}

and the sapwood differential as

.. math:: \frac{\delta  b_{sw}}{\delta b_{struc}} = f_{swh}\left( h_{coh} \frac{\delta b_{leaf}}{ \delta b_{struc}} + b_{leaf,coh}\frac{\delta h}{\delta b_{struc}} \right)

where

.. math:: \frac{\delta h}{\delta b_{struc}} =  \frac{1}{c_{str}\times e_{str,hite} h_{coh}^{e_{str,dbh}-1}   \mathrm{dbh}_{coh}^{e_{str,dbh}} \rho_{ft}^{e_{str,dens}}}

In all of the above terms, height in in m, :math:`\mathrm{dbh}` is in
cm, and all biomass pools are in KgCm\ :math:`^{-2}`. The allometric
terms for the growth trajectory are all taken from the ED1.0 model, but
could in theory be altered to accomodate alternative allometric
relationships. Critically, the non-linear relationships between live and
structural biomass pools are maintained in this algorithm, which
diverges from the methodology currently deployed in the CLM4.5.

Integration of allocated fluxes
-------------------------------

All of the flux calculations generate differential of the biomass state
variables against time (in years). To integrate these differential rates
into changes in the state variables, we use a simple simple forward
Euler integration. Other methods exist (e.g. ODEINT solvers, Runge Kutta
methods etc.), but they are more prone to errors that become difficult
to diagnose, and the typically slow rates of change of carbon pools mean
that these are less important than they might be in strongly non-linear
systems (soil drainage, energy balance, etc.)

.. math:: b_{alive,t+1} = \textrm{min}\left( 0,b_{alive,t} +  \frac{\delta b_{alive}}{\delta t}  \delta t \right)

.. math:: b_{struc,t+1} = \textrm{min}\left(0, b_{struc,t} +  \frac{\delta b_{struc}}{\delta t}  \delta t \right)

.. math:: b_{store,t+1} = \textrm{min}\left(0, b_{store,t} +  \frac{\delta b_{store}}{\delta t}  \delta t \right)

In this case, :math:`\delta t` is set to be one day
(:math:`\frac{1}{365}` years).

.. raw:: latex

   \bigskip

.. raw:: latex

   \captionof{table}{Parameters needed for allocation model. }

+-----------------+-----------------+-----------------+-----------------+
| Parameter       | Parameter Name  | Units           | indexed by      |
| Symbol          |                 |                 |                 |
+=================+=================+=================+=================+
| S               | Target stored   | none            | *ft*            |
|                 | biomass as      |                 |                 |
|                 | fraction of     |                 |                 |
|                 | :math:`b_{leaf}`|                 |                 |
|                 |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| f               | Minimum         | none            | *ft*            |
|                 | fraction of     |                 |                 |
|                 | turnover that   |                 |                 |
|                 | must be met     |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| R               | Fraction        | none            | *ft*            |
|                 | allocated to    |                 |                 |
|                 | seeds           |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| C               | Fraction        | none            | *ft*            |
|                 | allocated to    |                 |                 |
|                 | clonal          |                 |                 |
|                 | reproduction    |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`\textrm{ | Diameter at     | m               | *ft*            |
| max}_{dbh}`     | which maximum   |                 |                 |
|                 | height is       |                 |                 |
|                 | achieved        |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| P               | Does this       | 1=yes, 0=no     | *ft*            |
|                 | cohort have an  |                 |                 |
|                 | evergreen       |                 |                 |
|                 | phenological    |                 |                 |
|                 | habit?          |                 |                 |
+-----------------+-----------------+-----------------+-----------------+

.. raw:: latex

   \bigskip 

Control of Leaf Area Index
^^^^^^^^^^^^^^^^^^^^^^^^^^

The leaf area :math:`A_{leaf}` (m:math:`^{-2}`) of each cohort is
calculated from leaf biomass :math:`b_{leaf,coh}` (kgC
individual\ :math:`^{-1}`) and specific leaf area (SLA, m\ :math:`^2` kg
C\ :math:`^{-1}`)

.. math:: A_{leaf,coh} = b_{leaf,coh} \cdot SLA_{ft}

For a given tree allometry, leaf biomass is determined from basal area
using the function used by :ref:`Moorcroft et al. 2001<mc_2001>` where :math:`d_w`
is wood density in g cm\ :math:`^{-3}`.

.. math:: b_{leaf,coh} = c_{leaf} \cdot dbh_{coh}^{e_{leaf,dbh}} \rho_{ft}^{e_{leaf,dens}}

However, using this model, where leaf area and crown area are both
functions of diameter, the leaf area index of each tree in a closed
canopy forest is always the same (where :math:`S_{c,patch}` =
:math:`S_{c,min}` , irrespective of the growth conditions. To allow
greater plasticity in tree canopy structure, and for tree leaf area
index to adapt to prevailing conditions, we implemented a methodology
for removing those leaves in the canopy that exist in negative carbon
balance. That is, their total annual assimilation rate is insufficient
to pay for the turnover and maintenance costs associated with their
supportive root and stem tissue, plus the costs of growing the leaf. The
tissue turnover maintenance cost (KgC m\ :math:`^{-2} y^{-1}` of leaf is
the total maintenance demand divided by the leaf area:

.. math:: L_{cost,coh} = \frac{t_{md,coh}} {b_{leaf,coh} \cdot \textrm{SLA}}

The net uptake for each leaf layer :math:`U_{net,z}` in (KgC
m\ :math:`^{-2}` year\ :math:`^{-1}`) is

.. math:: U_{net,coh,z} = g_{coh,z}-r_{m,leaf,coh,z}

where :math:`g_{z}` is the GPP of each layer of leaves in each tree (KgC
m\ :math:`^{-2}` year\ :math:`^{-1}`), :math:`r_{m,leaf,z}` is the rate
of leaf dark respiration (also KgC m\ :math:`^{-2}`
year\ :math:`^{-1}`). We use an iterative scheme to define the cohort
specific canopy trimming fraction :math:`C_{trim,coh}`, on an annual
time-step, where

.. math:: b_{leaf,coh} =   C_{trim} \times 0.0419  dbh_{coh}^{1.56} d_w^{0.55}

If the annual maintenance cost of the bottom layer of leaves (KgC m-2
year-1) is less than then the canopy is trimmed by an increment :math:`\iota_l`\ (0.01), which is applied until the end of  next calander year. Because this is an optimality model, there is an
issue of the timescale over which net assimilation is evaluated, the
timescale of response, and the plasticity of plants to respond to these
pressures. These properties should be investigated further in future
efforts.

.. math::

   C_{trim,y+1}  = \left\{ \begin{array}{ll}
   \rm{max}(C_{trim,y}-\iota_l,1.0)&\rm{for} (L_{cost,coh} > U_{net,coh,nz})\\
   &\\
   \rm{min}(C_{trim,y}+\iota_l,L_{trim,min})&\rm{for} (L_{cost,coh} < U_{net,coh,nz})\\
   \end{array} \right.

We impose an arbitrary minimum value on the scope of canopy trimming of
:math:`L_{trim,min}` (0.5). If plants are able simply to drop all of
their canopy in times of stress, with no consequences, then tree
mortality from carbon starvation is much less likely to occur because of
the greatly reduced maintenance and turnover requirements.

.. raw:: latex

   \bigskip

.. raw:: latex

   \captionof{table}{Parameters needed for leaf area control model.  }

+-----------------+-----------------+-----------------+-----------------+
| Parameter       | Parameter Name  | Units           | indexed by      |
| Symbol          |                 |                 |                 |
+=================+=================+=================+=================+
| :math:`\iota_l` | Fraction by     | none            | -               |
|                 | which leaf mass |                 |                 |
|                 | is reduced next |                 |                 |
|                 | year            |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`L_{trim, | Minimum         | -               |                 |
| min}`           | fraction to     |                 |                 |
|                 | which leaf mass |                 |                 |
|                 | can be reduced  |                 |                 |
+-----------------+-----------------+-----------------+-----------------+

.. raw:: latex

   \bigskip 

Phenology
^^^^^^^^^^^^^^^^^^^^

Cold Deciduous Phenology
------------------------

Cold Leaf-out timing
~~~~~~~~~~~~~~~~~~~~

The phenology model of :ref:`Botta et al. 2000<botta2000>` is used in
FATES to determine the leaf-on timing. The Botta et al. model was
verified against satellite data and is one of the only globally verified
and published models of leaf-out phenology. This model differs from the
phenology model in the CLM4.5. The model simulates leaf-on date as a
function of the number of growing degree days (GDD), defined by the sum
of mean daily temperatures (:math:`T_{day}` :math:`^{o}`\ C) above a
given threshold :math:`T_{g}` (0 :math:`^{o}`\ C).

.. math:: GDD=\sum \textrm{max}(T_{day}-T_{g},0)

Budburst occurs when :math:`GDD` exceeds a threshold
(:math:`GDD_{crit}`). The threshold is modulated by the number of
chilling days experienced (NCD) where the mean daily temperature falls
below a threshold determined by `Botta et al. 2000<botta2000>` as
5\ :math:`^{o}`\ C. A greater number of chilling days means that fewer
growing degree days are required before budburst:

.. math:: GDD_{crit}=a+be^{c.NCD}

where a = -68, b= 638 and c=-0.01 `Botta et al. 2000<botta2000>`. In the
Northern Hemisphere, counting of degree days begins on 1st January, and
of chilling days on 1st November. The calendar opposite of these dates
is used for points in the Southern Hemisphere.

If the growing degree days exceed the critical threshold, leaf-on is
triggered by a change in the gridcell phenology status flag
:math:`S_{phen,grid}` where ‘2’ indicates that leaves should come on and
‘1’ indicates that they should fall.

.. math::

   \begin{array}{ll}
   S_{phen,grid} = 2
   &\textrm{ if } S_{phen,grid} = 1\textrm{ and } GDD_{grid} \ge GDD_{crit} \\
   \end{array}

Cold Leaf-off timing
~~~~~~~~~~~~~~~~~~~~

The leaf-off model is taken from the Sheffield Dynamic Vegetation Model
(SDGVM) and is similar to that for LPJ
:ref:`Sitch et al. 2003<sitch2003>` and IBIS
:ref:`Foley et al. 1996<Foley1996>` models. The average daily
temperatures of the previous 10 day period are stored. Senescence is
triggered when the number of days with an average temperature below
7.5\ :math:`^{o}` (:math:`n_{colddays}`) rises above a threshold values
:math:`n_{crit,cold}`, set at 5 days.

.. math::

   \begin{array}{ll}
   S_{phen,grid} = 1
   &\textrm{ if } S_{phen,grid} = 2\textrm{ and } n_{colddays} \ge n_{crit,cold} \\
   \end{array}

Global implementation modifications
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because of the global implementation of the cold-deciduous phenology
scheme, adjustments must be made to account for the possibility of
cold-deciduous plants experiencing situations where no chilling period
triggering leaf-off ever happens. If left unaccounted for, these leaves
will last indefinitely, resulting in highly unrealistic behaviour.
Therefore, we implement two additional rules. Firstly, if the number of
days since the last senescence event was triggered is larger than 364,
then leaf-off is triggered on that day. Secondly, if no chilling days
have occured during the winter accumulation period, then leaf-on is not
triggered. This means that in effect, where there are no cold periods,
leaves will fall off and not come back on, meaning that cold-deciduous
plants can only grow in places where there is a cold season.

Further to this rule, we introduce a ‘buffer’ time periods after leaf-on
of 30 days, so that cold-snap periods in the spring cannot trigger a
leaf senescence. The 30 day limit is an arbitrary limit. In addition, we
constrain growing degree day accumulation to the second half of the year
(Jult onwards in the Northern hemisphere, or Jan-June in the Southern)
and only allow GDD accumulation while the leaves are off.

Drought-deciduous Phenology: TBD
-------------------------------- 

In the current version of the model, a drought deciduous algorithm
exists, but is not yet operational, due to issue detected in the existing
CN and soil moisture modules, which also affect the behaviour of the
native ED drought deciduous model. This is a priority to address before
the science tag is released.

Carbon Dynamics of deciduous plants
----------------------------------- 

In the present version, leaf expansion and senescence happen over the
course of a single day. This is clearly not an empirically robust
representation of leaf behaviour, whereby leaf expansion occurs over a
period of 10-14 days, and senescence over a similar period. This will be
incorporated in later versions. When the cold or drought phenological
status of the gridcell status changes (:math:`S_{phen,grid}`) from ‘2’
to ‘1’, and the leaves are still on (:math:`S_{phen,coh}` =2 ), the leaf
biomass at this timestep is ’remembered’ by the model state variable
:math:`l_{memory,coh}`. This provides a ‘target’ biomass for leaf onset
at the beginning of the next growing season (it is a target, since
depletion of stored carbon in the off season may render achieving the
target impossible).

.. math:: l_{memory,coh} = b_{leaf,coh}

Leaf carbon is then added to the leaf litter flux :math:`l_{leaf,coh}`
(KgC individual\ :math:`^{-1}`)

.. math:: l_{leaf,coh} = b_{leaf,coh}

The alive biomass is depleted by the quantity of leaf mass lost, and the
leaf biomass is set to zero

.. math:: b_{alive,coh} = b_{alive,coh} - b_{leaf,coh}

.. math:: b_{leaf,coh} = 0

Finally, the status :math:`S_{phen,coh}` is set to 1, indicating that
the leaves have fallen off.

For bud burst, or leaf-on, the same occurs in reverse. If the leaves are
off (:math:`S_{phen,coh}`\ =1) and the phenological status triggers
budburst (:math:`S_{phen,grid}`\ =2) then the leaf mass is set the
maximum of the leaf memory and the available store

.. math:: b_{leaf,coh} =  \textrm{max}\left(l_{memory,coh}, b_{store,coh}\right.)

this amount of carbon is removed from the store

.. math:: b_{store,coh} = b_{store,coh}  - b_{leaf,coh}

and the new leaf biomass is added to the alive pool

.. math:: b_{alive,coh} = b_{alive,coh}  + b_{leaf,coh}

Lastly, the leaf memory variable is set to zero and the phenological
status of the cohort back to ‘2’. No parameters are currently required
for this carbon accounting scheme.

.. raw:: latex

   \bigskip

.. raw:: latex

   \captionof{table}{Parameters needed for phenology model.  }

+-----------------+-----------------+-----------------+-----------------+
| Parameter       | Parameter Name  | Units           | indexed by      |
| Symbol          |                 |                 |                 |
+=================+=================+=================+=================+
| :math:`n_{crit, | Threshold of    | none            | -               |
| cold}`          | cold days for   |                 |                 |
|                 | senescence      |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`T_{g}`   | Threshold for   | :math:`^{o}`\ C |                 |
|                 | counting        |                 |                 |
|                 | growing degree  |                 |                 |
|                 | days            |                 |                 |
+-----------------+-----------------+-----------------+-----------------+

.. raw:: latex

   \bigskip 

Seed Dynamics and Recruitment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


The production of seeds and their subsequent germination is a process
that must be captured explicitly or implicitly in vegetation models. FATES contains a seed bank model designed to allow the dynamics of
seed production and germination to be simulated independently. In the
ED1.0 model, seed recruitment occurs in the same timestep as allocation
to seeds, which prohibits the survival of a viable seed bank through a
period of disturbance or low productivity (winter, drought). In FATES, a plant functional type specific seed bank is tracked in
each patch (:math:`Seeds_{patch}` KgC m\ :math:`^{-2}`), whose rate of
change (KgC m\ :math:`^{-2}` y\ :math:`^{-1}`) is the balance of inputs,
germination and decay:

.. math:: \frac{\delta Seeds_{FT}}{\delta t } = Seed_{in,ft} - Seed_{germ,ft} - Seed_{decay,ft}

where :math:`Seed_{in}`, :math:`Seed_{germ}` and :math:`Seed_{decay}`
are the production, germination and decay (or onset of inviability) of
seeds, all in KgC m\ :math:`^{-2}` year\ :math:`^{-1}`.

Seeds are assumed to be distributed evenly across the site (in this
version of the model), so the total input to the seed pool is therefore
the sum of all of the reproductive output of all the cohorts in each
patch of the correct PFT type.

.. math:: Seed_{in,ft} =  \frac{\sum_{p=1}^{n_{patch}}\sum_{i=1}^{n_{coh}}p_{seed,i}.n_{coh}}{area_{site}}

Seed decay is the sum of all the processes that reduce the number of
seeds, taken from :ref:`Lischke et al. 2006<lischke2006>`. Firstly, the rate at
which seeds become inviable is described as a constant rate :math:`\phi`
(y:math:`^{-1}`) which is set to 0.51, the mean of the parameters used
by :ref:`Lischke et al. 2006<lischke2006>`.

.. math:: Seed_{decay,ft} = Seeds_{FT}.\phi

The seed germination flux is also prescribed as a fraction of the
existing pool (:math:`\alpha_{sgerm}`), but with a cap on maximum
germination rate :math:`\beta_{sgerm}`, to prevent excessive dominance
of one plant functional type over the seed pool.

.. math:: Seed_{germ,ft} = \textrm{max}(Seeds_{FT}\cdot \alpha_{sgerm},\beta_{sgerm})

.. raw:: latex

   \bigskip

.. raw:: latex

   \captionof{table}{Parameters needed for seed model.  }

+-----------------+-----------------+-----------------+-----------------+
| Parameter       | Parameter Name  | Units           | indexed by      |
| Symbol          |                 |                 |                 |
+=================+=================+=================+=================+
| :math:`K_s`     | Maximum seed    | kgC m\          |                 |
|                 | mass            | :math:`^{-2}`   |                 |
|                 |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`\alpha_{ | Proportional    | -               |                 |
| sgerm}`         | germination     |                 |                 |
|                 | rate            |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`\beta_{s | Maximum         | KgC             |                 |
| germ}`          | germination     | m\ :math:`^{-2}`|                 |
|                 | rate            |                 |                 |
|                 |                 | y\ :math:`^{-1}`|                 |
|                 |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`\phi`    | Decay rate of   | none            | FT              |
|                 | viable seeds    |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`R_{frac, | Fraction of     | none            | FT              |
| ft}`            | :math:`C_{bal}` |                 |                 |
|                 | devoted to      |                 |                 |
|                 | reproduction    |                 |                 |
+-----------------+-----------------+-----------------+-----------------+

.. raw:: latex

   \bigskip 

Litter Production and Fragmentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| The original CLM4.5 model contains streams of carbon pertaining to
  different chemical properties of litter (lignin, cellulose and labile
  streams, specifically). In FATES model, the fire simulation
  scheme in the SPITFIRE model requires that the model tracks the pools
  of litter pools that differ with respect to their propensity to burn
  (surface area-volume ratio, bulk density etc.). Therefore, this model
  contains more complexity in the representation of coarse woody debris.
  We also introduce the concept of ’fragmenting’ pools, which are pools
  that can be burned, but are not available for decomposition or
  respiration. In this way, we can both maintain above-ground pools that
  affect the rate of burning, and the lag between tree mortality and
  availability of woody material for decomposition.
| FATES recognizes four classes of litter. Above- and below-ground
  coarse woody debris (:math:`CWD_{AG}`, :math:`CWD_{BG}`) and leaf
  litter (:math:`l_{leaf}` and fine root litter :math:`l_{root}`). All
  pools are represented per patch, and with units of kGC
  m\ :math:`{^-2}`. Further to this, :math:`CWD_{AG}`, :math:`CWD_{BG}`
  are split into four litter size classes (:math:`lsc`) for the purposes
  of proscribing this to the SPITFIRE fire model (seed ’Fuel Load’
  section for more detail. 1-hour (twigs), 10-hour (small branches),
  100-hour (large branches) and 1000-hour(boles or trunks). 4.5 %, 7.5%,
  21 % and 67% of the woody biomass (:math:`b_{store,coh} + b_{sw,coh}`)
  is partitioned into each class, respectively.

:math:`l_{leaf}` and :math:`l_{root}` are indexed by plant functional
type (:math:`ft`). The rational for indexing leaf and fine root by PFT
is that leaf and fine root matter typically vary in their
carbon:nitrogen ratio, whereas woody pools typically do not.

Rates of change of litter, all in kGC m\ :math:`{^-2}`
year\ :math:`^{-1}`, are calculated as

.. math:: \frac{\delta CWD_{AG,out,lsc}}{ \delta t }= CWD_{AG,in,lsc} - CWD_{AG,out,lsc}

.. math:: \frac{\delta CWD_{BG,out,lsc}}{ \delta t } = CWD_{BG,in,lsc} - CWD_{BG,in,lsc}

.. math:: \frac{\delta l_{leaf,out,ft} }{ \delta t } = l_{leaf,in,ft} -  l_{leaf,out,ft}

.. math:: \frac{\delta l_{root,out,ft} }{ \delta t } = l_{root,in,ft} - l_{root,out,ft}

Litter Inputs
-------------

Inputs into the litter pools come from tissue turnover, mortality of
canopy trees, mortality of understorey trees, mortality of seeds, and
leaf senescence of deciduous plants.

.. math:: l_{leaf,in,ft} =\Big(\sum_{i=1}^{n_{coh,ft}} n_{coh}(l_{md,coh}  + l_{leaf,coh}) + M_{t,coh}.b_{leaf,coh}\Big)/\sum_{p=1}^{n_{pat}}A_{patch}

where :math:`l_{md,coh}` is the leaf turnover rate for evergreen trees
and :math:`l_{leaf,coh}` is the leaf loss from phenology in that
timestep (KgC :math:`m^{-2}`. :math:`M_{t,coh}` is the total mortality
flux in that timestep (in individuals). For fine root input.
:math:`n_{coh,ft}` is the number of cohorts of functional type
‘:math:`FT`’ in the current patch.

.. math:: l_{root,in,ft} =\Big(\sum_{i=1}^{n_{coh,ft}} n_{coh}(r_{md,coh} ) + M_{t,coh}.b_{root,coh}\Big)/\sum_{p=1}^{n_{pat}}A_p

where :math:`r_{md,coh}` is the root turnover rate. For coarse woody
debris input (:math:`\mathit{CWD}_{AG,in,lsc}` , we first calculate the
sum of the mortality :math:`M_{t,coh}.(b_{struc,coh}+b_{sw,coh})` and
turnover :math:`n_{coh}(w_{md,coh}`) fluxes, then separate these into
size classes and above/below ground fractions using the fixed fractions
assigned to each (:math:`f_{lsc}` and :math:`f_{ag}`)

.. math:: \mathit{CWD}_{AG,in,lsc} =\Big(f_{lsc}.f_{ag}\sum_{i=1}^{n_{coh,ft}}n_{coh}w_{md,coh}  + M_{t,coh}.(b_{struc,coh}+b_{sw,coh})\Big)/\sum_{p=1}^{n_{pat}}A_p

.. math:: \mathit{CWD}_{BG,in,lsc} =\Big(f_{lsc}.(1-f_{ag})\sum_{i=1}^{n_{coh,ft}}n_{coh}w_{md,coh}  + M_{t,coh}.(b_{struc,coh}+b_{sw,coh})\Big)/\sum_{p=1}^{n_{pat}}A_p

Litter Outputs
--------------

The fragmenting litter pool is available for burning but not for
respiration or decomposition. Fragmentation rates are calculated
according to a maximum fragmentation rate (:math:`\alpha_{cwd,lsc}` or
:math:`\alpha_{litter}`) which is ameliorated by a temperature and water
dependent scalar :math:`S_{tw}`. The form of the temperature scalar is
taken from the existing CLM4.5BGC decomposition cascade calculations).
The water scaler is equal to the water limitation on photosynthesis
(since the CLM4.5BGC water scaler pertains to the water potential of
individual soil layers, which it is difficult to meaningfully average,
given the non-linearities in the impact of soil moisture). The scaler
code is modular, and new functions may be implemented trivially. Rate
constants for the decay of the litter pools are extremely uncertain in
literature, as few studies either separate litter into size classes, nor
examine its decomposition under non-limiting moisture and temperature
conditions. Thus, these parameters should be considered as part of
sensitivity analyses of the model outputs.

.. math:: \mathit{CWD}_{AG,out,lsc} = CWD_{AG,lsc}. \alpha_{cwd,lsc}.S_{tw}

.. math:: \mathit{CWD}_{BG,out,lsc} = CWD_{BG,lsc} .\alpha_{cwd,lsc}.S_{tw}

.. math:: l_{leaf,out,ft} = l_{leaf,ft}.\alpha_{litter}.S_{tw}

.. math:: l_{root,out,ft} = l_{root,ft}.\alpha_{root,ft}.S_{tw}

Flux into decompsition cascade
------------------------------

Upon fragmentation and release from the litter pool, carbon is
transferred into the labile, lignin and cellulose decomposition pools.
These pools are vertically resolved in the biogeochemistry model. The
movement of carbon into each vertical layer is obviously different for
above- and below-ground fragmenting pools. For each layer :math:`z` and
chemical litter type :math:`i`, we derive a flux from ED into the
decomposition cascade as :math:`ED_{lit,i,z}` (kGC m\ :math:`^{-2}`
s\ :math:`^{-1}`)

where :math:`t_c` is the time conversion factor from years to seconds,
:math:`f_{lab,l}`, :math:`f_{cel,l}` and :math:`f_{lig,l}` are the
fractions of labile, cellulose and lignin in leaf litter, and
:math:`f_{lab,r}`, :math:`f_{cel,r}` and :math:`f_{lig,r}` are their
counterparts for root matter. Similarly, :math:`l_{prof}`,
:math:`r_{f,prof}`\ and :math:`r_{c,prof}` are the fractions of leaf,
coarse root and fine root matter that are passed into each vertical soil
layer :math:`z`, derived from the CLM(BGC) model.

.. raw:: latex

   \bigskip

.. raw:: latex

   \captionof{table}{Parameters needed for litter model.  }

+-----------------+-----------------+-----------------+-----------------+
| Parameter       | Parameter Name  | Units           | indexed by      |
| Symbol          |                 |                 |                 |
+=================+=================+=================+=================+
| :math:`\alpha_{ | Maximum         | y\ :math:`^{-1}`|                 |
| cwd,lsc}`       | fragmentation   |                 |                 |
|                 | rate of CWD     |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`\alpha_{ | Maximum         | y\ :math:`^{-1}`|                 |
| litter}`        | fragmentation   |                 |                 |
|                 | rate of leaf    |                 |                 |
|                 | litter          |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`\alpha_{ | Maximum         | y\ :math:`^{-1}`|                 |
| root}`          | fragmentation   |                 |                 |
|                 | rate of fine    |                 |                 |
|                 | root litter     |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`f_{lab,l | Fraction of     | none            |                 |
| }`              | leaf mass in    |                 |                 |
|                 | labile carbon   |                 |                 |
|                 | pool            |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`f_{cel,l | Fraction of     | none            |                 |
| }`              | leaf mass in    |                 |                 |
|                 | cellulose       |                 |                 |
|                 | carbon pool     |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`f_{lig,l | Fraction of     | none            |                 |
| }`              | leaf mass in    |                 |                 |
|                 | lignin carbon   |                 |                 |
|                 | pool            |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`f_{lab,r | Fraction of     | none            |                 |
| }`              | root mass in    |                 |                 |
|                 | labile carbon   |                 |                 |
|                 | pool            |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`f_{cel,r | Fraction of     | none            |                 |
| }`              | root mass in    |                 |                 |
|                 | cellulose       |                 |                 |
|                 | carbon pool     |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`f_{lig,r | Fraction of     | none            |                 |
| }`              | root mass in    |                 |                 |
|                 | lignin carbon   |                 |                 |
|                 | pool            |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`l_{prof, | Fraction of     | none            | soil layer      |
| z}`             | leaf matter     |                 |                 |
|                 | directed to     |                 |                 |
|                 | soil layer z    |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`r_{c,pro | Fraction of     | none            | soil layer      |
| f,z}`           | coarse root     |                 |                 |
|                 | matter directed |                 |                 |
|                 | to soil layer z |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`r_{f,pro | Fraction of     | none            | soil layer      |
| f,z}`           | fine root       |                 |                 |
|                 | matter directed |                 |                 |
|                 | to soil layer z |                 |                 |
+-----------------+-----------------+-----------------+-----------------+

.. raw:: latex

   \bigskip 

Plant Mortality
^^^^^^^^^^^^^^^^

Total plant mortality per cohort :math:`M_{t,coh}`, (fraction
year\ :math:`^{-1}`) is simulated as the sum of four additive terms,

.. math:: M_{t,coh}= M_{b,coh} + M_{cs,coh} + M_{hf,coh} + M_{f,coh},

where :math:`M_b` is the background mortality that is unaccounted by
any of the other mortality rates and is fixed at 0.014. :math:`M_{cs}`
is the carbon starvation derived mortality, which is a function of the
non-structural carbon storage term :math:`b_{store,coh}` and the
PFT-specific ‘target’ carbon storage, :math:`l_{targ,ft}`, as follows:

.. math:: M_{cs,coh}= \rm{max} \left(0.0, S_{m,ft} \left(0.5 -  \frac{b_{store,coh}}{l_{targ,ft}b_{leaf}}\right)\right)

where :math:`S_{m,ft}` is the `stress mortality` parameter, or the
fraction of trees in a landscape that die when the mean condition of a
given cohort triggers mortality. This parameter is needed to scale from
individual-level mortality simulation to grid-cell average conditions.

Mechanistic simulation of hydraulic failure is not undertaken on account
of it’s mechanistic complexity (see :ref:`McDowell et al. 2013<Mcdowelletal2013>` for
details). Instead, we use a proxy for hydraulic failure induced
mortality (:math:`M_{hf,coh}`) that uses a water potential threshold
beyond mortality is triggered, such that the tolerance of low water
potentials is a function of plant functional type (as expressed via the
:math:`\psi_c` parameter). For each day that the aggregate water
potential falls below a threshold value, a set fraction of the trees are
killed. The aggregation of soil moisture potential across the root zone
is expressed using the :math:`\beta` function. We thus determine plant
mortality caused by extremely low water potentials as

.. math::

   M_{hf,coh} = \left\{ \begin{array}{ll}
   S_{m,ft}& \textrm{for } \beta_{ft} < 10^{-6}\\
   &\\
   0.0& \textrm{for } \beta_{ft}>= 10^{-6}.\\
   \end{array} \right.

The threshold value of 10\ :math:`^{-6}` represents a state where the
average soil moisture potential is within 10\ :math:`^{-6}` of the
wilting point (a PFT specific parameter :math:`\theta_{w,ft}`).

:math:`M_{hf,coh}` is the fire-induced mortality, as described in the
fire modelling section.

.. raw:: latex

   \bigskip

.. raw:: latex

   \captionof{table}{Parameters needed for mortality model.  }

+---------------------+--------------------------------+-------+------------+
| Parameter Symbol    | Parameter Name                 | Units | indexed by |
+=====================+================================+=======+============+
| :math:`S_{m,ft}`    | Stress Mortality Scaler        | none  |            |
+---------------------+--------------------------------+-------+------------+
| :math:`l_{targ,ft}` | Target carbon storage fraction | none  | ft         |
+---------------------+--------------------------------+-------+------------+

.. raw:: latex

   \bigskip 

Fire (SPITFIRE)
^^^^^^^^^^^^^^^^^

The influence of fire on vegetation is estimated using the SPITFIRE
model, which has been modified for use in ED following it’s original
implementation in the LPJ-SPITFIRE model
(:ref:`Thonicke et al. 2010<thonickeetal2010>`, :ref:`Pfeiffer et al. 2013<pfeiffer2013>`). This model as
described is substantially different from the existing CLM4.5 fire model
:ref:`Li et al. 2012<Lietal2012a>`, however, further developments are
intended to increase the merging of SPITFIRE’s natural vegetation fire
scheme with the fire suppression, forest-clearing and peat fire
estimations in the existing model. The coupling to the ED model allows
fires to interact with vegetation in a size-structured manner, so
small fires can burn only understorey vegetation. Also, the patch
structure and representation of succession in the ED model allows the
model to track the impacts of fire on different forest stands, therefore
removing the problem of area-averaging implicit in area-based DGVMs. The
SPITFIRE approach has also been coupled to the LPJ-GUESS
individual-based model (Forrest et al. in prep) and so this is not the
only implementation of this type of scheme in existence.

The SPITFIRE model operates at a daily timestep and at the patch level,
meaning that different litter pools and vegetation charecteristics of
open and closed forests can be represented effectively (we omit the
`patch` subscript throughout for simplicity).

Properties of fuel load
-----------------------

Many fire processes are impacted by the properties of the litter pool in
the SPITFIRE model. There are one live (live grasses) and five dead fuel
categories (dead leaf litter and four pools of coarse woody debris).
Coarse woody debris is classified into 1h, 10h, 100h, and 1000h fuels,
defined by the order of magnitude of time required for fuel to lose
(or gain) 63% of the difference between its current moisture content and
the equilibrium moisture content under defined atmospheric conditions.
:ref:`Thonicke et al. 2010<thonickeetal2010>`. For the purposes of describing
the behaviour of fire, we introduce a new index 'fuel class' *fc*, the
values of which correspond to each of the six possible fuel categories
as follows.

+------------+------------------+-------------+
| *fc* index | Fuel type        | Drying Time |
+============+==================+=============+
| 1          | dead grass       | n/a         |
+------------+------------------+-------------+
| 2          | twigs            | 1h fuels    |
+------------+------------------+-------------+
| 3          | small branches   | 10h fuel    |
+------------+------------------+-------------+
| 4          | large branches   | 100h fuel   |
+------------+------------------+-------------+
| 5          | stems and trunks | 1000h fuel  |
+------------+------------------+-------------+
| 6          | live grasses     | n/a         |
+------------+------------------+-------------+

.. raw:: latex

   \bigskip 

Nesterov Index
---------------

Dead fuel moisture (:math:`\emph{moist}_{df,fc}`), and several other
properties of fire behaviour, are a function of the ‘Nesterov Index’
(:math:`N_{I}`) which is an accumulation over time of a function of
temperature and humidity (Eqn 5, :ref:`Thonicke et al. 2010<Thonickeetal2010>`),

.. math:: N_{I}=\sum{\textrm{max}(T_{d}(T_{d}-D),0)}

where :math:`T_{d}` is the daily mean temperature in :math:`^{o}`\ C and
:math:`D` is the dew point calculated as .

.. math::

   \begin{aligned}
   \upsilon&=&\frac{17.27T_{d}}{237.70+T_{d}}+\log(RH/100)\\
   D&=&\frac{237.70\upsilon}{17.27-\upsilon}\end{aligned}

where :math:`RH` is the relative humidity (%).

On days when the total precipitation exceeds 3.0mm, the Nesterov index accumulator is reset back to zero.

Fuel properties
---------------

Total fuel load :math:`F_{tot,patch}` for a given patch is the sum of
the above ground coarse woody debris and the leaf litter, plus the alive
grass leaf biomass :math:`b_{l,grass}` multiplied by the non-mineral
fraction (1-:math:`M_{f}`).

.. math:: F_{tot,patch}=\left(\sum_{fc=1}^{fc=5}  CWD_{AG,fc}+l_{litter}+b_{l,grass}\right)(1-M_{f})

Many of the model behaviours are affected by the patch-level weighted
average properties of the fuel load. Typically, these are calculated in
the absence of 1000-h fuels because these do not contribute greatly to
fire spread properties.


Dead Fuel Moisture Content
~~~~~~~~~~~~~~~~~~~~~~~~~~

Dead fuel moisture is calculated as

.. math:: \emph{moist}_{df,fc}=e^{-\alpha_{fmc,fc}N_{I}}

where :math:`\alpha_{fmc,fc}` is a parameter defining how fuel moisture
content varies between the first four dead fuel classes.

Live grass moisture Content
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The live grass fractional moisture content(\ :math:`\emph{moist}_{lg}`)
is a function of the soil moisture content. (Equation B2 in
:ref:`Thonicke et al. 2010<Thonickeetal2010>`)

.. math:: \emph{moist}_{lg}=\textrm{max}(0.0,\frac{10}{9}\theta_{30}-\frac{1}{9})

where :math:`\theta_{30}` is the fractional moisture content of the top
30cm of soil.

Patch Fuel Moisture
~~~~~~~~~~~~~~~~~~~

The total patch fuel moisture is based on the weighted average of the
different moisture contents associated with each of the different live
grass and dead fuel types available (except 1000-h fuels).

.. math:: F_{m,patch}=\sum_{fc=1}^{fc=4}  \frac{F_{fc}}{F_{tot}}\emph{moist}_{df,fc}+\frac{b_{l,grass}}{F_{tot}}\emph{moist}_{lg}

Effective Fuel Moisture Content
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Effective Fuel Moisture Content is used for calculations of fuel
consumed, and is a function of the ratio of dead fuel moisture content
:math:`M_{df,fc}` and the moisture of extinction factor,
:math:`m_{ef,fc}`

.. math:: E_{moist,fc}=\frac{\emph{moist}_{fc}}{m_{ef,fc}}

where the :math:`m_{ef}` is a function of surface-area to volume ratio.

.. math:: m_{ef,fc}=0.524-0.066\log_{10}{\sigma_{fc}}

Patch Fuel Moisture of Extinction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The patch ‘moisture of extinction’ factor (:math:`F_{mef}`) is the
weighted average of the :math:`m_{ef}` of the different fuel classes

.. math:: F_{mef,patch}=\sum_{fc=1}^{fc=5}  \frac{F_{fc}}{F_{tot}}m_{ef,fc}+\frac{b_{l,grass}}{F_{tot}}m_{ef,grass}

Patch Fuel Bulk Density
~~~~~~~~~~~~~~~~~~~~~~~

The patch fuel bulk density is the weighted average of the bulk density
of the different fuel classes (except 1000-h fuels).

.. math:: F_{bd,patch}=\sum_{fc=1}^{fc=4} \frac{F_{fc}}{F_{tot}}\beta_{fuel,fc}+\frac{b_{l,grass}}{F_{tot}}\beta_{fuel,lgrass}

where :math:`\beta_{fuel,fc}` is the bulk density of each fuel size
class (kG m\ :math:`^{-3}`)

Patch Fuel Surface Area to Volume
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The patch surface area to volume ratio (:math:`F_{\sigma}`) is the
weighted average of the surface area to volume ratios
(:math:`\sigma_{fuel}`) of the different fuel classes (except 1000-h
fuels).

.. math:: F_{\sigma}=\sum_{fc=1}^{fc=4}  \frac{F_{fc}}{F_{tot}}\sigma_{fuel,fc}+\frac{b_{l,grass}}{F_{tot}}\sigma_{fuel,grass}

Forward rate of spread
----------------------

For each patch and each day, we calculate the rate of forward spread of
the fire *ros*\ :math:`_{f}` (nominally in the direction of the wind).

.. math:: \emph{ros}_{f}=\frac{i_{r}x_{i}(1+\phi_{w})}{F_{bd,patch}e_{ps}q_{ig}}

:math:`e_{ps}` is the effective heating number
(:math:`e^{\frac{-4.528}{F_{\sigma,patch}}}`). :math:`q_{ig}` is the
heat of pre-ignition (:math:`581+2594F_{m}`). :math:`x_{i}` is the
propagating flux calculated as (see :ref:`Thonicke et al. 2010<Thonickeetal2010>`
Appendix A).

.. math::

   x_{i}= \frac{e^{0.792+3.7597F_{\sigma,patch}^{0.5}(\frac{F_{bd,patch}}{p_{d}}+0.1)}}{192+7.9095F_{\sigma,patch}}

:math:`\phi_{w}` is the influence of windspeed on rate of spread.

.. math:: \phi_{w}=cb_{w}^{b}.\beta^{-e}

Where :math:`b`, :math:`c` and :math:`e` are all functions of
surface-area-volume ratio :math:`F_{\sigma,patch}`:
:math:`b=0.15988F_{\sigma,patch}^{0.54}`,
:math:`c=7.47e^{-0.8711F_{\sigma,patch}^{0.55}}`,
:math:`e=0.715e^{-0.01094F_{\sigma,patch}}`. :math:`b_{w}=196.86W` where
:math:`W` is the the windspeed in ms\ :math:`^{-1}`, and
:math:`\beta=\frac{F_{bd}/p_{d}}{0.200395F_{\sigma,patch}^{-0.8189}}`
where :math:`p_{d}` is the particle density (513).

:math:`i_{r}` is the reaction intensity, calculated using the following
set of expressions (from :ref:`Thonicke et al. 2010<Thonickeetal2010>` Appendix A).:

.. math::

   \begin{aligned}
   i_{r}&=&\Gamma_{opt}F_{tot}Hd_{moist}d_{miner}\\
   d_{moist}&=&\textrm{max}\Big(0.0,(1-2.59m_{w}+5.11m_{w}^{2}-3.52m_{w}^{3})\Big)\\
   m_{w}&=&\frac{F_{m,patch}}{F_{mef,patch}}\\
   \Gamma _{opt}&=&\Gamma_{max}\beta^{a}\lambda\\
   \Gamma _{max}&=&\frac{1}{0.0591+2.926F_{\sigma,patch}^{-1.5}}\\
   \lambda&=&e^{a(1-\beta)}\\
   a&=&8.9033F_{\sigma,patch}^{-0.7913}\end{aligned}

:math:`\Gamma_{opt}` is the residence time of the fire, and
:math:`d_{miner}` is the mineral damping coefficient
(=0.174:math:`S_e^{-0.19}` , where :math:`S_e` is 0.01 and so =
:math:`d_{miner}` 0.41739).

Fuel Consumption
----------------

The fuel consumption (fraction of biomass pools) of each dead biomass
pool in the area affected by fire on a given day (:math:`f_{c,dead,fc}`)
is a function of effective fuel moisture :math:`E_{moist,fc}` and size
class *fc* (Eqn B1, B4 and B5, :ref:`Thonicke et al. 2010<Thonickeetal2010>`). The
fraction of each fuel class that is consumed decreases as its moisture
content relative to its moisture of extinction (:math:`E_{moist,fc}`)
increases.

.. math:: f_{cdead,fc}=\textrm{max}\left(0,\textrm{min}(1,m_{int,mc,fc}-m_{slope,mc,fc}E_{moist,fc})\Big)\right.

:math:`m_{int}` and :math:`m_{slope}` are parameters, the value of which
is modulated by both size class :math:`fc` and by the effective fuel
moisture class :math:`mc`, defined by :math:`E_{moist,fc}`.
:math:`m_{int}` and :math:`m_{slope}` are defined for low-, mid-, and
high-moisture conditions, the boundaries of which are also functions of
the litter size class following :ref:`Peterson and Ryan 1986 <Peterson1986>` (page
802). The fuel burned, :math:`f_{cground,fc}` (Kg m\ :math:`^{-2}`
day\ :math:`^{-1}`) iscalculated from :math:`f_{cdead,fc}` for each fuel
class:

.. math:: f_{cground,fc}=f_{c,dead,fc}(1-M_{f})\frac{F_{fc}}{0.45}

Where 0.45 converts from carbon to biomass. The total fuel consumption,
:math:`f_{ctot,patch}`\ (Kg m\ :math:`^{-2}`), used to calculate fire
intensity, is then given by

.. math:: f_{ctot,patch}=\sum_{fc=1}^{fc=4} f_{c,ground,fc} +  f_{c,ground,lgrass}

There is no contribution from the 1000 hour fuels to the patch-level
:math:`f_{ctot,patch}` used in the fire intensity calculation.

Fire Intensity
--------------

Fire intensity at the front of the burning area (:math:`I_{surface}`, kW
m\ :math:`^{-2}`) is a function of the total fuel consumed
(:math:`f_{ctot,patch}`) and the rate of spread at the front of the
fire, :math:`\mathit{ros}_{f}` (m min\ :math:`^{-1}`) (Eqn 15
:ref:`Thonicke et al. 2010<Thonickeetal2010>`)

.. math:: I_{surface}=\frac{0.001}{60}f_{energy} f_{ctot,patch}\mathit{ros}_{f}

where :math:`f_{energy}` is the energy content of fuel (Kj/Kg - the
same, 18000 Kj/Kg for all fuel classes). Fire intensity is used to define whether an
ignition is successful. If the fire intensity is greater than 50Kw/m
then the ignition is successful.

Fire Duration
-------------

Fire duration is a function of the fire danger index with a maximum
length of :math:`F_{dur,max}` (240 minutes in
:ref:`Thonicke et al. 2010<Thonickeetal2010>` Eqn 14, derived from Canadian Forest
Fire Behaviour Predictions Systems)

.. math:: D_{f}=\textrm{min}\Big(F_{dur,max},\frac{F_{dur,max}}{1+F_{dur,max}e^{-11.06fdi}}\Big)

Fire Danger Index
-----------------

Fire danger index (*fdi*) is a representation of the effect of
meteorological conditions on the likelihood of a fire. It is calculated
for each gridcell as a function of the Nesterov Index .
:math:`\emph{fdi}` is calculated from :math:`NI` as

.. math:: \emph{fdi}=1-e^{\alpha N_{I}}

where :math:`\alpha` = 0.00037 following
:ref:`Venevsky et al. 2002<venevsky2002>`.

Area Burned
-----------

Total area burnt is assumed to be in the shape of an ellipse, whose
major axis :math:`f_{length}` (m) is determined by the forward and
backward rates of spread (:math:`ros_{f}` and :math:`ros_{b}`
respectively).

.. math:: f_{length}=F_{d}(ros_{b}+ros_{f})

:math:`ros_{b}` is a function of :math:`ros_{f}` and windspeed (Eqn 10
:ref:`Thonicke et al. 2010<Thonickeetal2010>`)

.. math:: ros_{b}=ros_{f}e^{-0.012W}

The minor axis to major axis ratio :math:`l_{b}` of the ellipse is
determined by the windspeed. If the windspeed (:math:`W`) is less than
16.67 ms\ :math:`^{-1}` then :math:`l_{b}=1`. Otherwise (Eqn 12 and 13,
:ref:`Thonicke et al. 2010<Thonickeetal2010>`)

.. math:: l_{b}=\textrm{min}\Big(8,f_{tree}(1.0+8.729(1.0-e^{-0.108W})^{2.155})+(f_{grass}(1.1+3.6W^{0.0464}))\Big)

:math:`f_{grass}` and :math:`f_{tree}` are the fractions of the patch
surface covered by grass and trees respectively.

The total area burned (:math:`A_{burn}` in m\ :math:`^{2}`) is therefore
(Eqn 11, :ref:`Thonicke et al. 2010<Thonickeetal2010>`)

.. math:: A_{burn}=\frac{n_{f}\frac{3.1416}{4l_{b}}(f_{length}^{2}))}{10000}

where :math:`n_{f}` is the number of fires.

Crown Damage
-------------

:math:`c_{k}` is the fraction of the crown which is consumed by the
fire. This is calculated from scorch height :math:`H_{s}`, tree height
:math:`h` and the crown fraction parameter :math:`F_{crown}` (Eqn 17
:ref:`Thonicke et al. 2010<Thonickeetal2010>`):

.. math::

   c_{k} = \left\{ \begin{array}{ll}
   0 & \textrm{for $H_{s}<(h-hF_{crown})$}\\
   1-\frac{h-H_{s}}{h-F_{crown}}& \textrm{for $h>H_{s}>(h-hF_{crown})$}\\
   1 & \textrm{for $H_{s}>h$ }
   \end{array} \right.

The scorch height :math:`H_{s}` (m) is a function of the fire intensity,
following :ref:`Byram, 1959<byram1959>`, and is proportional to a plant
functional type specific parameter :math:`\alpha_{s,ft}` (Eqn 16
:ref:`Thonicke et al. 2010<Thonickeetal2010>`):

.. math:: H_{s}=\sum_{FT=1}^{NPFT}{\alpha_{s,p}\cdot f_{biomass,ft}} I_{surface}^{0.667}

where :math:`f_{biomass,ft}` is the fraction of the above-ground
biomass in each plant functional type.

Cambial Damage and Kill
-----------------------

The cambial kill is a function of the fuel consumed :math:`f_{c,tot}`,
the bark thickness :math:`t_{b}`, and :math:`\tau_{l}`, the duration of
cambial heating (minutes) (Eqn 8, :ref:`Peterson and Ryan 1986<peterson1986>`):

.. math:: \tau_{l}=\sum_{fc=1}^{fc=5}39.4F_{p,c}\frac{10000}{0.45}(1-(1-f_{c,dead,fc})^{0.5})

Bark thickness is a linear function of tree diameter :math:`dbh_{coh}`,
defined by PFT-specific parameters :math:`\beta_{1,bt}` and
:math:`\beta_{2,bt}` (Eqn 21 :ref:`Thonicke et al. 2010<Thonickeetal2010>`):

.. math:: t_{b,coh}=\beta_{1,bt,ft}+\beta_{2,bt,ft}dbh_{coh}

The critical time for cambial kill, :math:`\tau_{c}` (minutes) is given
as (Eqn 20 :ref:`Thonicke et al. 2010<Thonickeetal2010>`):

.. math:: \tau_{c}=2.9t_{b}^{2}

The mortality rate caused by cambial heating :math:`\tau_{pm}` of trees
within the area affected by fire is a function of the ratio between
:math:`\tau_{l}` and :math:`\tau_{c}` (Eqn 19,
:ref:`Thonicke et al. 2010<Thonickeetal2010>`):

.. math::

   \tau_{pm} = \left\{ \begin{array}{ll}
   1.0 & \textrm{for } \tau_{1}/\tau_{c}\geq \textrm{2.0}\\
   0.563(\tau_{l}/\tau_{c}))-0.125 & \textrm{for } \textrm{2.0} > \tau_{1}/\tau_{c}\ge \textrm{0.22}\\
   0.0 & \textrm{for } \tau_{1}/\tau_{c}< \textrm{0.22}\\
   \end{array} \right.

.. raw:: latex

   \bigskip

.. raw:: latex

   \captionof{table}{Parameters needed for fire model.  }

+-----------------+-----------------+-----------------+-----------------+
| Parameter       | Parameter Name  | Units           | indexed by      |
| Symbol          |                 |                 |                 |
+=================+=================+=================+=================+
| :math:`\beta_{1 | Intercept of    | mm              | *FT*            |
| ,bt}`           | bark thickness  |                 |                 |
|                 | function        |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`\beta_{2 | Slope of bark   | mm              | *FT*            |
| ,bt}`           | thickness       | cm\ :math:`^{-1 |                 |
|                 | function        | }`              |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`F_{crown | Ratio of crown  | none            | *FT*            |
| }`              | height to total |                 |                 |
|                 | height          |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`\alpha_{ | Fuel moisture   | :math:`{^o}`\ C | *fc*            |
| fmc}`           | parameter       | \ :math:`^{-2}` |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`\beta_{f | Fuel Bulk       | kG              | *fc*            |
| uel}`           | Density         | m\ :math:`^{-3}`|                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`\sigma_{ | Surface area to | cm              | *fc*            |
| fuel,fc}`       | volume ratio    | :math:`^{-1}`   |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`m_{int}` | Intercept of    | none            | :math:`fc`,     |
|                 | fuel burned     |                 | moisture class  |
+-----------------+-----------------+-----------------+-----------------+
| :math:`m_{slope | Slope of fuel   | none            | :math:`fc`,     |
| }`              | burned          |                 | moisture class  |
+-----------------+-----------------+-----------------+-----------------+
| :math:`M_f`     | Fuel Mineral    |                 |                 |
|                 | Fraction        |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`F_{dur,m | Maximum         | Minutes         |                 |
| ax}`            | Duration of     |                 |                 |
|                 | Fire            |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`f_{energ | Energy content  | kJ/kG           |                 |
| y}`             | of fuel         |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| :math:`\alpha_{ | Flame height    |                 | *FT*            |
| s}`             | parameter       |                 |                 |
+-----------------+-----------------+-----------------+-----------------+



