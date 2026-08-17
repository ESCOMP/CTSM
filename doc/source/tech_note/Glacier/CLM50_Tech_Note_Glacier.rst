.. _rst_Glaciers:

Glaciers
========

In CLM, glaciers are represented using glacier land units (Chapter :numref:`rst_Surface Characterization, Vertical Discretization, and Model Input Requirements`) that simulate snow, ice, surface energy balance, and runoff processes over permanently glaciated surfaces. 

CLM distinguishes between two major glacier categories:

#. Continental glaciers, i.e., the Greenland and Antarctic Ice Sheets.

#. Mountain (or alpine) glaciers, represented in the Randolph Glacier Inventory (including glaciers on the peripheries of the ice sheets).

Some glacier processes and model infrastructure are shared between these glacier types e.g., glacier land units. However, the fully developed glacier capabilities in CLM currently focus on continental ice sheets through coupling with the Community Ice Sheet Model (CISM; https://escomp.github.io/cism-docs/). Refer to section :numref:`Mountain Glaciers` for ongoing and future developments focusing on glaciers outside the ice sheets.


.. _Glacier regions:

Glacier regions and their behaviors
-----------------------------------

Within CLM, the world's glaciers and ice sheets are divided into three default glacier regions, each with distinct behaviors related to elevation classes, glacial meltwater treatment, and runoff generation. These default configurations are summarized in Table :numref:`Table Glacier region behaviors`.

.. _Table Glacier region behaviors:

.. table:: Glacier region behaviors

 +---------------+---------------+---------------+---------------+
 | Region        | Elevation     | Glacial melt  | Ice runoff    |
 |               | classes       |               |               |
 +===============+===============+===============+===============+
 | Greenland     | Virtual       | Replaced by   | Remains ice   |
 |               |               | ice           |               |
 +---------------+---------------+---------------+---------------+
 | Antarctica    | Multiple      | Replaced by   | Remains ice   |
 |               |               | ice           |               |
 +---------------+---------------+---------------+---------------+
 | Mountain      | Single        | Remains in    | Melted        |
 | glaciers      |               | place         |               |
 +---------------+---------------+---------------+---------------+


The glacier regions differ in three primary respects: 

#. Elevation class configuration

   a. Multiple elevation classes (section :numref:`Multiple elevation class scheme`)
   b. Multiple elevation classes plus virtual elevation classes 
   c. A single elevation class whose elevation matches the atmospheric grid-cell topography, such that no atmospheric downscaling is applied

#. Treatment of glacial melt water

   a. Replaced by ice: Meltwater runs off and is immediately replaced by ice, maintaining a permanently frozen glacier column. In the absence of a dynamic ice sheet model, this treatment implicitly assumes an unlimited ice reservoir available for melting, with additional adjustments applied to maintain mass and energy conservation. This behavior is discussed further in section :numref:`Computation of the surface mass balance`.
   b. Remains in place: Regions using this approach cannot compute surface mass balance (SMB), because negative SMB would be physically inconsistent in the presence of retained liquid water on top of the ice column. Although this treatment is less realistic physically, it avoids the persistent negative ice runoff required by the "replaced by ice" formulation to conserve mass and energy. This behavior is particularly useful for mountain glaciers, where atmospheric topographic smoothing can produce unrealistically warm conditions. In such cases, avoiding unrealistic negative runoff is often preferable to representing more realistic glacier physics. Section  :numref:`Mountain Glaciers` provides an overview of ongoing/future work to address these limitations with mountain glaciers.

#. Treatment of runoff from snow capping 

   a. Ice runoff from snow capping remains ice. This serves as a crude parameterization of iceberg calving, and is most appropriate for regions with substantial real-world calving.
   b. Ice runoff from snow capping is melted, generating a negative sensible heat flux, and then routed as liquid runoff. This behavior matches that of non-glacier land units and is more appropriate in regions with little iceberg calving. It can also help avoid unrealistic ocean cooling and runaway sea ice growth.

Further detail on snow capping is provided in section :numref:`Runoff from glaciers and snow-capped surfaces`. Note that these runoff treatments are irrelevant when using an evolving, two-way-coupled ice sheet model, because the snow capping flux is transferred directly to CISM rather than routed as runoff.


.. note::
    The combination of "Glacial melt = Replaced by ice" and "Ice runoff = Melted" produces strongly non-physical behavior. During glacier melt, the model generates negative ice runoff under the "Replaced by ice" treatment. Under the "Ice runoff = Melted" treatment, this negative ice runoff is converted into negative liquid runoff and a positive sensible heat flux. The resulting behavior produces zero net runoff but an artificial positive sensible heat flux associated with glacier melt. Because this behavior is physically unrealistic, CLM does not allow this combination of glacier-region settings.

.. note::
    Note that the CLM Greenland region extends only to the Greenland Ice Sheet boundary defined by CISM. As a result, SMB is not computed by default for grid cells that lie within the broader CISM domain but outside the Greenland Ice Sheet itself (i.e., the peripheral glaciers). Presently, these non-Greenland portions of the CISM domain are treated using the same configuration as the mountain glacier regions, rather than using the Greenland ice-sheet configuration. This choice helps avoid unrealistic runoff fluxes from the Canadian Arctic Archipelago that could otherwise contribute to excessive sea ice growth in the surrounding ocean.

.. note::
    Non-virtual, non-SMB-computing glacier regions can exist within the CISM domain, as is the case for portions of the Greenland CISM domain outside the Greenland Ice Sheet itself. However, these regions always provide zero SMB and cannot respond to CISM-driven changes in glacier extent. For this reason, it is generally preferable for as much of the CISM domain as possible to use virtual, SMB-computing glacier regions.



.. _Ice_Sheets:

Ice Sheets
----------

CLM computes and provides two quantities that are passed to the ice sheet model:

#. Surface mass balance (SMB) - the net annual accumulation and ablation of mass at the upper surface (section :numref:`Computation of the surface mass balance`)

#. Ground surface temperature, which serves as an upper boundary condition for CISM's temperature calculation. Ice sheet models are typically run at much higher spatial resolution than CLM (for example, :math:`\sim\ 5km` versus :math:`\sim\ 100km`). To improve the downscaling of atmospheric forcing from the CLM grid to the ice sheet grid, the glaciated portion of each CLM grid cell is divided into multiple elevation classes (section :numref:`Multiple elevation class scheme`). The CESM coupler then performs horizontal and vertical interpolation to generate high-resolution fields for CISM.
 
**Static ice sheet configuration**

In typical simulations, CISM is run in a non-evolving configuration. In this mode, CLM computes SMB and passes it to CISM, but the ice-sheet geometry remains fixed throughout the simulation. Under this configuration, CISM serves two primary roles:


#. Defining glacier extent and topography: Within the CISM domain (typically Greenland in CESM2), CISM specifies glacier area and topographic elevation, overriding the corresponding values in the CLM surface dataset. CISM also defines the elevation of non-glacier land units within its domain. Atmospheric downscaling over non-glacier land units is applied only within the CISM domain. If the stub glacier model configuration (`SGLC <https://escomp.github.io/cism-docs/cism-in-cesm/versions/master/html/clm-cism-coupling.html#stub-glc-model-cism-absent>`_) is used instead of CISM, glacier areas and elevations are taken entirely from the CLM surface dataset, and no atmospheric downscaling is applied over non-glacier land units. 


#. Providing the target grid for SMB downscaling: CISM provides the high-resolution grid onto which SMB fields are downscaled. When using SGLC, SMB is still computed in CLM, but it is not interpolated to a separate high-resolution ice sheet grid.


**Evolving ice sheet configuration**

CESM can also be run with an evolving ice sheet. In this mode, CLM responds dynamically to changes in ice-sheet geometry computed by CISM. As the ice sheet evolves, CLM updates the glacier land-unit area, elevation-class area fractions, and mean elevation of each elevation class. This ensures that glacier extent and surface topography remain consistent between CLM and CISM throughout the simulation. Conservation of mass and energy follows the same framework used for other land-cover transitions (Chapter :numref:`rst_Transient Landcover Change`).


.. _Multiple elevation class scheme:

Multiple elevation class scheme
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The glacier land unit contains multiple columns based on surface elevation. These are known as elevation classes, and the land unit is referred to as *glacier\_mec*. As described in section :numref:`Glacier regions`, some regions have only a single elevation class, but they are still referred to as *glacier\_mec* land units. The default is to have 10 elevation classes whose lower limits are 0, 200, 400, 700, 1000, 1300, 1600, 2000, 2500, and 3000 m. Each column is characterized by a fractional area and surface elevation that are read in during model initialization, and then possibly overridden by CISM as the run progresses. Each *glacier\_mec* column within a grid cell has distinct ice and snow temperatures, snow water content, surface fluxes, and SMB.

The atmospheric surface temperature, potential temperature, specific humidity, density, and pressure are downscaled from the atmosphere's mean grid cell elevation to the *glacier\_mec* column elevation using a specified lapse rate (typically 6.0 deg/km) and an assumption of uniform relative humidity. Longwave radiation is downscaled by assuming a linear decrease in downwelling longwave radiation with increasing elevation (0.032 W m\ :sup:`-2` m\ :sup:`-1`, limited to 0.5 - 1.5 times the gridcell mean value, then normalized to conserve gridcell total energy) :ref:`(Van Tricht et al., 2016)<VanTrichtetal2016>`. Total precipitation is partitioned into rain vs. snow as described in Chapter :numref:`rst_Surface Characterization, Vertical Discretization, and Model Input Requirements`. The partitioning of precipitation is based on the downscaled temperature, allowing rain to fall at lower elevations while snow falls at higher elevations.

This downscaling allows lower-elevation columns to undergo surface melting while columns at higher elevations remain frozen. This gives a more accurate simulation of summer melting, which is a highly nonlinear function of air temperature. Within the CISM domain, this same downscaling procedure is also applied to all non-urban land units. The elevation of non-glacier land units is taken from the mean elevation of ice-free grid cells in CISM. This is done in order to keep the glaciated and non-glaciated portions of the CISM domain as consistent as possible.

In contrast to most CLM subgrid units, *glacier\_mec* columns can be active (i.e., have model calculations run there) even if their area is zero. These are known as "virtual" columns. This is done because the ice sheet model may require a SMB for some grid cells where CLM has zero glacier area in that elevation range. Virtual columns also facilitate glacial advance and retreat in the two-way coupled case. Virtual columns do not affect energy exchange between the land and the atmosphere.

.. _Computation of the surface mass balance:

Surface mass balance
~~~~~~~~~~~~~~~~~~~~~

Computing SMB in CLM rather than directly within CISM provides several advantages:

#. Computational efficiency: SMB can be computed much more efficiently in CLM using a limited number of elevation classes than on the full high-resolution CISM grid. For example, consider a simulation with CLM at 50 km resolution and CISM at 5 km resolution. The Greenland Ice Sheet spans roughly 1000 x 2000 km, corresponding to approximately:

    - 20 x 40 x 10 (elevation classes) = 8,000 glacier elevation-class columns in CLM
    - 200 x 400 = 80,000 grid cells in CISM

#. Shared snow and surface physics: CLM already contains a sophisticated snow and surface energy balance parameterization. Computing SMB within CLM avoids the need to implement and maintain a separate SMB scheme within CISM, while ensuring that improvements to CLM physics are automatically applied to ice sheets.

#. Atmosphere - ice sheet feedbacks: Computing SMB in CLM allows the atmosphere model to respond interactively to changes in ice-sheet surface properties, even without fully evolving ice-sheet geometry. As shown by :ref:`Pritchard et al. (2008)<Pritchardetal2008>`, interactive albedo feedbacks are critical for realistic simulations of long-term ice-sheet retreat.

#. Consistency across glacier types: Improvements to SMB calculations in CLM are potentially applicable to all glaciated grid cells, including mountain glaciers, and not only to continental ice sheets.


**Computation of the surface mass balance**

This section describes the computation of SMB and associated runoff terms. The discussion applies only to glacier regions where meltwater runs off and the lost ice is immediately replaced, rather than to regions where meltwater remains within the glacier column. By default, this treatment applies to the Greenland and Antarctic Ice Sheets, but not to mountain glaciers. 

The SMB of a glacier or ice sheet is defined as the net annual mass gain or loss at the upper surface. Accumulation occurs primarily through snowfall and deposition, while ablation occurs primarily through melting and evaporation/sublimation. Ablation is defined here as the mass of water that ultimately runs off to the ocean. Not all surface meltwater contributes directly to runoff; some meltwater percolates into the snowpack and refreezes.

CLM computes SMB using a surface energy balance (SEB) approach, in which melt depends on the combined radiative, turbulent, and conductive energy fluxes at the surface. In glaciology, SMB is typically defined as the net balance of both snow and ice accumulation and ablation. However, the SMB flux passed from CLM to CISM represents the mass balance of the underlying ice only, excluding transient changes in snow storage. Conceptually, CLM can be viewed as owning the snowpack, while CISM owns the underlying glacier ice. As a result, fluctuations in snow depth between 0 and 10 m water equivalent are not reflected in the SMB passed to CISM. In transient simulations, this treatment can delay the onset of accumulation or ablation signals in a glacier column by several decades.

SMB is computed and sent to the CESM coupler regardless of whether and where CISM is operating. However, the effect of SMB terms on runoff fluxes differs depending on whether and where CISM is evolving in two-way-coupled mode. This is described by the variable *glc\_dyn\_runoff\_routing* (this is real-valued in the code to handle the edge case where a CLM grid cell partially overlaps with the CISM grid, but we describe it as a logical variable here for simplicity.) In typical cases where CISM is not evolving, *glc\_dyn\_runoff\_routing* will be false everywhere; in these cases, CISM's mass is not considered to be part of the coupled system. In cases where CISM is evolving and sending its own calving flux to the coupler, *glc\_dyn\_runoff\_routing* will be true over the CISM domain and false elsewhere.

Any snow capping (section :numref:`Runoff from glaciers and snow-capped surfaces`) is added to :math:`q_{ice,frz}`. Any liquid water (i.e., melted ice) below the snow pack in the glacier column is added to :math:`q_{ice,melt}`, then is converted back to ice to maintain a pure-ice column. Then the total SMB is given by :math:`q_{ice,tot}`:

.. math::
   :label: 13.1

   q_{ice,tot} = q_{ice,frz} - q_{ice,melt}

CLM is responsible for generating glacial surface melt, even when running with an evolving ice sheet. Thus, :math:`q_{ice,melt}` is always added to liquid runoff (:math:`q_{rgwl}`), regardless of *glc\_dyn\_runoff\_routing*. However, the ice runoff flux depends on *glc\_dyn\_runoff\_routing*. If *glc\_dyn\_runoff\_routing* is true, then CISM controls the fate of the snow capping mass in :math:`q_{ice,frz}` (e.g., eventually transporting it to lower elevations where it can be melted or calved). Since CISM will now own this mass, the snow capping flux does *not* contribute to any runoff fluxes generated by CLM in this case.

If *glc\_dyn\_runoff\_routing* is false, then CLM sends the snow capping flux as runoff, as a crude representation of ice calving. However, this ice runoff flux is reduced by :math:`q_{ice,melt}`. This reduction is needed for conservation; its need is subtle, but can be understood with either of these explanations:

- When ice melts, we let the liquid run off and replace it with new ice. That new ice needs to come from somewhere to keep the coupled system in water balance. We "request" the new ice from the ocean by generating a negative ice runoff equivalent to the amount we have melted.

- Ice melt removes mass from the system, as it should. But the snow capping flux also removes mass from the system. The latter is a crude parameterization of calving, assuming steady state - i.e., all ice gain is balanced by ice loss. This removal of mass due to both accumulation and melt represents a double-counting. Each unit of melt indicates that one unit of accumulation should not have made it to the ocean as ice, but instead melted before it got there. So we need to correct for this double-counting by removing one unit of ice runoff for each unit of melt.

For a given point in space or time, this reduction can result in negative ice runoff. However, when integrated over space and time, for an ice sheet that is near equilibrium, this just serves to decrease the too-high positive ice runoff from snow capping. (The treatment of snow capping with *glc\_dyn\_runoff\_routing* false is based on this near-equilibrium assumption - i.e., that ice accumulation is roughly balanced by :math:`calving + melt`, integrated across space and time. For glaciers and ice sheets that violate this assumption, either because they are far out of equilibrium with the climate or because the model is being run for hundreds of years, there are two ways to avoid the unrealistic ice runoff from snow capping: by running with an evolving, two-way-coupled ice sheet or by changing a glacier region's ice runoff behavior as described in section :numref:`Glacier regions`.)

In regions where SMB is computed for glaciers, SMB is also computed for the natural vegetated land unit. Because there is no ice to melt in this land unit, it can only generate a zero or positive SMB. A positive SMB is generated once the snow pack reaches its maximum depth. When running with an evolving ice sheet, this condition triggers glacial inception.

.. _Mountain Glaciers:

Mountain Glaciers
-----------------

Beginning with CLM5.2, mountain glacier outlines are derived from version 6 of the Randolph Glacier Inventory (RGI).

As discussed earlier in this chapter, the current representation of mountain glaciers in CLM remains relatively limited compared to the treatment of continental ice sheets. By default, mountain glaciers are treated using a simplified single-elevation-class configuration, and SMB is not computed for these regions. This limitation arises from the glacier-region configuration described in section :numref:`Glacier regions`, where meltwater is retained within the glacier column rather than running off and being replaced by ice.

This treatment is primarily intended to avoid unrealistic runoff fluxes in regions where coarse atmospheric topography produces climates that are too warm to sustain glaciers realistically. Such issues are particularly important for mountain glaciers, where strong local elevation gradients are poorly resolved at typical climate-model resolutions.

Ongoing development efforts are focused on extending CLM's mountain glacier capabilities. In particular, the hillslope hydrology configuration (`Chapter Hillslope_Hydrology <https://escomp.github.io/CTSM/tech_note/Hillslope_Hydrology/CLM50_Tech_Note_Hillslope_Hydrology.html>`_) is being adapted to support SMB calculations for mountain glaciers . These developments are intended to enable coupling between CLM and CISM for mountain glacier applications in a manner similar to the current treatment of ice sheets, leveraging recent CISM capabilities for simulating mountain glacier dynamics :ref:`(Minallah and Lipscomb et al., 2025) <Minallah2025>`.
