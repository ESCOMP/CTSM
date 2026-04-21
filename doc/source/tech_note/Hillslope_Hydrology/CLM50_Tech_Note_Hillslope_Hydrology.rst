.. _rst_Hillslope_Hydrology:

Hillslope Hydrology
=====================

CLM is typically operated using a single column to represent the hydrologic state of the vegetated land
unit.  The CLM 'Hillslope Hydrology' configuration extends the standard configuration by instantiating multiple columns per gridcell, each performing the standard vertical model physics :ref:`Swenson et al. 2019b<Swensonetal2019b>`.  A set of connected columns in a gridcell is referred to as a hillslope, and lateral runoff fluxes can be passed between hillslope columns.  In addition to calculating lateral moisture fluxes, meteorological downscaling is enabled in the hillslope configuration (:numref:`Figure Hillslope Hydrology Schematic`).  These processes require a description of the columns' topographic characteristics, which are described by specifying six geomorphic parameters: the column area :math:`A`, the mean column height :math:`h`, the width of the column at its downslope interface :math:`w`, the mean distance of the column from the channel :math:`d`, the mean slope of the column :math:`\alpha`, and the column's aspect (defined with respect to North) :math:`\beta` (:numref:`Figure Hillslope Hydrology Geomorphic Parameters`).  A method for estimating the geomorphic parameters from a digital elevation model (DEM) is described by :ref:`Swenson and Lawrence 2025<SwensonLawrence2025>`.

.. _Figure Hillslope Hydrology Schematic:

.. Figure:: image1.png

   Processes enabled in Hillslope Hydrology  
	    
.. _Figure Hillslope Hydrology Geomorphic Parameters:

.. Figure:: image2.png

   Hillslope geomorphic parameters
	    
Saturated subsurface flow along the hillslope profile is described by Darcy's law

.. math::
   :label: 999.1

   q = - K_{s} \frac{d\Psi}{dx} ,

where :math:`q` is the moisture flux (m\ :sup:`3` m\ :sup:`-2` s\ :sup:`-1`), :math:`K_{s}` is the saturated hydraulic conductivity (m/s), :math:`\Psi` is the hydraulic head (m), and :math:`x` is the distance from the base of the hillslope (m).  In unconfined saturated flow, the hydraulic head is the absolute water table height, :math:`h`.  To calculate the hydraulic head gradient, each column's water table depth must be subtracted from the relative hillslope elevation.

Downscaling of incoming solar radiation uses topographic slope and aspect to redistribute insolation between the gridcell's hillslope columns, while conserving the total energy flux.  Incoming direct solar radiation, :math:`S_{gridcell}`, is divided by the gridcell cosine of the solar zenith angle, :math:`cosz_{gridcell}`, then multiplied by the column's cosine of the solar zenith angle, :math:`cosz_{column}`, which is calculated using the column's topographic slope and aspect

.. math::
   :label: 999.2

   S_{column} = S_{gridcell} \frac{cosz_{column}}{cosz_{gridcell}} ,

Other meteorological inputs, e.g. temperature and precipitation, can be adjusted using elevation‐based downscaling.  For such variables, the gridcell-level value is modified by an anomly defined by the lapse rate :math:`\lambda`

.. math::
   :label: 999.3

   X_{column} = X_{gridcell} + \lambda * \Delta Z_{column}

where :math:`X` represents a meteorological variable and :math:`\Delta Z_{column}` represents the column elevation anomaly with respect to the gridcell mean elevation.

