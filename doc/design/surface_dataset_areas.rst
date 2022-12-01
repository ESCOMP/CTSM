.. sectnum::

.. contents::

==================================
 Overview of this design document
==================================

This document gives a high-level overview of the specification of sub-grid areas on surface datasets and the raw datasets used to create them.

See also https://github.com/ESCOMP/CTSM/issues/1716 for related discussion.

=================================================
 Specification of landunit areas in raw datasets
=================================================

In the high-resolution raw datasets used as input to the mksurfdata tool, landunit areas should be specified as percent of the grid cell, **not** percent of the land area. For example, on the urban raw dataset, if there is a grid cell that is 50% land and 50% ocean, and 60% of the land area is urban, the urban area in that grid cell should be given as 30%.

One reason for this is that it makes reconciling the different landunit areas more straightforward if different raw datasets disagree about the land mask. In addition, this convention makes it easier to map raw datasets to the model resolution. For example, consider averaging two grid cells into a destination grid cell: one with 2% land area, of which 50% is glacier; and one with 100% land area, of which none is glacier. If we encode these as percent of the land area, we would have 50% glacier and 0% glacier, and then the final glacier cover would be 25%, suggesting that 25% of the land area is glacier, but this is incorrect. If we instead encode these as percent of the total grid cell area, we would have 1% glacier and 0% glacier, and then the final glacier cover would be 0.5%, suggesting that 0.5% of the total grid cell is glacier, which is correct.

=====================================================
 Specification of landunit areas in surface datasets
=====================================================

In contrast to the raw datasets, landunit areas in surface datasets are specified as percent of the land area, **not** as percent of the total grid cell. This is because we don't know what the model's actual land fraction will be at the time when the surface datasets are created: instead, this land fraction is determined at runtime based on the ocean grid.

===========================================================================================
 Procedure for converting landunit areas from percent of grid cell to percent of land area
===========================================================================================

There are a few important aspects to how we determine final landunit areas in mksurfdata:

When mapping landunit areas from the raw data resolution to the model resolution, we initially want to maintain areas as percent of the total grid cell area. To achieve this, we set ``norm_by_fracs=.false.`` in the call to ``create_routehandle``, resulting in the use of ``ESMF_NORMTYPE_DSTAREA`` rather than ``ESMF_NORMTYPE_FRACAREA`` as a ``normType`` when computing mapping weights. Using ``FRACAREA`` normalization is appropriate when you want to treat areas outside of the source mask as missing values that do not contribute in any way to the final destination value. Using ``DSTAREA`` normalization, in contrast, essentially treats areas outside of the source mask as zeroes. ``FRACAREA`` normalization is appropriate for many surface dataset fields, but ``DSTAREA`` is appropriate for areas specified as percent of the grid cell. For example, if a source grid cell is entirely ocean, then we want to treat glacier area in that source grid cell as 0%.

The conversion from percent of the grid cell area to percent of the land area happens in the subroutine ``normalize_and_check_landuse``. An important piece of doing this conversion is to determine an estimate of the land fraction for each model grid cell. This is not straightforward given the disparate land masks used for each raw dataset. We start by using the land fraction from the vegetation (PFT) raw dataset, with the assumption that that is probably the most reliable land mask. However, there are areas where using that land fraction is problematic, particularly where the areas of other landunits extend beyond the PFT's land mask. This is especially an issue for glaciers, where floating ice shelves can extend beyond the land-ocean border. To deal with this problem, if the sum of the areas of special landunits and crops exceeds the land fraction from the PFT data, we instead use that sum as an estimate of the land fraction. Exactly which landunits to include in this sum is a bit arbitrary. Arguably, the natural vegetated landunit should also be included in this sum. However, we ideally want to avoid changes in this estimated land fraction through time in transient datasets, which means that we generally want to use the PFT data's land fraction, only falling back on this landunit sum in exceptional cases. By excluding the natural vegetated area from this sum, we are more likely to use the PFT's land fraction. An implication of this choice is that we are more likely to replace natural vegetated areas with special landunits in cases where there are disagreements between the different raw datasets in coastal grid cells. For more detailed explanation and rationale, see some comments in ``normalize_and_check_landuse``.

In grid cells where the estimated land fraction is essentially zero, we set the land cover to wetland, as a rough parameterization of ocean. This situation will only arise if the areas of all landunits on the raw datasets are essentially zero for the given grid cell, so we would have no information to choose any particular land cover for the grid cell. (This wetland area may end up being changed to bare ground at runtime, depending on the value of the ``convert_ocean_to_land`` namelist flag.)

In grid cells where the estimated land fraction is greater than zero, we fill any unclaimed land area with the natural vegetated landunit. We then normalize all landunit areas based on the estimated land fraction so that they now specify areas as percent of the land area rather than as percent of the grid cell.

===================
 Example scenarios
===================

The following example scenarios illustrate the operation of ``normalize_and_check_landuse``; in the following, any landunit not explicitly mentioned has 0% area:

(a) With pctlnd_pft = 0% and all initial landunit areas 0%: wetland area = 100%

(b) With pctlnd_pft = 0% and initial glacier area 1%: glacier area = 100%

(c) With pctlnd_pft > 0% and all initial landunit areas 0%: natural vegetated area = 100%

(d) With pctlnd_pft = 40%, initial crop area 20%, natural vegetated area 10%: crop area = 50%, natural vegetated area = 50%

(e) With pctlnd_pft = 40%, initial crop area 20%, natural vegetated area 10%, glacier area 10%: crop area = 50%, natural vegetated area = 25%, glacier area = 25%

(f) With pctlnd_pft = 40%, initial crop area 20%, natural vegetated area 10%, glacier area 15%: crop area = 50%, natural vegetated area = 12.5%, glacier area = 37.5%

(g) With pctlnd_pft = 40%, initial crop area 20%, natural vegetated area 10%, glacier area 20%: crop area = 50%, glacier area = 50%

(h) With pctlnd_pft = 40%, initial crop area 20%, natural vegetated area 10%, glacier area 30%: crop area = 40%, glacier area = 60%

(i) With pctlnd_pft = 40%, initial crop area 0%, natural vegetated area 40%, glacier area 40%: glacier area = 100%

(j) With pctlnd_pft = 2%, initial natural vegetated area 1%, glacier area 1%: natural vegetated area = 50%, glacier area = 50%

(k) With pctlnd_pft = 2%, initial natural vegetated area 0%, glacier area 1%: natural vegetated area = 50%, glacier area = 50%

(l) With pctlnd_pft = 2%, initial natural vegetated area 2%, glacier area 1%: natural vegetated area = 50%, glacier area = 50%
