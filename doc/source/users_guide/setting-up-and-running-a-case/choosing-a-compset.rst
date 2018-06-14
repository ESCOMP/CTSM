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


Here is the entire list of compsets available. 
Note that using the "-user_compset" option even more combinations are possible. 
In the list below we give the alias name and then the long name which describes each component in parenthesis. 
Alias (Long-name with time-period and each component)

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
