.. sectnum::

.. contents::

==================================
 Overview of this design document
==================================

This documents some of the high-level design decisions in the implementation of water
tracers, which are used for water isotopes and potentially other tracers.

=============================================================
 Distinction between bulk-only vs. bulk-and-tracer variables
=============================================================

Options considered for how to set up the water data structures
==============================================================

Some water variables need to have separate instances for each water tracer, whereas some
apply only to bulk water. An example of the latter is snow depth.

Broadly speaking, we considered three ways to store information for each tracer:

1. Single instance of ``waterstate_type`` (etc.), with an extra dimension on arrays for
   which we need values for each tracer. Bulk water is the first index in these arrays.

2. Like (1), single instance of ``waterstate_type`` (etc.), with an extra dimension on
   arrays for which we need values for each tracer. Bulk water is a separate variable, so
   we have things like ``h2osoi_liq_col(c,j)`` for bulk, and
   ``h2osoi_liq_col_tracer(c,j,m)`` for tracers).

3. Multiple instances of ``waterstate_type`` (etc.), with no extra dimension required for
   individual variables. 

We decided to go with option (3) for a number of reasons:

1. This is consistent with what is done for carbon isotopes. It has felt like this design
   has worked well for carbon isotopes, and we felt it made sense to keep these two
   consistent unless there was a good reason to treat them differently.

2. With the object-oriented design laid out in `Object-oriented design for water tracer
   types`_, it is easy to change a variable from bulk-only to bulk-and-tracer, or vice
   versa.

3. In contrast to (1), science routines that operate on bulk water can be written without
   needing to worry about the extra dimension. Similarly, unit tests remain relatively
   easier to set up.

4. In contrast to (2), tracers get their history and restart variables added
   automatically, avoiding duplication and the potential of getting out-of-sync with the
   bulk.

5. We can explicitly pass the bulk instance of (e.g.) ``waterstate_type`` to science
   routines that expect to be operating on this bulk instance, making the intent of code
   more clear.

6. Compared with (1), it's easier to search the code for where tracers are handled.

There were some downsides to this option, though, both short-term and long-term:

1. This required more up-front work to set up the object-oriented infrastructure

2. As discussed in more depth in `All water-related variables should go in one of the
   water* types`_, this required moving all water-related fluxes into one of the main
   water types, rather than having water variables distributed amongst science
   modules. However, even with one of the other possible designs, we might have wanted to
   move these variables, because with any of these designs: Any derived type with water
   variables needs to have extra infrastructure for dealing with water tracers. So if we
   kept water variables distributed amongst various types, this tracer-handling would also
   be distributed in many places.

Object-oriented design for water tracer types
=============================================

The object-oriented design for water tracer types was fleshed out in discussions here:
https://github.com/ESCOMP/ctsm/pull/395.

One item of particular note is that we use inheritance to define the full list of
variables used for bulk water: The set of bulk water variables is the set of
bulk-and-tracer variables plus some additional variables that just apply to bulk water; we
avoid duplication by having bulk types that extend the basic types (e.g.,
``waterstatebulk_type`` extends ``waterstate_type``). Some benefits of this
inheritance-based design are:

1. Having a separate derived type for bulk water allows science routines to declare that
   they expect an argument of the bulk water type.

2. This allows us to easily migrate variables back and forth between being bulk-only or
   being bulk-and-tracer: You can simply move the variable (along with its associated
   allocations, history addfld calls, etc.) between (for example) ``waterstatebulk_type``
   and ``waterstate_type`` without needing to change all of the science modules that
   reference that variable (because they generally reference it through the bulk
   instance).

All water-related variables should go in one of the water* types
================================================================

An implication of this design (discussed in greater detail here:
https://github.com/ESCOMP/ctsm/pull/395#issuecomment-392299340) is that all water
variables should be defined in one of the centralized water* types. This is in contrast to
our general desire to have variables live in the modules that compute them. This is needed
to avoid having to split up a bunch of science modules into two types: a type that needs
multiple instances (for each water tracer) and a type that only needs a single instance
(for bulk water, or variables that are not water-related at all).

We would still like to keep some variables local to the science modules – especially
variables that can be made private to those modules. But the fluxes were never handled in
a completely (textbook) object-oriented fashion anyway, since they were accessed publicly
by other modules, so we're not too sad to move these public flux variables back out to the
generic flux types. And while this has the disadvantage of losing the attachment of a
variable to the module that computes that variable, it has some advantages, too, in terms
of understanding the suite of water flux variables, and reducing the number of instances
that need to get passed all around (e.g., we won't need to pass irrigation_inst to as many
places).

========================
 Loops over all tracers
========================

Initially, I was hoping that we could keep loops over tracers in WaterType, for the
following reasons:

1. To keep this complexity out of other modules

2. To make it easier to change the details of how the bulk and tracer instances are
   stored, if we ever need to: By keeping as many of the loops as possible in WaterType,
   we reduce the number of places that would need to be changed

However, it was starting to get too awkward to require all loops over tracers to happen in
WaterType (or some other centralized location): I had originally imagined that we wouldn't
need too many loops over tracers, but it turns out that we need loops over tracers in a
lot of places. Requiring all of these loops over tracers to be in WaterType both (a)
bloats that module, and (b) adds extra indirection (which makes it harder to understand
the code, because you're bouncing back and forth between more modules, and has possible
performance implications as we break routines into tiny pieces for this purpose).

So we allow loops over tracers (or bulk plus tracers) anywhere in the code. See comments
at the top of WaterType.F90 for example code showing how to write these loops.

Note that the bulk instances (``waterfluxbulk_inst``, etc.) can be obtained in two ways:

1. Using ``water_inst%water*bulk_inst``

2. As one of the indices in ``water_inst%bulk_and_tracers(:)%water*_inst``

Method (2) is just meant to be used when you are doing the same operation on bulk water
and all water tracers. Reasons why it is better, or necessary, to use method (1) when you
are really just working with bulk water are:

- This makes it more explicit and clear that you are working with bulk water

- When passing an argument to a subroutine whose dummy argument is of type
  ``water*bulk_type``, method (2) only works if you surround the call with a ``select
  type`` statement, which is awkward, to say the least. (Subroutines that expect bulk
  water should generally declare their dummy arguments to be of type ``water*bulk_type``,
  as discussed in `Object-oriented design for water tracer types`_.)

==============================================
 Infrastructure for looping through variables
==============================================

Some operations (for now just the tracer consistency check) need to be run on every water
tracer. In order to facilitate this, and to reduce the likelihood that a variable will
accidentally be left out of the tracer consistency check, we have introduced
infrastructure for holding a list of pointers to all tracer-related variables.

Variables in all of the types that apply to water tracers as well as bulk should be
allocated via the routines ``AllocateVar1d`` or ``AllocateVar2d``, rather than simple
``allocate`` statements. These routines add a pointer to the given variable in the
appropriate water tracer container structure. This allows us to later loop through all
such variables in order to call the tracer consistency check (and later, possibly other
routines as well) on each variable.
