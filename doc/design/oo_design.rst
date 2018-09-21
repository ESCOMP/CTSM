.. sectnum::

.. contents::

=============================================================
 Use of Init method rather than constructor for most classes
=============================================================

Most of the object-oriented classes in CTSM (and particularly the science-focused classes)
are initialized with a method named ``Init``, rather than the more standard
object-oriented pattern of using a constructor. This is largely for historical reasons:
Object initialization was done this way when object orientation was first introduced to
CESM (possibly because of compiler bugs that prevented the general use of constructors for
this purpose?). As more object orientation was added, we continued to use an ``Init``
method for this purpose to remain consistent with existing code.

We could probably refactor this to use constructors. However, we still occasionally run
into trouble with some compilers in the assignment that results from writing::

  foo_inst = foo_type(...)

Some components seem to get set to garbage occasionally with some compilers when doing
that. (This was encountered 2018-09-11 with intel17.0.1). It's possible that this is user
error rather than a compiler bug, but in any case, it's harder to get that to work
robustly. So for now, we're generally sticking with using::

  call foo_inst%Init(...)

rather than::

  foo_inst = foo_type(...)
