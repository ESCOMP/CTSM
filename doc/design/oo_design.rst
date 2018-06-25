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

At this point, we could probably refactor this to use constructors.
