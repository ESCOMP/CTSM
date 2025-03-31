.. include:: ../substitutions.rst

.. _comparing-history-files:

=========================
 Comparing History Files
=========================

``cprnc`` is a tool shared by |cesmrelease| to compare two NetCDF history files. It differences every field that is shared on both files, and reports a summary of the difference. The summary includes the three largest differences, as well as the root mean square (RMS) difference. It also gives some summary information on the field as well. You have to enter at least one file, and up to two files. With one file it gives you summary information on the file, and with two it gives you information on the differences between the two. At the end it will give you a summary of the fields compared and how many fields were different and how many were identical.

Options:

-m = do NOT align time-stamps before comparing

-v = verbose output

-ipr

-jpr

-kpr

See the ``cprnc`` `README <https://github.com/ESMCI/cime/blob/master/CIME/non_py/cprnc/README>`_ file for more details.

