.. comparing-history-files:

=========================
 Comparing History Files
=========================

**cprnc** is a tool shared by both CAM and CLM to compare two NetCDF history files. 
It differences every field that has a time-axis that is also shared on both files, and reports a summary of the difference. 
The summary includes the three largest differences, as well as the root mean square (RMS) difference. 
It also gives some summary information on the field as well. 
You have to enter at least one file, and up to two files. 
With one file it gives you summary information on the file, and with two it gives you information on the differences between the two. 
At the end it will give you a summary of the fields compared and how many fields were different and how many were identical.

Options:

-m = do NOT align time-stamps before comparing

-v = verbose output

-ipr

-jpr

-kpr

See the **cprnc** `README <CLM-URL>`_ file for more details.

.. note:: To compare files with OUT a time axis you can use the **cprnc.ncl** NCL script in ``$CTSMROOT/tools/ncl_scripts``. It won't give you the details on the differences but will report if the files are identical or different.
