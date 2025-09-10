This is a more complicated sequence. And I didn't make this as friendly as example 1. It would be interesting to see how we could reconfigure this use scenario into an actual set of tools. My plan is to write this a bit more generically with a newer tag, using the AD/SASU/pSASU/transient sequence. This current example is for actual work that I'm doing right now.

It feels like we could pull out things like CIMEROOT COMPSET GRID etc into a config file, and then have a utility that would cp the shell scripts and yamls into a user-specified WDIR, and give some instructions on how to proceed and the various options for customization. We could also add a line in the script for user_mods so user hopefully wouldn't need to edit the main script very often.

The general sequence is:
 - Run AD for 100 years, from an existing restart
 - Test for spinup stability and run 20 more years, as needed
 - Run postAD for 100 years, picking up from the AD case
 - Test for spinup stability and run 20 more years, as needed
 - Run a transient case from 1880-1949, picking up from the postAD case
