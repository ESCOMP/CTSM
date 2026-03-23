$CTSMROOT/README.CHECKLIST.new\_case                       03/01/2021

This is a check list of things to do when setting up a new case in order to help ensure everything is correct. There are lots of tiny details that need to be right and it's easy to get something wrong. So the first screening to make sure it's right is for you to carefully check through your case and make sure it's right.

The following assumes you have created a new case and are in it's case directory.

General Checklist to always do:

- Make sure CLM\_ env settings are correct: ./xmlquery \-p CLM  
- Make sure you are using the correct CLM\_PHYSICS\_VERSION: ./xmlquery \-p CLM\_PHYSICS\_VERSION  
- Make sure you are running the appropriate overall CLM vegetation model, i.e. the "-bgc" option of either Satellite Phenology (sp) or Full BioGeoChemistry (bgc) or FATES (fates): ./xmlquery \-p CLM\_BLDNML\_OPTS
- If you are running the bgc model, check to see if you should be running the prognostic crop model: option \-crop in CLM\_BLDNML\_OPTS
- Make sure the LND\_TUNING\_MODE is correct: ./xmlquery LND\_TUNING\_MODE
- For an "I compset" make sure you are running over the correct forcing years: usually ./xmlquery \-p DATM\_YR  
- For an "I compset" make sure the DATM streams are operating over the correct years: look at the CaseDocs/datm.streams.xml file
- First and align year for streams should be the start year of a historical simulation: ./xmlquery RUN\_STARTDATE; grep stream\_year\_first CaseDocs/lnd\_in; grep model\_year\_align CaseDocs/lnd\_in
- Last year for streams should be the last year you are going to run to (or beyond it): grep stream\_year\_last CaseDocs/lnd\_in
- Make sure you are starting from appropriate spunup initial conditions:
    - Check the run-type with: ./xmlquery RUN\_TYPE
    - Check finidat for a startup or hybrid simulation: grep finidat CaseDocs/lnd\_in
    - Check nrevsn for a branch simulation: grep nrevsn CaseDocs/lnd\_in  
- Run for a month (or some short period) and go over the log files and especially the settings and files read in them: for an I case you especially want to look at the lnd.log and atm.log files

Some other suggestions on things that can be done:

- Compare namelist files to an existing case if you are doing something almost the same as a previous simulation
- Ask another collaborator to look over your case directory