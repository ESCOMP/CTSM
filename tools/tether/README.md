# Tether
## a utility to automate the tethering between multi-segment CTSM simulations

## Specific tether examples:

1. example1
    - a very simple example that runs two one month 4x5 simulations tethered together
    - this is merely illustrative
    - not as deeply configurable as sasu_spinup

2. sasu_spinup
    - runs a standard spinup sequence of AD->SASU->ND
    - checks for spinup stability periodically, resubmitting until stability, then proceeding to the next segment
    - I tried to prototype some config files to make this a bit more user friendly and easy to update
    - It's currently configured to run 2deg CLM60Bgc global, with loose spinup criteria

3. old_spinup
    - safe to ignore for now
    - I needed this for stuff I'm doing presently, and can eventually make this more like sasu_spinup if desired
    

## Generic tether information

### Basic tether call sequence:
1. run the commands specified in commands.txt
2. submit the case segment specified in case.txt
3. call tether.sh to presubmit your next case segment

### tether.sh syntax
 - ```./tether.sh $WDIR commands.txt $template```

### tether.sh tips
 - I prefer to write a stub shell script with a PBS header to run tether, in lieu of calling it from the command line, because:
    1. I am often running scripts that are not well-suited to the login nodes
    2. I can automatically pipe stdout and stderr to an appropriately named log file
 - As such, I will create a file: **segment001.job**, which can be a one-liner:
    - ```$TDIR/tether.sh $WDIR commands.txt $template```
    - where ```$TDIR``` provides the full path to tether.sh
    - with the addition of some standard PBS magic in the header
    - see tether/example1/segment001.job
 - then to run tether, I will opt for:
    - ```qsub segment001.job```
  
### $WDIR
 - the working directory within which all your cases and log files will exist

### commands.txt
 - a text file that will be run as bash commands line by line
 - example:
```
./part1.sh
```

### commands.txt tips
 - commands.txt should exist within ```$WDIR```
 - ```./``` refers to ```$WDIR```, any other directories should utilize full path
 - my intuition generally is to call a single script here, but I could see how eventually we may have smaller reusable scripts that are called in sequence, so commands.txt can have as many lines as needed
 - within the commands, you should create, setup, build, and customize whatever case you are planning to submit
 - somewhere within the commands, you need to write case.txt to indicate which case will be submitted in the current call to tether
 - all cases should be created in ```$WDIR```
 - you should **NOT** call ```case.submit``` within your commands, tether.sh is responsible for case submission
 - within the final command, you need to rewrite commands.txt for your next segment
    - e.g. the last line of part1.sh is: ```echo "./part2.sh" > commands.txt```
 - if there are no future segments, ```rm commands.txt``` will signal tether to exit after the current job submission
    - this is pretty important, because as currently written it's pretty easy to accidentally submit an infinite job sequence to the queue (**not advised**)

### $template
 - a stub shell script with PBS magic that is used to presubmit the next tethered segment
 - see ```tether/derecho.template``` for an example
 - you should only need to edit the project code within this template
 - this will end up looking very similar to segment001.job, apart from an extra PBS command with the ```afterok``` flag

### what it might look like to configure a tethered sequence:
 1. mkdir $WDIR
 2. write the necessary case configuration scripts
     - taking care to include the tethering step of editing commands.txt
     - and writing the casename to case.txt
 3. write the initial commands.txt
     - typically just a one-liner, e.g.: ./part1.sh
 4. prepare segment001.job
 5. qsub segment001.job

### what it looks like when tether is running:
 - if you run ```qstat -u $USER```, you'll probably see segment001 queued or running
 - as segment001 is running it will write log files with stdout and stderr in ```$WDIR```
 - it will also prepare segment002.job, which you can inspect but won't be very interesting
 - more interesting would be ```cat commands.txt```
    - this should eventually update to reflect the commands required for your second segment
 - then after segment001 is finished, if you rerun ```qstat -u $USER``` you may see something like:

                                      
| Job id   | Name     | User     | Time Use | S | Queue
| :------- | :------: | :------: | :------: |:------: |-------: |
| 2542575.desched1 | run.mycase | djk2120  | | Q | cpu
| 2542576.desched1 | st_archive.m | djk2120  | | H | cpu
| 2542577.desched1 | segment002 | djk2120  | | H | cpudev

 - The Q means queued, and the H's mean on hold. in this case job 2542576 is waiting for 2542575 to finish, and 2542577 is waiting for 2542576 to finish
 - In this situation once mycase finishes running, and the short-term archiver completes, segment002.job will submit the next case in the tethered sequence

### when things go wrong
 - you may end up needing to frantically qdel
 - e.g. ```qdel 2542577```
 - useful but very committing is:
     - ```qselect -u $USER | xargs qdel```
     - which will kill all your derecho jobs
 - alternatively you can pipe to tail and head as appropriate
     - ```qselect -u $USER | head -n 5 | tail -n 3``` etc.
     - until you isolate the set of jobs that you want to delete
 - and then call:
     - ```qselect -u $USER | head -n 5 | tail -n 3 | xargs qdel```


### restarting tether, after a debugging fix
 - clean house as needed, e.g. you may need to rm -r mycase and rm -r /glade/derecho/scratch/USER/mycase     
 - edit commands.txt to reference the script you want to start with
 - edit case.txt to point to the case that will be submitted
 - remove the afterok line from the last segment0XX.job
 - qsub segment0XX.job

