#!/bin/bash

# creates a jobname
get_jobname () {
    if [ -f segment001.job ]; then
        segment=$(ls segment*.job | wc -l)
        ((segment++))
    else
        segment="1"
    fi
    jobname="segment"$(printf %03d $segment)
    echo $jobname
}

# adds newline if commands file does not end in newline
prep_commands () {
    if ! [[ $(tail -c1 $commands | wc -l) -gt 0 ]]; then
        echo "" >> $commands
    fi
}

# logfile header
header () {
    echo "-----------------------------"
    echo "TETHER COMMANDS:"
    echo "-----------------------------"
}

# updates logfile for a given command
command_log () {
    echo $line
    echo "-----------------------------"
    echo "with stdout:"
    echo $OUT
    echo "-----------------------------"
    echo "and stderr:"
    cat errFile
    echo "-----------------------------"
}

main () {
    wdir=$1
    commands=$2
    template=$3
    cd $wdir
    if [ -f $commands ]; then
        prep_commands
        header
        while read -r line; do
            OUT=$(command ${line} 2> errFile) #run command
            command_log #log output
        done< $commands
        rm errFile

	#read case
	case=$(<case.txt)
	
	# announce in logfile
	echo "SUBMITTING "$(pwd)"/"$case
	echo "-----------------------------"
	cd $case
	
	# submit the case, capturing jobid
	./case.submit --resubmit-immediate
	X=$(./xmlquery JOB_IDS)
	arrX=(${X//:/ })
	jobid=${arrX[-1]} #save the last jobid
	
	#presubmit tethered case
        cd $wdir
	nextjob=$(get_jobname)
	qj=$nextjob".job"
	sed 's:jobname:'$nextjob':g' $template > $qj
        sed -i 's:jobid:'$jobid':g' $qj
	sed -i 's:wdir:'$wdir':g' $qj
	sed -i 's:commands:'$commands':g' $qj
	sed -i 's:template:'$template':g' $qj
        qsub $qj

    fi
}

main $1 $2 $3
