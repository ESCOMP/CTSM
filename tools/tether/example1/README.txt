This is a simple demonstration example of a two-segment tethered simulation.

This will run a 1-month 4x5deg CLM-SP simulation, then spawn a hybrid case that picks up and continues that for another month.

Steps to make this work for you:
1) open and edit change_user.sh
2) ./change_user.sh
3) manually update tdir and wdir in derecho.template and segment001.job 
4) qsub segment001.job
