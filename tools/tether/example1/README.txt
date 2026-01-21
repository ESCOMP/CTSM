This is a simple demonstration example of a two-segment tethered simulation.

This will run a 1-month 4x5deg CLM-SP simulation, then spawn a hybrid case that picks up and continues that for another month.

Steps to make this work for you:
1) ./create_config.sh
2) ./setup_config.sh
3) cd to tethered_sims/c2025xxxx  directory
4) qsub segment001.job
