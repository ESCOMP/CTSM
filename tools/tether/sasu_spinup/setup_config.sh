

# apply defaults to any unset config fields
python setup_config.py
rm spinup.config  # old unset file

# source variables from config
source <(grep = config.tmp)

# mk WDIR and mv config to WDIR
mkdir -p $WDIR
mv config.tmp $WDIR"/spinup.config"
echo "config file now located at: "$WDIR"/spinup.config"


# make sure scratch directory exists
SD='/glade/derecho/scratch/'$USER
if [ ! -d $SD ]; then
    echo $SD" is not a directory, $USER may be set badly"
    echo "examine "$WDIR"/spinup.config"
    exit 1
fi

# cp scripts to WDIR
cp *AD.sh $WDIR
cp *SASU.sh $WDIR
cp *ND.sh $WDIR


# add a path comment to the namelist files
#  it probably makes sense to move these to WDIR
segments=('AD' 'SASU' 'ND')
cd namelists
for segment in ${segments[@]}; do
    cd $segment
    sed 's:nlpath:'$(pwd)':g' unset_user_nl_clm > user_nl_clm
    cd ..
done
cd ..


# configure the PBS template
sed "s:TDIR:"$TDIR":g" derecho.template > $WDIR"/derecho.template"
sed -i "s:project:"$PROJECT":g" $WDIR"/derecho.template"


# configure yamls
HD=$SD"/archive/"$CASE_AD"/lnd/hist"
sed 's/CASE_AD/'$CASE_AD'/g' AD.yml > $WDIR"/AD.yml"
sed -i 's:HIST_AD:'$HD':g' $WDIR"/AD.yml"
HD=$SD"/archive/"$CASE_SASU"/lnd/hist"
sed 's/CASE_SASU/'$CASE_SASU'/g' SASU.yml > $WDIR"/SASU.yml"
sed -i 's:HIST_SASU:'$HD':g' $WDIR"/SASU.yml"
HD=$SD"/archive/"$CASE_ND"/lnd/hist"
sed 's/CASE_ND/'$CASE_ND'/g' ND.yml > $WDIR"/ND.yml"
sed -i 's:HIST_ND:'$HD':g' $WDIR"/ND.yml"


# create the initial qsub job
sed '/afterok/d' $WDIR"/derecho.template" > $WDIR"/segment001.job"
sed -i 's/jobname/segment001/g' $WDIR"/segment001.job"
sed -i 's:wdir:'$WDIR':g' $WDIR"/segment001.job"
sed -i 's/commands/commands.txt/g' $WDIR"/segment001.job"
sed -i 's/template/derecho.template/g' $WDIR"/segment001.job"


# create the initial commands.txt
echo "./setupAD.sh spinup.config" > $WDIR"/commands.txt"

