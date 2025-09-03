

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

# cp other files to WDIR
cp setupAD.sh $WDIR

# configure some various files with sed replacements
sed -i "s:nldir:"$NLDIR":g" ./namelists/AD/user_nl_clm

# configure yamls with sed replacements
HD=$SD"/archive/"$CASE_AD"/lnd/hist"
sed 's/CASE_AD/'$CASE_AD'/g' AD.yml > $WDIR"/AD.yml"
sed -i 's:HIST_AD:'$HD':g' $WDIR"/AD.yml"
HD=$SD"/archive/"$CASE_SASU"/lnd/hist"
sed 's/CASE_SASU/'$CASE_SASU'/g' SASU.yml > $WDIR"/SASU.yml"
sed -i 's:HIST_SASU:'$HD':g' $WDIR"/SASU.yml"
HD=$SD"/archive/"$CASE_ND"/lnd/hist"
sed 's/CASE_ND/'$CASE_ND'/g' ND.yml > $WDIR"/ND.yml"
sed -i 's:HIST_ND:'$HD':g' $WDIR"/ND.yml"




