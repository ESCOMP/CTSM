

# apply defaults to any unset config fields
python setup_config.py
rm example1.config  # old unset file

# source variables from config
source <(grep = config.tmp)

# mk WDIR and mv config to WDIR
mkdir -p $WDIR
mv config.tmp $WDIR"/example1.config"
echo "config file now located at: "$WDIR"/example1.config"


# make sure scratch directory exists
SD='/glade/derecho/scratch/'$USER
if [ ! -d $SD ]; then
    echo $SD" is not a directory, $USER may be set badly"
    echo "examine "$WDIR"/spinup.config"
    exit 1
fi


# cp scripts to WDIR
cp part*.sh $WDIR


# configure the PBS template
sed "s:TDIR:"$TDIR":g" derecho.template > $WDIR"/derecho.template"
sed -i "s:project:"$PROJECT":g" $WDIR"/derecho.template"


# create the initial qsub job
sed '/afterok/d' $WDIR"/derecho.template" > $WDIR"/segment001.job"
sed -i 's/jobname/segment001/g' $WDIR"/segment001.job"
sed -i 's:wdir:'$WDIR':g' $WDIR"/segment001.job"
sed -i 's/commands/commands.txt/g' $WDIR"/segment001.job"
sed -i 's/template/derecho.template/g' $WDIR"/segment001.job"


# create the initial commands.txt
echo "./part1.sh example1.config" > $WDIR"/commands.txt"

