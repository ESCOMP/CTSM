
TAG=$(git describe --tags | cut -d- -f1)
NLDIR=$(pwd)"/namelists"


sed "s/ustr/"$USER"/g" .pre_config > spinup.config
sed -i "s:pwd:"$(pwd)":g" spinup.config
sed -i "s/tagstr/"$TAG"/g" spinup.config
sed -i "s:nldir:"$NLDIR":g" spinup.config
