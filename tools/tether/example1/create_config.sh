
TAG=$(git describe --tags | cut -d- -f1)
NLDIR=$(pwd)"/namelists"


sed "s/ustr/"$USER"/g" .pre_config > example1.config
sed -i "s:pwd:"$(pwd)":g" example1.config
sed -i "s/tagstr/"$TAG"/g" example1.config
sed -i "s:nldir:"$NLDIR":g" example1.config
