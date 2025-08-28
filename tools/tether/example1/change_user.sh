
# here are some lines to apply your username in the scripts
# replace djk2120 with your username
sed -i 's/USER/djk2120/g' part1.sh
sed -i 's/USER/djk2120/g' part2.sh
sed -i 's/USER/djk2120/g' reset.sh

# here are some lines to apply your project code in the scripts
# replace abcdefg with your code
# or comment out if you have access to the LMWG code
sed -i 's/P93300041/abcdefg/g' part1.sh
sed -i 's/P93300041/abcdefg/g' part2.sh
sed -i 's/P93300041/abcdefg/g' segment001.job
sed -i 's/P93300041/abcdefg/g' derecho.template
