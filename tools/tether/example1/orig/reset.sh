#this will reset example1 to before change_user, etc.
echo "bash part1.sh" > commands.txt
rm -rf I1850.f45_g37.txd.*
rm -rf /glade/derecho/scratch/USER/I1850.f45_g37.txd.*
rm -rf /glade/derecho/scratch/USER/archive/I1850.f45_g37.txd.*
rm case.txt
rm segment002.job
rm segment003.job

cp orig/* ./
