#!/bin/sh
#SBATCH -e err
#SBATCH -o out
#SBATCH --account=gsienkf
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=30
#SBATCH --job-name="stoch_unit_tests"

source ./module-setup.sh
module purge
module use $( pwd -P )
module load modules.stoch

# compile codes
sh compile_standalone.hera_intel
sh compile_compare.sh

# test 3 different domain decompositions and compare to baseline
#layout 1x4
cp input.nml.template input.nml
sed -i -e "s/NPX/1/g" input.nml
sed -i -e "s/NPY/4/g" input.nml
srun -n 24 standalone_stochy.x
mkdir layout_1x4
mv workg* layout_1x4

#layout 2x2
cp input.nml.template input.nml
sed -i -e "s/NPX/2/g" input.nml
sed -i -e "s/NPY/2/g" input.nml
srun -n 24 standalone_stochy.x
mkdir layout_2x2
mv workg* layout_2x2
#layout 1x4
cp input.nml.template input.nml
sed -i -e "s/NPX/4/g" input.nml
sed -i -e "s/NPY/1/g" input.nml
srun -n 24 standalone_stochy.x
mkdir layout_4x1
mv workg* layout_4x1

compare_output
if [ $? -ne 0 ];then
   echo "unit tests failed"
else
   echo "unit tests successful"
   rm -rf layout_*
   rm logfile*
fi
