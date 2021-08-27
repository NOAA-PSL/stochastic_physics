#!/bin/sh
#SBATCH -e err
#SBATCH -o out
#SBATCH --account=gsienkf
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=20
#SBATCH --job-name="stoch_unit_tests"
DO_CA_SGS=.true.
DO_CA_GLOBAL=.true.
source ./module-setup.sh
module purge
module use $( pwd -P )
module load modules.stoch
EXEC=standalone_ca.x
# compile codes
sh compile_standalone_ca.hera_intel
if [ ! -f $EXEC ];then
  echo "compilation errors"
  exit 1
fi
sh compile_compare_ca.sh

# copy input directory
if [ ! -d INPUT ]; then
   cp -r /scratch2/BMC/gsienkf/Philip.Pegion/stochastic_physics_unit_tests/input_data INPUT
fi
mkdir -p RESTART

#layout 1x1
cp input.nml.ca_template input.nml
sed -i -e "s/NPX/1/g" input.nml
sed -i -e "s/NPY/1/g" input.nml
sed -i -e "s/CA_SGS/${DO_CA_SGS}/g" input.nml
sed -i -e "s/CA_GLOBAL/${DO_CA_GLOBAL}/g" input.nml
sed -i -e "s/WARM_START/.false./g" input.nml
export OMP_NUM_THREADS=1
time srun --label -n 6 $EXEC >& stdout.1x1
ls -l RESTART >> stdout.1x1
exit
mkdir ca_layout_1x1
mv ca_out* ca_layout_1x1
ct=1
#mv RESTART/ca_data.nc RESTART/run1_end_ca_data.nc
while [ $ct -le 6 ];do
   mv RESTART/ca_data.tile${ct}.nc RESTART/run1_end_ca_data.tile${ct}.nc
   ct=`expr $ct + 1`
done

# test 3 different domain decompositions and compare to baseline
#layout 1x4
cp input.nml.ca_template input.nml
sed -i -e "s/NPX/1/g" input.nml
sed -i -e "s/NPY/4/g" input.nml
sed -i -e "s/CA_SGS/${DO_CA_SGS}/g" input.nml
sed -i -e "s/CA_GLOBAL/${DO_CA_GLOBAL}/g" input.nml
sed -i -e "s/WARM_START/.false./g" input.nml
export OMP_NUM_THREADS=1
time srun --label -n 24 $EXEC  >& stdout.1x4
ls -l RESTART >> stdout.1x4
mkdir ca_layout_1x4
mv ca_out* ca_layout_1x4
ct=1
#mv RESTART/ca_data.nc RESTART/run2_end_ca_data.nc
while [ $ct -le 6 ];do
   mv RESTART/ca_data.tile${ct}.nc RESTART/run2_end_ca_data.tile${ct}.nc
   ct=`expr $ct + 1`
done
#layout 2x2
export OMP_NUM_THREADS=2
cp input.nml.ca_template input.nml
sed -i -e "s/NPX/2/g" input.nml
sed -i -e "s/NPY/2/g" input.nml
sed -i -e "s/CA_SGS/${DO_CA_SGS}/g" input.nml
sed -i -e "s/CA_GLOBAL/${DO_CA_GLOBAL}/g" input.nml
sed -i -e "s/WARM_START/.false./g" input.nml
time srun --label -n 24 $EXEC   >& stdout.2x2
ls -l RESTART >> stdout.2x2
mkdir ca_layout_2x2
mv ca_out* ca_layout_2x2
ct=1
#mv RESTART/ca_data.nc RESTART/run3_end_ca_data.nc
while [ $ct -le 6 ];do
   mv RESTART/ca_data.tile${ct}.nc RESTART/run3_end_ca_data.tile${ct}.nc
   ct=`expr $ct + 1`
done
##layout 1x4
export OMP_NUM_THREADS=1
cp input.nml.ca_template input.nml
sed -i -e "s/NPX/4/g" input.nml
sed -i -e "s/NPY/1/g" input.nml
sed -i -e "s/CA_SGS/${DO_CA_SGS}/g" input.nml
sed -i -e "s/CA_GLOBAL/${DO_CA_GLOBAL}/g" input.nml
sed -i -e "s/WARM_START/.false./g" input.nml
time srun --label -n 24 $EXEC   >& stdout.4x1
ls -l RESTART >> stdout.4x1
mkdir ca_layout_4x1
mv ca_out* ca_layout_4x1
ct=1
#mv RESTART/ca_data.nc RESTART/run4_end_ca_data.nc
while [ $ct -le 6 ];do
   mv RESTART/ca_data.tile${ct}.nc RESTART/run4_end_ca_data.tile${ct}.nc
   ct=`expr $ct + 1`
done
# restart run
ct=1
#mv RESTART/mid_run.ca_data.nc INPUT/ca_data.nc
while [ $ct -le 6 ];do
   mv RESTART/mid_run.ca_data.tile${ct}.nc INPUT/ca_data.tile${ct}.nc
   ct=`expr $ct + 1`
done
export OMP_NUM_THREADS=1
cp input.nml.ca_template input.nml
sed -i -e "s/NPX/1/g" input.nml
sed -i -e "s/NPY/1/g" input.nml
sed -i -e "s/CA_SGS/${DO_CA_SGS}/g" input.nml
sed -i -e "s/CA_GLOBAL/${DO_CA_GLOBAL}/g" input.nml
sed -i -e "s/WARM_START/.true./g" input.nml
time srun --label -n 6 $EXEC >& stdout.1x1_restart
mkdir ca_layout_1x1_restart
mv ca_out* ca_layout_1x1_restart

compare_ca_output
if [ $? -ne 0 ];then
   echo "unit tests failed"
else
   ct=1
   while [ $ct -le 6 ];do
      diff RESTART/ca_data.tile${ct}.nc RESTART/run1_end_ca_data.tile${ct}.nc
      if [ $? -ne 0 ];then
         echo "restart test failed"
         exit 1
      fi   
      ct=`expr $ct + 1`
   done
   if [ $? -eq 0 ];then
      echo "unit tests successful"
exit
      rm -rf ca_layout_*
      rm logfile*
      rm ../*.o ../*.mod
      rm $EXEC
   fi
fi
