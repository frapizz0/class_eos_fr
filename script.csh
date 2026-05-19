#! /bin/csh -f

source /mirror/scratch/pace/Work/CODES/plc-2.0/bin/clik_profile.csh

setenv OMP_NUM_THREADS 12

echo "Program started at: `date`"

./class LCDM_Comparison.ini

echo "Program finished with exit code $? at: `date`"
