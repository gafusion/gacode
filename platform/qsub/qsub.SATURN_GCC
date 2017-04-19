
bfile=$SIMDIR/batch.src
echo "#!/bin/bash -l" > $bfile
echo "#SBATCH -J $LOCDIR" >> $bfile
echo "#SBATCH -o $SIMDIR/batch.out" >> $bfile
echo "#SBATCH -e $SIMDIR/batch.err" >> $bfile
echo "#SBATCH -t $WALLTIME" >> $bfile
echo "#SBATCH -n $cores_used" >> $bfile
echo "#SBATCH --mem 45G" >> $bfile
if [ "$QUEUE" = "null_queue" ]
then
  echo "#SBATCH -p medium" >> $bfile
else
  echo "#SBATCH -p $QUEUE" >> $bfile
fi
echo "$CODE -e $LOCDIR -n $nmpi -nomp $nomp -numa $numa -mpinuma $mpinuma -p $SIMROOT" >> $bfile
