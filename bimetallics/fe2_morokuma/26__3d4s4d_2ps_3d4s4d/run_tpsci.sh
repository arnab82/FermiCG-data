#!/bin/bash
#SBATCH --nodes=1
#SBATCH -p preemptable_q
##SBATCH -p normal_q
#SBATCH --cpus-per-task=1
#SBATCH --mem=230GB
#SBATCH -t 12:00:00
#SBATCH --account=personal
#SBATCH --job-name cmx_bst.jl
##SBATCH --exclusive

export NTHREAD=1
export JULIAENV=/home/arnab22/bst_cmx.jl

# allow time to load modules, this takes time sometimes for some reason
sleep 10

hostname

module reset
module load site/tinkercliffs-rome/easybuild/setup
module load site/tinkercliffs/easybuild/setup
module load Anaconda3/2020.07
module load gcc/8.2.0
#module load Python/3.8.6-GCCcore-10.2.0
#module load Julia/1.7.2-linux-x86_64

source activate bst_cmx

echo "Usage: sbatch submit.sh {input file} {data file} {data file}"

# set these to 1 if you are using julia threads heavily to avoid oversubscription
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export INFILE=$1
export OUTFILE="${INFILE}.out"
export WORKDIR=$(pwd)

echo $INFILE
echo $OUTFILE
echo $WORKDIR
echo $TMPDIR

cp $INFILE $TMPDIR/
if [ "$2" ]
then
        cp $2 $TMPDIR/
fi
if [ "$3" ]
then
        cp $3 $TMPDIR/
fi
cd $TMPDIR


#Start an rsync command which runs in the background and keeps a local version of the output file up to date
touch $OUTFILE
while true; do rsync -av $OUTFILE $WORKDIR/"${INFILE}.${SLURM_JOB_ID}.out"; sleep 60; done &


julia --project=$JULIAENV -t $NTHREAD $INFILE >& $OUTFILE

cp $OUTFILE $WORKDIR/"${INFILE}.out"
rm $WORKDIR/"${INFILE}.${SLURM_JOB_ID}.out"

mkdir $WORKDIR/"${INFILE}.${SLURM_JOB_ID}.scr"

cp -r * $WORKDIR/"${INFILE}.${SLURM_JOB_ID}.scr"
rm -r *

#moving standard output slurm file to specific job directory
mv $WORKDIR/"slurm-${SLURM_JOB_ID}.out" $WORKDIR/"${INFILE}.${SLURM_JOB_ID}.scr"

exit
