#$ -S /bin/sh
#$ -cwd
#$ -q all*

# Email at beginning and end of job.
#$ -m n 
#$ -M rocio@berkeley.edu 

# Run the simulation only on one of the following nodes: compute-0-1 compute-0-10 compute-0-7 compute-0-13 compute-0-8
#$ -l h="compute-0-1|compute-0-10|compute-0-7|compute-0-13|compute-0-8|compute-0-11|compute-0-12|compute-0-14|compute-0-2|compute-0-4|compute-0-5|compute-0-6|compute-0-9"

# Copy files to scratch 
echo "Job ID ${JOB_ID} "
echo "Running on host: "
echo "    $(hostname)"
echo "    with ${num_cpu} cpus."

JOB_DIR=/scratch/$USER/${JOB_ID}
echo "in directory: "
echo "    $JOB_DIR"
echo "submission directory: "
echo "    $SGE_O_WORKDIR"

mkdir -p $JOB_DIR
cp * $JOB_DIR/.

cd $JOB_DIR

# Run interpenetrated structures program
module load python
python2.7 findInterpenetratedStructures.py

cd $SGE_O_WORKDIR

# Copy files back to working directory
for i in $JOB_DIR/*
do
  mv $i .
done

