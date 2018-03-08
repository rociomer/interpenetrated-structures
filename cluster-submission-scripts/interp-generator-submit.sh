#$ -S /bin/sh
#$ -cwd
#$ -q all*

# Email at beginning and end of job.
#$ -m n 
#$ -M 

# Run the simulation only on the following nodes: 
#$ -l h="compute-0-1|compute-0-10|compute-0-7|compute-0-13"

# Copy files to scratch 
echo "Job ID ${JOB_ID} "
echo "Running on host: "
echo "    $(hostname)"

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

