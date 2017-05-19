#$ -S /bin/sh
#$ -cwd
#$ -q all*

# Email at beginning and end of job.
#$ -m n 
#$ -M rocio@berkeley.edu 

JOBSDIR=jobs

# Test for jobs directory. If not available, create it.

if [ ! -d ${HOME}/${JOBSDIR} ]; then
    mkdir ${HOME}/${JOBSDIR}
fi
if [ ! -d ${HOME}/${JOBSDIR}/running ]; then
    mkdir ${HOME}/${JOBSDIR}/running
fi
if [ ! -d ${HOME}/${JOBSDIR}/completed ]; then
    mkdir ${HOME}/${JOBSDIR}/completed
fi
if [ ! -e ${HOME}/${JOBSDIR}/jobs.log ]; then
    touch ${HOME}/${JOBSDIR}/jobs.log
fi

# Record the current location.
CURDIR=`pwd`

TIMESTAMP=`date '+%m-%d-%y %H:%M:%S'`

# Link to the running jobs directory

ln -s "${CURDIR}" ${HOME}/${JOBSDIR}/running/${JOB_ID}-${JOB_NAME}

echo "${TIMESTAMP} [SGE] Job ${JOB_ID}: ${JOB_NAME} started on ${HOSTNAME}" >> ${HOME}/${JOBSDIR}/jobs.log

# Run the simulation
#export OMP_NUM_THREADS=1

module load python
python2.7 findInterpenetratedStructures.py


ENDTIME=`date '+%m-%d-%y %H:%M:%S'`

echo "${ENDTIME} [SGE] Job ${JOB_ID}: ${JOB_NAME} completed" >> ${HOME}/${JOBSDIR}/jobs.log

rm -f ${HOME}/${JOBSDIR}/running/${JOB_ID}-${JOB_NAME}
ln -s "${CURDIR}" ${HOME}/${JOBSDIR}/completed/${JOB_ID}-${JOB_NAME}

