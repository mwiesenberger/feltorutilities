#!/bin/bash

FILE="${2%.*}" # outfilename without extension
JOBNAME=${FILE: -8} # last 8 characters

# help simplesimdb recognize that data is being generated
touch $2

if [ $# -eq 2 ] # note the whitespaces
then
    # $@ forwards all arguments
    sbatch -o ${FILE}.out -J $JOBNAME  job.sh "$@"
else
    # do something else for restart
    PREVIOUS="${3%.*}" # previous outfile without extension
    PREVNAME=${PREVIOUS: -8} # last 8 characters
    # check if a job with that name exists
    JOBID=$(squeue --me --noheader --format="%i" --name "$PREVNAME")
    # --me only display owned jobs
    # --noheader suppress the header
    # --format="%i" only show the jobid
    # --name only show jobs with specified name
    if [ -z $JOBID]
    then
        sbatch -o ${FILE}.out -J $JOBNAME  job.sh "$@"
    else
        sbatch --dependency=afterany:${JOBID} -o ${FILE}.out \
               -J $JOBNAME job.sh "$@"
    fi
fi
