#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o job_out/PRN_$JOB_ID.txt
#$ -j y
## Edit the line below as needed: h_data = 32G
#$ -l h_rt=24:00:00,h_data=8G
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 1
# Email address to notify
#$ -M $USER@mail 
# Notify when
#$ -m bea

echo $PATH

echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

source ~/.bashrc
module load anaconda3
which python


echo 'python run_PageRank_Nibble.py'
python run_PageRank_Nibble.py