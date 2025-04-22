#!/bin/bash

#SBATCH --output outputs/log/mtp-oncho.out-%A_%a
#SBATCH --array 370,373,474,498,499,500,502,503,504,505,506,507,511,512,513,514,515,516,517,518,519,566,570,589,590,591,592,593,596,600,601,602,603,604,617,624,631,633,642,644,645,646,647,202,307,309,371,372,496,497,501,509,510,580,583,584,585,588,594,598,606,607,608,609,610,612,613,620,621,622,623,625,626,628,630,634,635,637,639,640,641,643,648,310,473,508,537,586,595,597,605,615,618,627,629,636,638,308,461,539,632,611,614,616,582
#SBATCH --nodes=6
#SBATCH --cpus-per-task 12
#SBATCH --time=24:00:00

#####################
# Your task is here #
#####################
# Add extra commands here to load a recent version of R
module purge
module load GCC/11.3.0 OpenMPI/4.1.4 R/4.2.1

stdbuf -i0 -o0 -e0 command
####################
# End of your task #
####################

# Now that you have loaded R above, we can run our R script
Rscript oncho-endgame-multipletimepts.R ${SLURM_ARRAY_TASK_ID} 

