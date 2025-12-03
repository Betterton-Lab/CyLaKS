#!/bin/bash

# === List of SBATCH arguments ===
#SBATCH --account=ucb370_asc2		# allocation account
#SBATCH --partition=amilan		# Type of node to run on
#SBATCH --qos=normal #long #normal        		# quality of service/queue
#SBATCH --nodes=1			# number of cluster nodes, abbrv. by -N
#SBATCH --ntasks=4			# number of parallel process
#SBATCH --cpus-per-task=1           	# cpu allocation per process
#SBATCH --time=24:00:00 #168:00:00			# walltime, abbrv. by -t
#SBATCH --job-name=shep			# name of job
#SBATCH --output=run_shep_%j.out 	# name of the stdout redirection file
#SBATCH --error=run_shep_%j.err	    	# name of the stderr redirection file

# === Purge all modules and load needed ones ===
module purge
#module load singularity
module load apptainer

# === Run the program ===
echo "Current variables CC=${CC} and CXX=${CXX}"
echo "Running on $(hostname --fqdn)"

#apptainer exec --bind $PWD cylaks_latest.sif cylaks.exe k401.yaml test # ${PARAM_FILE} ${SIM_NAME} &
bash scan.sh $1 $2
