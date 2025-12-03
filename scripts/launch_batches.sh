#!/bin/bash
SCAN_NAME="shepherding"
BASE_PARAMS="k401.yaml"
echo "Launching ${SCAN_NAME} scan with ${BASE_PARAMS} parameter file."

I_BATCH=0

for XLINK_CONC in 1 0.1
do
	for MOT_CONC in 100 10
	do
		for N_PFS in 8
		do
			for N_SITES in 1000 # 250 500 750 1000 1250
			do
				BASE_NAME="shep_${XLINK_CONC}nM_${MOT_CONC}nM_${N_PFS}_${N_SITES}"
				echo ${BASE_NAME}
				PARAM_FILE="params_temp_${BASE_NAME}.yaml"
				cp ${BASE_PARAMS} ${PARAM_FILE}
				yq eval -i ".xlinks.c_bulk = ${XLINK_CONC}" ${PARAM_FILE}
				yq eval -i ".motors.c_bulk = ${MOT_CONC}" ${PARAM_FILE}
				yq eval -i ".filaments.n_subfilaments = ${N_PFS}" ${PARAM_FILE}
				yq eval -i ".filaments.n_sites[0] = ${N_SITES}" ${PARAM_FILE}
				I_BATCH=$(( $I_BATCH + 1 ))
				JOB_NAME="shep${I_BATCH}"
				sbatch -J ${JOB_NAME} SubmitJob.slurm ${BASE_NAME} ${PARAM_FILE}
			done
		done
	done
done
