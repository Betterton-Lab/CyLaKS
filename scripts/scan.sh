#!/bin/bash
BASE_NAME=$1
BASE_PARAMS=$2 # "k401.yaml"
echo "Starting ${BASE_NAME} scan"
echo "Base parameter file is ${BASE_PARAMS}"

calc() { awk "BEGIN{ printf \"%.6f\n\", $* }"; }

BASE_SEED=198261346419

BASE_DIFF=0.131
BASE_KON=0.000238
BASE_KOFF=0.03

BASE_KON_MOT=0.000358
BASE_KOFF_I_MOT=10
BASE_KOFF_II_MOT=260

for E_INT in 0.2 0.4 0.6 0.8
do
	for BIND_AFF in 3 # 2 3 # 4 6 8 
	do
		for I_SEED in 0 # 1 2 3 4
		do
			SIM_NAME="${BASE_NAME}_${E_INT}kT_${BIND_AFF}x_${I_SEED}"
			PARAM_FILE="params_temp_${SIM_NAME}.yaml"
			cp ${BASE_PARAMS} ${PARAM_FILE}
			yq eval -i ".xlinks.neighb_neighb_energy = ${E_INT}" ${PARAM_FILE}
		#	yq eval -i ".xlinks.d_i = $(calc ${BASE_DIFF}/${BIND_AFF})" ${PARAM_FILE}
		#	yq eval -i ".xlinks.d_side = $(calc ${BASE_DIFF}/${BIND_AFF})" ${PARAM_FILE}
		#	yq eval -i ".xlinks.k_off_i = $(calc ${BASE_KOFF}/${BIND_AFF})" ${PARAM_FILE}
		#	yq eval -i ".xlinks.k_on = $(calc ${BASE_KON}*${BIND_AFF})" ${PARAM_FILE}
			#yq eval -i ".motors.k_on = $(calc ${BASE_KON_MOT}*${BIND_AFF})" ${PARAM_FILE}
			#yq eval -i ".motors.k_off_i = $(calc ${BASE_KOFF_I_MOT}/${BIND_AFF})" ${PARAM_FILE}
			#yq eval -i ".motors.k_off_ii = $(calc ${BASE_KOFF_II_MOT}/${BIND_AFF})" ${PARAM_FILE}
		    	yq eval -i ".seed = $(( ${BASE_SEED} + ${I_SEED} ))" ${PARAM_FILE}
			# Run sim for these parameter values
			echo "Launching sim ${SIM_NAME} with parameter file ${PARAM_FILE}"
			#singularity exec --bind $PWD cylaks_latest.sif cylaks.exe ${PARAM_FILE} ${SIM_NAME} &
			apptainer exec --bind $PWD cylaks_app.sif cylaks.exe ${PARAM_FILE} ${SIM_NAME} &
		done
	done
done
wait
#rm params_temp_${BASE_NAME}_*
