#!/bin/bash
BASE_NAME="shep"
BASE_PARAMS="params/shepherding.yaml"
echo "Starting ${BASE_NAME} scan"
echo "Base parameter file is ${BASE_PARAMS}"

calc() { awk "BEGIN{ printf \"%.6f\n\", $* }"; }

BASE_SEED=198261346419
BASE_DIFF=0.131

for SCALE_LONG in 0.03 0.1 0.3 1 3 10 30
do
	for SCALE_LAT in 0 0.03 0.1 0.3 1 3 10 30
	do
		for I_SEED in 0 # 1 2 3 4
		do
			SIM_NAME="${BASE_NAME}_diffusion_${SCALE_LONG}xLong_${SCALE_LAT}xLat_${I_SEED}"
			PARAM_FILE="params_temp_${SIM_NAME}.yaml"
			cp ${BASE_PARAMS} ${PARAM_FILE}
			yq eval -i ".xlinks.d_i = $(calc ${BASE_DIFF}*${SCALE_LONG})" ${PARAM_FILE}
			yq eval -i ".xlinks.d_side = $(calc ${BASE_DIFF}*${SCALE_LAT})" ${PARAM_FILE}
			yq eval -i ".xlinks.c_bulk = 0.1" ${PARAM_FILE}
			yq eval -i ".motors.c_bulk = 10" ${PARAM_FILE}
			yq eval -i ".t_run = 600" ${PARAM_FILE}
			yq eval -i ".seed = $(( ${BASE_SEED} + ${I_SEED} ))" ${PARAM_FILE}
			# Run sim for these parameter values
			echo "Launching sim ${SIM_NAME} with parameter file ${PARAM_FILE}"
			./cylaks.exe ${PARAM_FILE} ${SIM_NAME} &
		done
	done
done
wait
rm params_temp_${BASE_NAME}_*
