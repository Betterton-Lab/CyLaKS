#!/bin/bash
BASE_NAME="shep"
BASE_PARAMS="params/shepherding.yaml"
echo "Starting ${BASE_NAME} scan"
echo "Base parameter file is ${BASE_PARAMS}"

calc() { awk "BEGIN{ printf \"%.6f\n\", $* }"; }

BASE_SEED=198261346419
BASE_C=10
BASE_HYDRO=75
BASE_KOFF_I_MOT=20

i_var=0
j_var=0
weights=(0.08 0.11 0.20 0.29 0.52 1 1.8 5 12 24 35)

for HYDRO_SCALE in 0.1 0.3 1 3 10
do
	j_var=0
	for BIND_AFF in 0.1 0.3 1 3 10 30
	do
		for I_SEED in 0 # 1 2 3 4
		do
			echo "i, j"
			echo "${i_var}, ${j_var}"
			index=$((${i_var} - ${j_var} + 5))
			echo "${weights[index]}"
			SIM_NAME="${BASE_NAME}_motorMotilityConstNum_${HYDRO_SCALE}xVel_${BIND_AFF}xProc_${I_SEED}"
			PARAM_FILE="params_temp_${SIM_NAME}.yaml"
			cp ${BASE_PARAMS} ${PARAM_FILE}
			yq eval -i ".xlinks.c_bulk = 0.1" ${PARAM_FILE}
			yq eval -i ".motors.c_bulk = $(calc ${BASE_C}*${weights[index]})" ${PARAM_FILE}
			yq eval -i ".motors.k_hydrolyze = $(calc ${BASE_HYDRO}*${HYDRO_SCALE})" ${PARAM_FILE}
			yq eval -i ".motors.k_off_i = $(calc ${BASE_KOFF_I_MOT}/${BIND_AFF})" ${PARAM_FILE}
			yq eval -i ".seed = $(( ${BASE_SEED} + ${I_SEED} ))" ${PARAM_FILE}
			# Run sim for these parameter values
			echo "Launching sim ${SIM_NAME} with parameter file ${PARAM_FILE}"
			./cylaks.exe ${PARAM_FILE} ${SIM_NAME} &
		done
		j_var=$((${j_var} + 1))
	done
	i_var=$((${i_var} + 1))
done
wait
rm params_temp_${BASE_NAME}
