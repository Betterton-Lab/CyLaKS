#!/bin/bash
BASE_NAME="shep"
BASE_PARAMS="params/shepherding.yaml"
echo "Starting ${BASE_NAME} scan"
echo "Base parameter file is ${BASE_PARAMS}"

calc() { awk "BEGIN{ printf \"%.6f\n\", $* }"; }

BASE_SEED=198261346419
ENERGIES=(0.6 0.6 0.6 0.8 1.2)
for MOT_CONC in 50 10
do
	I_ENERGY=1
	for N_PFS in 8 5 3 2 1
	do
		E_INT=${ENERGIES[I_ENERGY]}
		for I_SEED in 0 # 1 2 3 4
		do
			SIM_NAME="${BASE_NAME}_protoNum_${MOT_CONC}nM_${N_PFS}_${I_SEED}"
			PARAM_FILE="params_temp_${SIM_NAME}.yaml"
			cp ${BASE_PARAMS} ${PARAM_FILE}
			yq eval -i ".filaments.n_subfilaments = ${N_PFS}" ${PARAM_FILE}
			yq eval -i ".xlinks.neighb_neighb_energy = ${E_INT}" ${PARAM_FILE}
			yq eval -i ".motors.c_bulk = ${MOT_CONC}" ${PARAM_FILE}
			yq eval -i ".xlinks.c_bulk = 0.1" ${PARAM_FILE}
			yq eval -i ".t_run = 600" ${PARAM_FILE}
			yq eval -i ".seed = $(( ${BASE_SEED} + ${I_SEED} ))" ${PARAM_FILE}
			# Run sim for these parameter values
			echo "Launching sim ${SIM_NAME} with parameter file ${PARAM_FILE}"
			./cylaks.exe ${PARAM_FILE} ${SIM_NAME} &
		done
	I_ENERGY=$(( $I_ENERGY + 1 ))
	done
done
wait
rm params_temp_${BASE_NAME}_*
