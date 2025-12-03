#!/bin/bash
BASE_NAME=$1
BASE_PARAMS=$2 # "k401.yaml"
echo "Starting ${BASE_NAME} scan"
echo "Base parameter file is ${BASE_PARAMS}"

calc() { awk "BEGIN{ printf \"%.6f\n\", $* }"; }

BASE_SEED=198261346419

BASE_C=10

BASE_DIFF=0.131
BASE_KON=0.00238
BASE_KOFF=0.01

BASE_HYDRO=75

BASE_KON_MOT=0.00358
BASE_KOFF_I_MOT=20
BASE_KOFF_II_MOT=260

i_var=0
j_var=0
#weights=(0.2 0.25 0.3 0.4 0.6 1 1.8 5 12 24 35)
weights=(0.08 0.11 0.20 0.29 0.52 1 1.8 5 12 24 35)

#for E_INT in 0.2 0.4 0.6 0.8
for HYDRO_SCALE in 0.1 0.3 1 3 10 30
#for BIND_AFF_SIDE in 1 # 0.0 0.001 0.003 0.01 #0.03 0.1 0.3 1 3 10 30
do
	j_var=0
	for BIND_AFF in 0.1 0.3 1 3 10 30 #0.03 0.1 0.3 1 3 10 30
	do
		for I_SEED in 0 # 1 2 3 4
		do
			echo "i, j"
			echo "${i_var}, ${j_var}"
			index=$((${i_var} - ${j_var} + 5))
			echo "${weights[index]}"

		#	SIM_NAME="${BASE_NAME}_3x_${BIND_AFF_SIDE}x_${I_SEED}"
		#	SIM_NAME="${BASE_NAME}_${E_INT}kT_${BIND_AFF}x_${I_SEED}"
		#	SIM_NAME="${BASE_NAME}_${I_SEED}_motorVelConstNum_${HYDRO_SCALE}x_${BIND_AFF}x"
			SIM_NAME="${BASE_NAME}_${I_SEED}_motorVelWeighted_${HYDRO_SCALE}x_${BIND_AFF}x"
		#	SIM_NAME="${BASE_NAME}_${I_SEED}_xlinkDiffNorm_${BIND_AFF}x_${BIND_AFF_SIDE}x"
			PARAM_FILE="params_temp_${SIM_NAME}.yaml"
			cp ${BASE_PARAMS} ${PARAM_FILE}
		#	yq eval -i ".xlinks.neighb_neighb_energy = ${E_INT}" ${PARAM_FILE}
		#	yq eval -i ".xlinks.d_i = $(calc ${BASE_DIFF}*${BIND_AFF})" ${PARAM_FILE}
		#	yq eval -i ".xlinks.d_side = $(calc ${BASE_DIFF}*${BIND_AFF_SIDE})" ${PARAM_FILE}
		#	yq eval -i ".xlinks.k_off_i = $(calc ${BASE_KOFF}/${BIND_AFF})" ${PARAM_FILE}
		#	yq eval -i ".xlinks.k_on = $(calc ${BASE_KON}*${BIND_AFF})" ${PARAM_FILE}
		#	yq eval -i ".motors.k_on = $(calc ${BASE_KON_MOT}*${BIND_AFF})" ${PARAM_FILE}
			yq eval -i ".motors.c_bulk = $(calc ${BASE_C}*${weights[index]})" ${PARAM_FILE}
			yq eval -i ".motors.k_hydrolyze = $(calc ${BASE_HYDRO}*${HYDRO_SCALE})" ${PARAM_FILE}
			yq eval -i ".motors.k_off_i = $(calc ${BASE_KOFF_I_MOT}/${BIND_AFF})" ${PARAM_FILE}
		#	yq eval -i ".motors.k_off_ii = $(calc ${BASE_KOFF_II_MOT}/${BIND_AFF})" ${PARAM_FILE}
		    	yq eval -i ".seed = $(( ${BASE_SEED} + ${I_SEED} ))" ${PARAM_FILE}
			# Run sim for these parameter values
			echo "Launching sim ${SIM_NAME} with parameter file ${PARAM_FILE}"
			#singularity exec --bind $PWD cylaks_latest.sif cylaks.exe ${PARAM_FILE} ${SIM_NAME} &
			apptainer exec --bind $PWD cylaks_app.sif cylaks.exe ${PARAM_FILE} ${SIM_NAME} &
		done
		j_var=$((${j_var} + 1))
	done
	i_var=$((${i_var} + 1))
done
wait
#rm params_temp_${BASE_NAME}
