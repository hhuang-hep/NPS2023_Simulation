#!/bin/bash

simDir=$(pwd) # simulation home directory
macro=bashrun.mac # macro file to run GEANT4 simulation

gen_type=0 # index for different reactions (0: DVCS (default), 1: exclusive pi0)
# This index should corresponds to the src and include files in the Geant4 package (for DVCS or pi0)

i_job=$4 # index for different output file names for job submission
if [ -z "$i_job" ]; then
    i_job=1
fi

# ------------ Environment setup ------------
CONFIG_FILE="$(dirname "$0")/source_simu.sh"

if [ ! -f "$CONFIG_FILE" ]; then # check if source_simu.sh exist
    echo "Error: Configuration file not found '$CONFIG_FILE'."
    exit 1
fi

if [ -z "$folder_g4" ]; then
    echo "Setting environment for simulation... Sourcing configuration file: '$CONFIG_FILE'."
    source "$CONFIG_FILE"

    if [ -z "$folder_g4" ]; then
        echo "Error: $folder_g4 variable is still unset after sourcing. Please check $CONFIG_FILE"
        exit 2
    else
        echo "Environment setup completed successfully."
    fi
else
    echo "Environment variables found. Configuration load skipped."
fi

# NPS photon reconstruction environment
export NPS_SOFT=$folder_NPS_SOFT
export LD_LIBRARY_PATH=${NPS_SOFT}:${LD_LIBRARY_PATH}
export PATH=${NPS_SOFT}:${PATH}

# Event generator
export DVCS_EVENT_GEN=$folder_DVCS_EVENT_GEN
export LD_LIBRARY_PATH=${DVCS_EVENT_GEN}:${LD_LIBRARY_PATH}
export PATH=${DVCS_EVENT_GEN}:${PATH}

# ------ input parameter of kinematics (based on run numbers) and target for simulation ------
KINEMATICS_FILE="$folder_hms/KinC_list.txt"

if [ ! -f "$KINEMATICS_FILE" ]; then # the file does not exist
    echo "Error: can't find kinematics list file in $KINEMATICS_FILE" >&2
    exit 1
fi

KINEMATICS_LIST=$( # get the list of runs and corresponding kinematics from the list
    awk 'NR > 1 && $1 ~ /^KinC/ {print $1 "\t" $2}' "$KINEMATICS_FILE"
)

if [ -z "$KINEMATICS_LIST" ]; then # the file should not be empty 
    echo "Error: $KINEMATICS_FILE is empty" >&2
    exit 1
fi

# For print out information
declare -A TARGET1_MAP=( ["0"]="LH2" ["1"]="LD2" )
declare -A TARGET2_MAP=( ["0"]="proton" ["1"]="neutron" ["2"]="deuteron" )

# Validation Functions
validate_run_number() { # $1: Run Number
    if [ -z "$1" ]; then return 1; fi # Error: empty
    if echo "$KINEMATICS_LIST" | awk '{print $2}' | grep -q -w "$1"; then
        return 0 # Succeed
    else
        return 2 # Error: not in list
    fi
}

validate_target_1() { # Target 1 (0 or 1)
    if [ -z "$1" ]; then return 1; fi # Error: empty
    if [[ "$1" == "0" || "$1" == "1" ]]; then
        return 0
    else
        return 2 # Error: invalid value
    fi
}

validate_target_2() { # Target 2 (0, 1, or 2)
    if [ -z "$1" ]; then return 1; fi # Error: empty
    if [[ "$1" == "0" || "$1" == "1" || "$1" == "2" ]]; then
        return 0
    else
        return 2 # Error: invalid value
    fi
}

# General input and validation of the parameters
get_validated_input() {
    local param_name="$1" # name of parameters (run number, target type 1, target type 2)
    local param_value="$2" # input parameter value
    local validator="$3" # functions for validation
    local error_prompt="$4" # error message
    local input_var_name="$5" 
    
    # try input parameter first
    if [ -n "$param_value" ]; then
        if $validator "$param_value"; then # valid parameter
            eval "$input_var_name='$param_value'" # dynamic assignment
            return 0
        else # invalid parameter
            VALIDATION_CODE=$?
            echo "Error: $param_name $param_value are invalid. Please re-enter manually" >&2
        fi
    fi

    # Interactive mode if input parameters are invalid
    
    # Print out the list of kinematics
    if [ "$param_name" == "run number" ]; then
        echo "========================================================="
        echo "Valid run numbers are:"
        echo "---------------------------------------------------------"
        echo "$KINEMATICS_LIST" | column -t
        echo "---------------------------------------------------------"
    fi

    while true; do
        read -p "Input $param_name ($error_prompt): " INPUT_VAL
        
        if $validator "$INPUT_VAL"; then
            eval "$input_var_name='$INPUT_VAL'" # dynamic assignment
            break
        else
            VALIDATION_CODE=$?
            local base_prompt="Error: $param_name $INPUT_VAL"
            if [ $VALIDATION_CODE -eq 1 ]; then
                 echo "$base_prompt is empty. Please enter a value." >&2
            elif [ $VALIDATION_CODE -eq 2 ]; then
                 echo "$base_prompt is invalid. $error_prompt" >&2
            fi
        fi
    done
    return 0
}

# 1. Run Number
get_validated_input "run number" "$1" validate_run_number \
    "Please choose a run number from the list" run_number

# 2. Target 1 (target type)
get_validated_input "target1 (Target Type)" "$2" validate_target_1 \
    "0 for hydrogen; 1 for deuterium" target1

# 3. Target 2 (nucleon type)
get_validated_input "target2 (Nucleon Type)" "$3" validate_target_2 \
    "0 for proton; 1 for neutron; 2 for deuteron" target2

# Print out parameters
KIN_NAME=$(echo "$KINEMATICS_LIST" | grep -w "$run_number" | awk '{print $1}')

echo "========================================================="
echo "Run number: $run_number ($KIN_NAME)"
echo "Target 1 (Target type): ${target1} (${TARGET1_MAP[$target1]})"
echo "Target 2 (Nucleon): ${target2} (${TARGET2_MAP[$target2]})"
echo "========================================================="

rm -rf ${simDir}/input_g4/input_${i_job}.txt
file_name=generated_events/${KIN_NAME}_${TARGET1_MAP[$target1]}_${i_job}.root
echo $file_name >> ${simDir}/input_g4/input_${i_job}.txt
echo $i_job >> ${simDir}/input_g4/input_${i_job}.txt
echo $i_job >> ${simDir}/input_g4/input_${i_job}.txt
echo $i_job >> ${simDir}/input_g4/input_${i_job}.txt

# ------------ NPS Geant4 simulation ------------
source /etc/profile.d/modules.sh
module purge
module use /group/halla/modulefiles
module load geant4/11.2.1
module load root/6.30.04

./DVCS $run_number $macro $target1 $target2 < $simDir/input_g4/input_$i_job.txt
echo "Geant4 simulation complete."
echo ""

# ------------ HMS mc-single-arm simulation ------------
if [ $? -eq 0 ]; then
    
    # environment setting
    module purge
    module use /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/modulefiles
    module load root/6.30.04-gcc11.4.0
    export ROOTSYS=/cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/root/6.30.04-gcc11.4.0

    # start the simulation
    echo "Move to $folder_hms and start HMS simulation"
    (
        cd "$folder_hms" || { echo "Error: can't find $folder_hms"; exit 1; }

        # copy the input file so the script can find
        cp infiles/hms_${TARGET1_MAP[$target1]}_${KIN_NAME}.inp infiles/hms_${TARGET1_MAP[$target1]}_${KIN_NAME}_${i_job}.inp
        # print out the setting for recording
        cat infiles/hms_${TARGET1_MAP[$target1]}_${KIN_NAME}_${i_job}.inp
        # run simulation
        ./run_mc_single_arm_tree hms_${TARGET1_MAP[$target1]}_${KIN_NAME}_${i_job} ${KIN_NAME}_${TARGET1_MAP[$target1]}_${i_job}
        # remove the duplicated input file
        rm infiles/hms_${TARGET1_MAP[$target1]}_${KIN_NAME}_${i_job}.inp
    )

    if [ $? -eq 0 ]; then
        echo "HMS simulation complete, move back to $simDir"
    else
        echo "Error: HMS simulation failed."
    fi

else
    echo "Error: Geant4 simulation failed. Skipping HMS simulation."
fi
echo ""

# ------------ vertex z and photon reconstruction ------------
module purge
unset ROOTSYS

module use /group/halla/modulefiles
module load geant4/11.2.1
module load root/6.30.04

cp ./Reconstruction.C ./temp/Reconstruction_${i_job}.C
fileNAME=$simDir/temp/Reconstruction_${i_job}
root -b <<EOF
    .L $fileNAME.C+
        Reconstruction("$KIN_NAME", $target1, $i_job, $gen_type);
    .q
EOF

rm $simDir/temp/Reconstruction_${i_job}.C