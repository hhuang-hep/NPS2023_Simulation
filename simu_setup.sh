#!/bin/bash

# Load configuration to setup soft links for simulation
CONFIG_FILE="$(dirname "$0")/source_simu.sh"

# check if source_simu.sh exist
if [ ! -f "$CONFIG_FILE" ]; then
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

# make folders and soft links
mkdir ./input_g4 # for Geant4 input files
ln -s $folder_g4/DVCS ./DVCS

# link of folders from mc-single-arm
ln -s $folder_hms/generated_events ./generated_events
ln -s $folder_hms/worksim ./worksim

# directory for reconstructed events
mkdir ./rootfiles
mkdir ./temp
