# NPS2023_Simulation
## Introduction
This is the simulation scripts for 2023-2024 NPS experiment.\
It is based on 4 packages with different environment configuration requirements:

- Event generator
- NPS-Geant4 package 
- NPS photon reconstruction software

**Requires:** root-6.30.04 and geant4-11.2.1 in /group/halla/modulefiles

- HMS MC-single-arm package

**Requires:** root-6.30.04 and gcc11.4.0 in /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/modulefiles

After downloading and compiling these packages, one can easily setup and run the simulation with the script in this repository

## Download and compiling of required packages
1. Generator, Geant4 and NPS software
    - Setup environment\
    `module purge`\
    `module use /group/halla/modulefiles`\
    `module load geant4/11.2.1`\
    `module load root/6.30.04`
    - Compile event generator\
    `git clone https://github.com/hhuang-hep/dvcs_gen.git`\
    `cd dvcs_gen`\
    `make lib`\
    `cd ..`

    - Compile NPS software\
    `git clone https://github.com/hhuang-hep/NPS_SOFT.git`\
    `cd NPS_SOFT`\
    `make biglib`\
    `cd ..`

    - Compile NPS-Geant4 package\
    `git clone https://github.com/hhuang-hep/HallC_NPS.git`\
    `cd HallC_NPS/DVCS_evt_gen/`\
    `mkdir build`\
    `cd build`\
    `ccmake ../DVCS/` (then type 'c' to configure and 'g' to generate the 'Makefile')\
    `make`
2. HMS MC-single-arm package
    - Setup environment\
    `module purge`\
    `module use /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/modulefiles`\
    `module load root/6.30.04-gcc11.4.0`\
    `setenv ROOTSYS /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/root/6.30.04-gcc11.4.0`

    - Compile MC-single-arm package\
    (method from Yaopeng: https://github.com/YaopengZhang/mc-single-arm/blob/NPS/README.md) \
    `git clone --branch NPS https://github.com/YaopengZhang/mc-single-arm.git`\
    `cd mc-single-arm/`\
    `cd src`\
    `make`\
    `cd ../util/ntuple/`\
    `make`\
    `cd ../root_tree/`\
    `make`

## Setup for simulation
1. Download this repository\
    `git clone https://github.com/hhuang-hep/NPS2023_Simulation.git`

2. Modify the configuration file "source_simu.sh"
    - folder_g4: path to the "build" directory of NPS-Geant4 package, which contains the executable 'DVCS'
    - folder_hms: path to the mc-single-arm directory
    - folder_NPS_SOFT: path to the NPS_SOFT directory
    - folder_DVCS_EVENT_GEN: path to the dvcs_gen directory

3. Create folders and soft links for running simulation\
    `./simu_setup.sh`\
    this should create the folders "input_g4", "rootfiles", and "temp",\
    as well as the soft links of "generated_events", "worksim" from same forlder in mc-single-arm package, and the executable 'DVCS' from NPS-Geant4

## Run a simulation
If the path in "source_simu.sh" are setup correctly, the simulation can be performed by "NPS_HMS_simu.sh"\
This script automates the environment setup and execution of the packages above\
It requires 3 (or 4) parameters:
1. Run number of the corresponds kinematics (e.g., 3728 for simulating KinC_x36_5_3)
2. Target type (0 for liquid hydrogen, 1 for liquid deuterium)
3. Nucleon type (0 for proton, 1 for neutron, 2 for deuteron)
4. The 4th parameter is set to 1 by default and not always required to input. It plays two roles:\
    (1) a suffix for the name of output files, so that they won't be overwritten when running with jobs\
    (2) a random seed for NPS-Geant4 package (DO NOT set 0 to this parameter. It's not going to work)

To run a simulation, one can simply do\
`./NPS_HMS_simu.sh`

This will list all the kinematics and their run numbers, then ask to input from one of them.\
After input the 2nd and 3rd parameters, the simulation will start.

Alternatively, one can also do the following if the parameters are known\
`./NPS_HMS_simu.sh <run number> <target type> <nucleon type> <job index(optional)>`

The output files can be found in the following folders after the simulation:
- generated_events: files from Geant4
- worksim: files from HMS mc-single-arm
- rootfiles: files of the final results of simulation, including the reconstructed vertex z and clusters/photons in NPS

The function of the rest of files are:
- bashrun.mac: the macro file for Geant4, including the number of events to simulate.
- MakeRunMacro.sh: a script to easily change the number of events in bashrun.mac by `./MakeRunMacro.sh <number of events>`
- Reconstruction.C: reconstruction of vertex z and photons in NPS and output the final results of simulation\
The details of the simulation and output variables are in https://indico.jlab.org/event/946/contributions/16514/attachments/12609/20085/20250506_DVCS_simulation_Hao_Huang.pdf

## Generate exclusive pi0 events
1. Go to 'HallC_NPS/DVCS_evt_gen/DVCS'. There is a folder "pi0Gen" contains the source and header files for the simulation of exclusive pi0 events (in 'pi0Gen/src' and 'pi0Gen/include').
2. Replace the files in 'HallC_NPS/DVCS_evt_gen/DVCS/src' and 'HallC_NPS/DVCS_evt_gen/DVCS/include' with the files in 'pi0Gen' (You could make a backup before this step).
3. Recompile the package.
4. Modify 'gen_type' in NPS_HMS_simu.sh (please set to 1)
5. Run the simulation

## Potential issue when setting or running the simulation
1. Failed to execute 'simu_setup.sh'
    - solution: check the permition and make it executable using 'chmod'

2. Conflict error when loading root/6.30.04 if other modules are loaded (e.g., with 'module use /group/nps/modulefiles' in ~/.cshrc)
    - solution: comment out that line in ~/.cshrc and start the simulation in a clean shell
    - potential solution: add 'module 'unuse /group/nps/modulefiles' after every 'module purge' in NPS_HMS_simu.sh

3. The option 'g' doesn't show up after 'ccmake ../DVCS/' and 'c' when installing the Geant4 package.
    - solution: after 'ccmake ../DVCS/' and 'c', do another 'c'. Now you should see 'g' at the bottom and could generate the Makefile. 


