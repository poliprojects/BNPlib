#!/bin/bash
#PBS -S /bin/bash


#-----------------SETTING THE REQUEST FOR THE HARDWARE ALLOCATION----#

#PBS -l nodes=1:ppn=20,walltime=24:00:00 -q gigat

# Set the job name
#PBS -N many_populations

# Set the output file and merge it to the sterr
#PBS -o out-hostname-XyZ-N1x1-qsub.txt
#PBS -j oe
#PBS -e out-hostname-XyZ-N1x1.txt

# Send emails
#PBS -M mario.beraha@polimi.it
#PBS -m abe


#------------------SETTING THE ENVIRONMENT----------------------------#
# Start the job in the current directory (PBS starts in the home folder)
cd ${PBS_O_WORKDIR}

export OMP_NUM_THREADS=20

export mkPrefix=/u/sw
source $mkPrefix/etc/profile
module load gcc-glibc/6
module load openblas

cd scripts/build
cmake ..
make many_population

#-------------------RUN THE EXECUTABLE---------------------------------#

./many_population
