#!/bin/tcsh
#
# Exemple: parallelc3.x  6  99  0.02  -5.0  WDIR
#
# Input parameters
set NP = $1
set NJP = $2
set DW = $3
set WI = $4
set DIRNAME = $5
#
# Inicilize varibles
set PATH = `pwd`
set N = 1
set NI = 1
if (null$NP == null) set NP = 1
if (null$NJP == null) set NJP = 1
if (null$DW == null) set DW = 0
if (null$WI == null) set WI = 0
if (null$DIRNAME == null) set DIRNAME = WDIR
#
# Calculate total number of jobs and number of jobs per processor
@ NT = $NP * $NJP
@ NF = $NJP
#
# Write parameters in screen
echo '   Main directory: '$PATH
echo '   Total number of jobs:' $NT
echo '   Total number of processors:' $NP
echo '   Total number of jobs per processor:' $NJP
#
# Loop of jobs submition
while ( $N <= $NP) 
   mkdir ${DIRNAME}${N}
   cp -rf espec_v04.x ginp2.x calc3.x dipol3.inp $PATH/${DIRNAME}${N} 
   chdir $PATH/${DIRNAME}${N}
   ln -s espec_v04.x espec.x  
#
#generate run.pbs to be submit the jobs in queue
cat > run.pbs <<EOF
#PBS -S /bin/bash
#PBS -N eSPec-Job$N
#PBS -l nodes=1
#PBS -l cput=24:00:00
#PBS -m ae
#
cd \$PBS_O_WORKDIR
#
./calc3.x $DW $NI $NF $WI > Jobs_P${N}.out
EOF
   qsub run.pbs
   chdir $PATH
   echo '      Jobs '$NI' - '$NF' submitted in processor: '$N
   @ N ++
   @ NI = $NF + 1
   @ NF = $NF +  $NJP + 1
end
# End of the program.
echo
echo '   <<<>>> Jobs submitted in parallel! <<<>>> '
echo

