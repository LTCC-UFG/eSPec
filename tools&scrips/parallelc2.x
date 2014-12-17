#!/bin/tcsh
#
# Exemple: parallelc2.x -0.73250 4 4096 Pdata2
#
# Input parameters
set DT = $1
set NP = $2
set NR = $3
set PATH_INDAT = $4
set DIRNAME = P1DIR
#
# Inicilize varibles
set PATH = `pwd`
set N = 1
if (null$NP == null) set NP = 1
if (null$NR == null) set NR = 1
if (null$PATH_INDAT == null) set PATH_INDAT = Pdata
set PATH_IN = $PATH/$PATH_INDAT
#
# Write parameters in screen
echo '   Main directory: '$PATH
echo '   Total number of jobs:' $NR
echo '   Total number of processors:' $NP
#
# Calculate number of jobs per processor
@ NJP = $NR / $NP
@ NFA = $NR - ($NP * $NJP)
if ($NFA == 0) set NFA = $NP
@ NF = $NR - $NFA - 1 
echo '   Total number of jobs per processors:' $NJP $NF $NFA
#
# Loop of jobs submition
while ( $N <= $NP) 
   mkdir ${DIRNAME}${N}
   cp -rf potential.inp ginp.x calc2.x $PATH/${DIRNAME}${N} 
   chdir $PATH/${DIRNAME}${N}
   ln -s $PATH/espec.x espec.x
#
#generate run.pbs to be submit the jobs in queue
cat > run.pbs <<EOF
#PBS -S /bin/bash
#PBS -N eSPec-Job$N
#PBS -l cput=720:05:00
#PBS -m ae
#
cd \$PBS_O_WORKDIR
#
./calc2.x $DT $N $NF $NP $PATH_IN > Jobs_P${N}.out
EOF
   @ NS = $N + $NP
   qsub run.pbs
   chdir $PATH
   echo '      Jobs '$N', '$NS', ..., '$NF' submitted in processor: '$N
   @ N ++
   @ NF ++
   if ($NF > $NR) @ NF = $NF - $NP  
end
# End of the program.
echo
echo '   <<<>>> Jobs submitted in parallel! <<<>>> '
echo

