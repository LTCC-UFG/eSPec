#!/bin/tcsh
#
# Exemple: parallel.x -0.73250 4 4096 Pdata1
# Input parameters
set DT = $1
set NP = $2
set NR = $3
set PATH_INDAT = $4
#
# Inicilize varibles
set PATH = `pwd`
set N = 1
set NI = 1
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
# Calculate neber of jobs per processor
@ NJP = $NR / $NP
echo '   Total number of jobs per processors:' $NJP
#
# Initiate NF varible
@ NF = $NJP
#
# Loop of jobs submitions
while ( $N <= $NP) 
   mkdir FDIR${N}
   cp -rf potential.inp calc.x ginp.x $PATH/FDIR${N} 
   chdir $PATH/FDIR$N
   ln -s $PATH/espec.x espec.x
   calc.x $DT $NI $NF $PATH_IN > & Jobs_P${N} &
   chdir $PATH
   echo '      Jobs '$NI'-'$NF' submitted in processor: '$N
   @ NI = $NI + $NJP
   @ NF = $NF + $NJP
   @ N ++
end
# End of the program.
echo
echo '   <<<>>> Jobs submitted in parallel! <<<>>> '
echo

