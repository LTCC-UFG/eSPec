#!/bin/tcsh
#
# Exemple: parallelc.x -0.73250 4 4096 Pdata2
#
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
set PATH_INPUT = $PATH_IN
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
# Loop of jobs submition
while ( $N <= $NP) 
   mkdir FDIR${N}
   cp -rf potential.inp ginp.x $PATH/FDIR${N} 
   chdir $PATH/FDIR$N
   ln -s $PATH/espec.x espec.x
#
#generate calc.pbs to be submit the jobs in queue
cat > calc.pbs <<EOF
#!/bin/tcsh
#PBS -S /bin/bash
#PBS -N eSPec-Job$N
#PBS -l cput=720:05:00
#PBS -m ae
#
cd \$PBS_O_WORKDIR
#
# input parameter
#       PATH_INPUT = path of initial eigenfunctions
#       NI = initial file
#       NF = final file
#
set DT = $DT
set N = $NI
set NF = $NF
set PATH_INPUT = $PATH_INPUT
set MAC = \`hostname\`
echo "---------------------------------------------------------------"
echo "      Machine name: \$MAC "
echo "      Input parameters directory: \$PATH_INPUT "
echo "      Initial job: \$N "
echo "      Final job: \$NF "
echo "      Time step: \$DT "
echo "---------------------------------------------------------------"
#
# Run eSPec and normf programs to obtain the norm (<\psi(\ni)|\psi(\ni)>),
# where |\psi(\ni)> = FFT(\psi(t)>), for different initial conditions.
#
while(\$N <= \$NF); #do
   echo "Processing job: \$N"
#generate input file for eSPec using program ginp.x
cat > fort.1 <<EOFI
   \${N}  \${DT}
EOFI
   ginp.x > input.spc
#
   if(\$N < 10)then
      cat -s $PATH_INPUT/ReIm_000\${N}.dat potential.inp > & data.inp
   else if(\$N < 100)then
      cat -s $PATH_INPUT/ReIm_00\${N}.dat potential.inp > & data.inp
   else if(\$N < 1000)then
      cat -s $PATH_INPUT/ReIm_0\${N}.dat potential.inp > & data.inp
   else
      cat -s $PATH_INPUT/ReIm_\${N}.dat potential.inp > & data.inp
   endif
#
   echo 'eSPec time: '
   time espec.x > & xray_\${N}.out
   rm xray_\${N}.out
#
   if(\$N < 10)then
      mv ReIm_fwp.dat wp_000\${N}.dat
   else if(\$N < 100)then
      mv ReIm_fwp.dat wp_00\${N}.dat
   else if(\$N < 1000)then
      mv ReIm_fwp.dat wp_0\${N}.dat
   else
      mv ReIm_fwp.dat wp_\${N}.dat
   endif
#
   @ N++
end
echo 'All jobs processed!'

EOF
   qsub calc.pbs
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

