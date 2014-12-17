#!/bin/tcsh
#
# input parameter 
# 	PATH_INPUT = path of initial eigenfunctions
#	NI = initial file
#	NF = final file
#
set DT = $1
set N = $2
set NF = $3
set PATH_INPUT = $4
echo "---------------------------------------------------------------"
echo "      Input parameters directory: $PATH_INPUT "
echo "      Initial job: $N "
echo "      Final job: $NF "
echo "      Time step: $DT"
echo "---------------------------------------------------------------"
echo
#
# Run eSPec and normf programs to obtain the norm (<\psi(\ni)|\psi(\ni)>),
# where |\psi(\ni)> = FFT(\psi(t)>), for different initial conditions. 
#
while($N <= $NF); #do
   echo "Processing job: $N"
#generate input file for eSPec using program ginp.x
cat > fort.1 <<EOF
   ${N}  ${DT}
EOF
   ginp.x > input.spc

   if($N < 10)then
      cat -s $PATH_INPUT/ReIm_000${N}.dat potential.inp > & data.inp
   else if($N < 100)then
      cat -s $PATH_INPUT/ReIm_00${N}.dat potential.inp > & data.inp
   else if($N < 1000)then
      cat -s $PATH_INPUT/ReIm_0${N}.dat potential.inp > & data.inp
   else
      cat -s $PATH_INPUT/ReIm_${N}.dat potential.inp > & data.inp
   endif  
#
   echo 'eSPec time: '
   time espec.x > & xray_${N}.out
   rm xray_${N}.out
#
   if($N < 10)then
      mv ReIm_fwp.dat wp_000${N}.dat
   else if($N < 100)then
      mv ReIm_fwp.dat wp_00${N}.dat
   else if($N < 1000)then
      mv ReIm_fwp.dat wp_0${N}.dat
   else
      mv ReIm_fwp.dat wp_${N}.dat
   endif
#
   @ N++
end
echo 'All jobs processed!'





