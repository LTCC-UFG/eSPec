#!/bin/tcsh
#
# input parameter 
#       DW = frequency step
#	NI = initial file
#	NF = final file
#
set DW = $1
set N = $2
set NF = $3
set WI = $4
set NAUX = 0
echo "---------------------------------------------------------------"
echo "      Initial job: $N "
echo "      Final job: $NF "
echo "      Frequency step: $DW "
echo "---------------------------------------------------------------"
echo
#
# Run eSPec and normf programs to obtain the norm (<\psi(\ni)|\psi(\ni)>),
# where |\psi(\ni)> = FFT(\psi(t)>), for different initial conditions. 
#
while($N <= $NF); #do
   @ NAUX ++
   echo "Processing job: $NAUX ('$N')"
#generate input file for eSPec using program ginp2.x
cat > fort.1 <<EOF
   ${N}  ${DW}  ${WI}
EOF
   ./ginp2.x > input.spc
#
   echo 'eSPec time: '
#
   if($N < 10)then
      time ./espec.x > & ixpp_000${N}.out
   else if($N < 100)then
      time ./espec.x > & ixpp_00${N}.out
   else if($N < 1000)then
      time ./espec.x > & ixpp_0${N}.out
   else
      time ./espec.x > & ixpp_${N}.out
   endif
#
   @ N ++
end
echo 'All jobs processed!'





