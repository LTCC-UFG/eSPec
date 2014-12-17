#!/bin/sh
#
set DIR = `pwd`
set INPUT = $1
#
if [ -f $INPUT ]; then
   echo
   echo "<<<>>> Input file error <<<>>>"
   echo "type: "
   echo "      runspc.x <input_file>"
   exit 
fi
cp ${INPUT} input.spc
echo "eSPec job number:" 
time nice -0 ./espec.x > & ${INPUT}.out &
sleep 1 
echo "submitted."
echo
sleep 2 
rm input.spc
echo "Output is being written in file:"
echo "      ${DIR}/${INPUT}.out"


