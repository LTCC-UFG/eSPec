#!/bin/sh
#
DIR=`pwd`
MAC=`uname -n`
DATES=`date '+%d%m%Y-%H:%M'`
#
#default values for options
WDIR=$DIR
OPTC=0
OPTT=0
OPTO=0
OPTB=0
OPTBI=0
OPTBO=0
OPTK=0
NICE=0
#
# Define usage message
#
usage (){
	  echo
	  echo "Usage:     runspc.x [OPTIONS] <INPUT_FILE>"
	  echo
	  echo 'Options:'
	  echo ' -c        Compress the output data (no compress by default)'
#	  echo ' -c        Compress the output data'
	  echo ' -t        Run the calculations at tmp (working dir is the default)'
	  echo ' -b        Do not backup old input and output data (backup by default)'
	  echo ' -bi       Do not backup old input data (backup by default)'
	  echo ' -bo       Do not backup old input data (backup by default is)'
	  echo ' -k        Do not clean the working dir'
          echo ' -w dir    Write the compressed output in dir  (working dir is the default)' 
	  echo ' -o file   Rename the output to the new file name '
          echo '              (the input file name with extension .out by default)'
	  echo ' -n nice   Job priority (default 0)'
	 }

while [ -n "`echo $1 | grep '^-'`" ]; do
   case $1 in
      -c ) OPTC=1;;
      -t ) OPTT=1;;
      -b ) OPTB=1;;
      -bi ) OPTBI=1;;
      -bo ) OPTBO=1;;
      -k ) OPTK=1;; 
      -w ) WDIR=$2; shift;;
      -o ) OPTO=1; OUTPUT=$2; shift;;
      -n ) NICE=$2; shift;;
       * ) usage; exit 1;;
   esac
   shift
done
#
# Get input file name
INPUT=$1
#
# Check if the input file exist
if [ ! $INPUT ]; then
   usage
   exit 1 
elif [ ! -f $INPUT ]; then
   echo; echo
   echo "    File \"${INPUT}\" does not exist in the directory \"${DIR}\""
   echo
   usage
   exit 1
fi
#
# Generate the default name for the output file
if [ $OPTO -eq 0 ]; then OUTPUT=${INPUT%.*}.out; fi
#
# Backup the old output files if necessary
if [ $OPTB -eq 0 -a $OPTBO -eq 0 ]; then
   if [ -f $OUTPUT ]; then echo; echo "Backup output file(s):"; fi
   for i in 9  8  7  6  5  4  3  2  1; do
      let j=$i-1
      if [ -f $OUTPUT -a -f "${OUTPUT%.*}-${j}_out.old" ]; then 
         cp -v ${OUTPUT%.*}-${j}_out.old ${OUTPUT%.*}-${i}_out.old; 
      fi
   done
   if [ -f $OUTPUT ]; then cp -v $OUTPUT ${INPUT%.*}-0_out.old; fi
fi
#
# Backup the old input files if necessary
if [ $OPTB -eq 0 -a $OPTBI -eq 0 ]; then
   if [ -f "input.spc" ]; then echo; echo "Backup input file(s):"; fi
   for i in 9  8  7  6  5  4  3  2  1; do
      let j=$i-1
      if [ -f "input.spc" -a -f "input-${j}_spc.old" ]; then 
         cp -v input-${j}_spc.old input-${i}_spc.old; 
      fi
   done
   if [ -f "input.spc" ]; then cp -v input.spc input-0_spc.old; fi
fi
#
# Generate the input file coping the input file to the standard 
# eSPec input file name
if [ "$INPUT" != "input.spc" ]; then cp $INPUT input.spc; fi
#
# Check and get the external input file(s) 
let kf=0
for i in `grep ".FILE" $INPUT`; do let kf=kf+1; done
for i in `grep ".READ" $INPUT`; do let kf=kf+1; done
let ki=0
FILEI=`grep "\.inp" $INPUT; grep "\.dat" $INPUT`
for i in $FILEI; do let ki=ki+1; done
#
if [ "$ki" -ne "$kf" -a "$ki" -ne "0" ]; then
   echo
   echo ' ->  The external input files must end with extensions ".inp" or ".dat"'
   exit 1
fi 
#
# Check, create and copy files to the tmp directory to run eSPec
if [ $OPTT -eq 1 ]; then
   USERNAME=`whoami`
   TDIRN=/tmp/$USERNAME
   if [ ! -d $TDIRN ]; then mkdir $TDIRN; fi
   TDIR=$TDIRN/$INPUT
   if [ ! -d $TDIR ]; then mkdir $TDIR; fi
   if [ "$INPUT" != "input.spc" ]; then cp input.spc $TDIR; fi
   cp $INPUT $FILED $FILEI $TDIR
fi
#
# Output the name of the machine that the job is being submitted
echo
echo "Start running eSPec in" $MAC
echo `date`
#
# Run eSPec
if [ $OPTT -eq 1 ]; then
   cd $TDIR; time nice -$NICE /home/freddy/espec_v0.5/espec.x > $OUTPUT
   if [ $OPTC -eq 0 ]; then 
      cd $TDIR; cp $OUTPUT  $WDIR; 
      if [ -f $TDIR/ReIm_0001.dat  ]; then cp ReIm*  $WDIR; fi
      if [ -f $TDIR/veff_0001.dat  ]; then cp veff*  $WDIR; fi
      if [ -f $TDIR/eigvc_0001.dat ]; then cp eigvc* $WDIR; fi
      if [ -f $TDIR/pulse.dat ]; then cp pulse.dat $WDIR; fi
   fi
else
   time nice -$NICE /home/freddy/espec_v0.5/espec.x > $OUTPUT
fi
echo
echo "Calculations finished!"
echo `date`
#
# Return to the place were the output transcript was written
if [ "$INPUT" != "input.spc" -a $OPTT -eq 0 ]; then rm input.spc; fi
echo;echo
echo "Output files were written in \"${DIR}/${OUTPUT}\""
#
# Compress output data if asked 
if [ $OPTC -eq 1 ]; then
#   DATEF=`date '+%d%m%Y-%l:%M'`
   COUT=${WDIR}/${OUTPUT%.*}_${DATES}.tar.gz
   echo
   echo "Compressing the output data to \"${COUT}\""
   if [ $OPTT -eq 1 ]; then
      cd $TDIR; tar cf - ReIm* veff* eigvc* movie.gplt pulse.dat $FILEI $OUTPUT $INPUT \
	  | gzip > ${COUT}
      if [ $OPTK -eq 0 ]; then if [ -f $COUT -o $OUTPUT ]; then rm -rf $TDIR; fi; fi
   else
      tar cf - ReIm* veff* eigvc* movie.gplt pulse.dat $FILEI $OUTPUT $INPUT \
	  | gzip > ${COUT}
      if [ -f $COUT ]; then rm -f ReIm*dat veff*dat eigvc*dat movie.gplt pulse.dat $OUTPUT; fi 
#	 if [ -f ReIm_0001.dat  ]; then rm -f ReIm*; fi
#	 if [ -f veff_0001.dat  ]; then rm -f veff*; fi
#	 if [ -f eigvc_0001.dat ]; then rm -f eigvc*; fi
#	 if [ -f movie.gplt ]; then rm -f movie.gplt; fi
#	 if [ -f pulse.dat ]; then rm -f pulse.dat; fi
#      fi
   fi
fi
#
exit 0
