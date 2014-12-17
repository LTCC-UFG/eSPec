#!/bin/sh
#
# This script run a set of tests for eSPec program comparing the corrent 
# computed results with previus calculated results. In the cases where 
# experimetal results were found, the comparison between the experimental 
# and computed values are in agreement.
# 
#
# Tests description:
# ==================
#    	test1:  Generation of estationary wave functions (WF), in a morse 
#               empirical potential. Propagate the ground state WF with PSOD 
#               method for 1 dimensinal case also in a morse empirical potential. 
#               The spectra profile, considereing a SG window, is computated. 
#		Some of outputs are: print potential and print eingenvectors. 
#               (The spectrum show good agreement with the experimental results)
#		
#       test2:  Test the generation of estationary wave functions, 
#		the PPSOD propagation method (considering 1 dimensinal case)
#               without laser pulse and the diople moment, the OHS empirical 
#               potential, the computation of spectra profile (considereing 
#               the life time of the system HWHM), print potential, and print 
#               eingenvectors. (The spectrum show good agreement with the 
#               experimental results related with the line possitions.
#               The intensities are wrong because of the harmonic potential 
#               used in calculations)
#
#       test3:  Test the generation of estationary wave functions, 
#		the PPSOD propagation method (considering 1 dimensinal case) 
#               with laser pulse and the diople moment (the dipole moment is
#               read from a external file) and the morse empirical potential.
#
#       test4:  Read the wave function and potential from a external file. 
#		Propagate WF with PSIL method for 1 dimensinal case. 
#               Compute the spectra profile.
#               (The spectrum show good agreement with the experimental results)
#
#       test5:  Read the wave function and potential from a external file. 
#		Propagate WF with PFSOD method for 1 dimensinal case. 
#               Compute the spectra profile.
#               (The spectrum show good agreement with the experimental results)
#
#       test6:  Read the wave function and potential from a external file. 
#		Propagate WF with PFSPO method for 1 dimensinal case. 
#               Compute the spectra profile.
#               (The spectrum show good agreement with the experimental results)
#
#       test7:  Read the 2-dimentional potentials from a external file. 
#		Calculate the initial WFs with LANCZSG procedure. 
#               The propagation is perfomed by PSOD method and a Hamiltoniana 
#               within the cross term (Dimension = .2DC). 
#
#       test8:  Read the correlation function from a external file, considering 
#               the life time of the system (applying the HWHM) and compute the 
#               spectra doing the Fast Fourier Transform (FFT).
#	
#    test9-14:  Test respectively the propagation methods: .PSOD, .PFSOD, .PPSOD, 
#               .PFSPO, .PSIL, and .PPSIL. For the 1-dimensional case. The computed 
#               spectrum must be the same for all propagation methods.
#       
#      
#
#
### tests
DIR=`pwd`
for i in 1  2  3  4  5  6  7  8  9  10  11  12  13  14
#
do
#
   if [ ! -e dipol2.inp -a "${i}" -eq "3" ]; then      
      ln -s ttests/dipol2.inp dipol2.inp
   elif [ ! -e eigpot_test3.inp -a "${i}" -ge "4" -a "${i}" -le "6" ]; then 
      ln -s ttests/eigpot_test3.inp eigpot_test3.inp
   elif [ ! -e water_2.pot -a "${i}" -eq "7" ]; then
      ln -s ttests/water_2.pot water_2.pot 
   elif [ ! -e correl8.inp -a "${i}" -eq "8" ]; then
      ln -s ttests/correl8.inp correl8.inp 
   fi
#
   echo
   echo 'Performing eSPec test #'${i}
#
   if [ -e input.spc ];then rm -f input.spc; fi
   cp -f ttests/input.test${i} input.spc
#
   if [ -e compar.test${i} ];then rm -f compar.test${i}; fi
   time ./espec.x > compar.test${i}
#
   if [ -e diffs_test${i} ];then rm -f diffs_test${i}; fi
   diff compar.test${i} ttests/results.test${i} > diffs_test${i}
#   
   if [ -s diffs_test${i} ]; then
      echo
      echo '   Some differeces were found in test #'${i}
      echo
      echo '   The differences between the corrent calculation and the default '
      echo '   result are outputed in the file:'
      echo '      '${DIR}/diffs_test${i}
      echo '   The results of this test are written in the file:'
      echo '      '${DIR}/compar.test${i}
   else 
      echo
      echo '   Test #'${i}' finished without any erros!'
      rm -f diffs_test${i}
      rm -f compar.test${i}
   fi
#
   if [ "${i}" -eq "3" ]; then
      rm -f dipol2.inp 
   elif [ "${i}" -ge "4" -a "${i}" -le "6" ]; then 
      rm -f eigpot_test3.inp 
   elif [ "${i}" -eq "7" ]; then
      rm -f water_2.pot 
   elif [ "${i}" -eq "8" ]; then
      rm -f correl8.inp 
   fi
done
#
rm input.spc
#
echo
echo 'eSPec tests finished!'


