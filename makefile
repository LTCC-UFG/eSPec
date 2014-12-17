VERSION = 06
TARGET = Linux-i686
#
PROG = espec_v$(VERSION).x
# Choose system 
# Linux: 
#        -malign-double uses 64 bit data alignment. This flag does not 
#                       work in Opteron machines, since this machines 
#                       are already 64 bits.
#   FFLAGS = -O3 -i4 -ident -Wall -static -malign-double
   FFLAGS = -O3 -Wall -std=gnu -fdefault-double-8 -fdefault-real-8 -static -march=native -ffloat-store -ffixed-line-length-none -fd-lines-as-comments
   FC = gfortran
   FFC = gfortran
# Digital 
#   FFLAGS = -O0 -i4 -ident -o $@ 
#   FC = f77 
# Other 
#   FFLAGS = 
#   FC = f77 
#
START    = $(shell date)
WORKDIR  = $(shell pwd)
DATE     = $(shell date +%Y-%m-%d)
#
# libraries
#
XLIBS = -L/usr/X11R6/lib -lX11 
PGPLOT_LIB = -L/usr/local/lib/pgplot -lpgplot
BLAS_LIB = -L/usr/lib/sse -lblas -L/usr/lib/sse/atlas -lblas
#
#
FILES = makefile input.example
#
#
MAIN_OBJ = src/espec.o 
#
ESPEC_OBJS = src/rdinput.o src/compar.o src/chlength.o src/rdpt.o \
	src/gvpot.o src/potfuncs.o src/mtrxdiag.o src/lanczs.o src/lnz.o \
	src/reort.o src/ecnorm.o src/au.o src/au1d.o src/au2di.o \
	src/au2dit.o src/au2dil.o \
        src/au2dc.o src/au3d.o src/rans.o src/lanczsc.o src/rdpte.o \
	src/dfdxi.o src/rif.o src/sod.o src/psod.o src/sil.o src/plnz.o \
	src/fcorrel.o src/eigenerg.o src/spectrumtd.o src/spectrumti.o \
	src/makecst.o src/getcorr.o src/prtpot.o src/prteigvc.o src/ef.o \
	src/chpot.o src/ppsod.o src/pplnz.o src/hwhm.o src/dipl.o src/prpt2.o \
	src/rdpt2.o src/sodfft.o src/psodfft.o src/ddmtv.o src/ddmtvtn.o \
	src/dvtn.o src/spofft.o src/pspofft.o src/dkeft.o src/eigenergfft.o \
	src/eigenergfft1.o src/absorb.o src/s2ppsod.o src/abm1.o \
	src/s2ppabm.o src/sod1.o src/s2ppabm2.o src/init_cond.o src/gauss.o 
#
XESPEC_OBJS =	
#
NESPEC_OBJS =
#
LAPACK_OBJS = src/dlagtf.o src/dlartg.o src/dlassq.o src/dlaebz.o \
	src/dlasrt.o src/dlarnv.o src/dlagts.o src/dlaruv.o src/dstevx.o \
	src/dlamch.o src/dlapy2.o src/dlasr.o src/dlascl.o src/dlae2.o \
	src/dlanst.o src/dlaset.o src/dlaev2.o src/ieeeck.o src/ilaenv.o \
	src/dstebz.o src/dsterf.o src/dstein.o src/dsteqr.o src/dspevx.o \
	src/dlansp.o src/dsptrd.o src/dopgtr.o src/dopmtr.o src/dlarfg.o \
	src/dorg2l.o src/dorg2r.o src/dlarf.o 
#
FFTS_OBJS = src/cfft.o src/fft3d.o 
#
BLAS_OBJS = src/lsame.o src/dasum.o src/dscal.o src/dswap.o src/daxpy.o \
	src/ddot.o src/dcopy.o src/idamax.o src/xerbla.o src/dnrm2.o \
	src/dgemv.o src/dspmv.o src/dspr2.o src/dger.o 
# 
OBJS = $(MAIN_OBJ) $(ESPEC_OBJS) 
#
all:	
	echo 'Compiling program eSPec ...'
	make espec
	echo 'Program eSPec compiled!'

	echo

	echo 'Compiling program normf ...'
	make normf
	echo 'Program normf compiled!'

	echo

	echo 'Compiling program normfa ...'
	make normfa
	echo 'Program normfa compiled!'	

	echo

	echo 'Compiling program normfb ...'
	make normfb
	echo 'Program normfb compiled!'

	echo

	echo 'Compiling program fluxp ...'
	make fluxp
	echo 'Program fluxp compiled!'

#
espec:  $(OBJS)
	make blas
	make ffts
	make lapack
	$(FFC) $(FFLAGS) -o $(PROG) $(OBJS) $(BLAS_OBJS) $(LAPACK_OBJS) $(FFTS_OBJS)
	size $(PROG)
	ln -sf $(PROG) espec.x
	chmod a+xr $(PROG) espec.x
#
normf:  src/normfs.o src/ecnorm.o src/spline.o src/chlength.o $(BLAS_OBJS)
	make ffts
	$(FC) $(FFLAGS) -o normf.x src/normfs.o src/cfft.o src/ecnorm.o \
		src/spline.o src/chlength.o src/dlamch.o $(BLAS_OBJS)
	size normf.x
	chmod a+xr normf.x
#
normfa:  src/normfsa.o src/ecnorm.o src/spline.o src/chlength.o $(BLAS_OBJS)
	make ffts
	$(FC) $(FFLAGS) -o normfa.x src/normfsa.o src/cfft.o src/ecnorm.o \
		src/spline.o src/chlength.o src/dlamch.o  $(BLAS_OBJS)
	size normfa.x
	chmod a+xr normfa.x
#
normfb:  src/normfsb.o src/ecnorm.o src/spline.o src/chlength.o $(BLAS_OBJS)
	make ffts
	$(FC) $(FFLAGS) -o normfb.x src/normfsb.o src/cfft.o src/ecnorm.o \
		src/spline.o src/chlength.o src/dlamch.o  $(BLAS_OBJS)
	size normfb.x
	chmod a+xr normfb.x
#
fluxp:  src/flux.o src/splinem2.o src/rdpt2.o src/mtrxdiag.o src/prteigvc.o \
		src/chlength.o $(BLAS_OBJS) $(LAPACK_OBJS)
	$(FC) $(FFLAGS) -o fluxp.x src/flux.o src/splinem2.o src/rdpt2.o \
		src/chlength.o src/mtrxdiag.o src/prteigvc.o $(BLAS_OBJS) \
		$(LAPACK_OBJS)
	size fluxp.x
	chmod a+xr fluxp.x
#
blas:   $(BLAS_OBJS)
#
ffts:   $(FFTS_OBJS)
#
lapack: $(LAPACK_OBJS) 
#
clean:	 
	rm -f $(OBJS) $(XESPEC_OBJS) $(NESPEC_OBJS) $(BLAS_OBJS) $(LAPACK_OBJS) $(FFTS_OBJS) \
	core a.out src/normfs.o src/ecnorm.o src/spline.o
#
veryclean:  
	make clean
	rm -f $(PROG) src/*.o core a.out normf.x normfa.x normfb.x 
#
distclean:
	make veryclean 
	rm -rf  *.* *~ core "*#" "#*" "*#" "#*" tes* \
		src/cab.f src/*~ src/*.log src/*.dat src/*.ptc src/*.spc \
		src/*.xmgr src/*.x src/tes* src/LIXO src/Lixo src/lixo \
		src/core src/a.out doc/*~ doc/core src/"*#" src/"#*" \
		doc/"*#" doc/"#*" 
#
distrib:
	rm -rf /tmp/eSPec
	mkdir /tmp/eSPec
	cp -rf src tools\&scrips ttests doc makefile /tmp/eSPec
	cd /tmp/eSPec ; \
	   rm -f /tmp/eSPec/*,v /tmp/eSPec/src/*,v \
		/tmp/eSPec/tools/*,v /tmp/eSPec/doc/*,v /tmp/eSPec/*.modif \
		/tmp/eSPec/src/*.modif /tmp/eSPec/tools/*.modif \
		/tmp/eSPec/doc/*.modif ; \
	   rm -rf  /tmp/eSPec/RCS /tmp/eSPec/src/RCS ; \
	   make distclean ; \
	   make ; \
	   make distclean 
	cp -rf runspc.x runtests.x /tmp/eSPec
	cd /tmp ;\
	   tar cf - eSPec | gzip > $(WORKDIR)/eSPec_v$(VERSION).tar.gz 
	rm -rf /tmp/eSPec

