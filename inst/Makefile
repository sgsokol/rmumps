
topdir = MUMPS_5.0.1
libdir = $(topdir)/lib
metisdir = metis-5.1.0

include $(topdir)/Makefile.inc


#SRC = mumps_cmplx_p.f90  mumps_p.f90 call_mumps_cmplx_p.f90 call_mumps_p.f90  
LIB = -L$(topdir)/lib -ldmumps -lmumps_common -lpord  -L$(LMETIS)  $(LIBSEQ)  $(LIBBLAS)


all:
	($(MAKE) -C $(metisdir) config )
	($(MAKE) -C $(metisdir) )
	($(MAKE) -C $(topdir) d)
	#$(FL) $(OPTL)  $(SRC) -o libMUMPS.so -I$(topdir)/include $(LIB)
clean:
	(cd $(metisdir); $(MAKE) distclean)
	(cd $(topdir); $(MAKE) clean)
	rm -f libMUMPS$(LIBEXT)
