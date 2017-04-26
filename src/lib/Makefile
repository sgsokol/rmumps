#DMUMPS=MUMPS_5.0.1
DMUMPS=MUMPS_5.1.1
DMETIS=metis-5.1.0
DSCOTCH=scotch_6.0.4
DLIB=../lib
ARITH = d
FLAGS.MPI=-I../libseq
FLAGS.MUMPS=-Dmetis -Dpord -Dscotch -I../include
FLAGS.PORD=-I../PORD/include
FLAGS.SCOTCH=-I../../$(DSCOTCH)/include
FLAGS.GKlib=-I$(DMETIS)/GKlib -DUSE_GKREGEX
ifeq "$(INC_OS)" ".win"
 FLAGS.GKlib += -DMINGW
endif
FLAGS.metis=-I. -I$(DMETIS)/libmetis -I$(DMETIS)/include -w

#PKG_CXXFLAGS += $(FLAGS.MUMPS) #-DINTSIZE64
#PKG_CFLAGS += -DAdd_ #-DINTSIZE64
#PKG_FFLAGS = $(FLAGS.MUMPS) $(FLAGS.MPI) -w #-fcray-pointer -fdefault-integer-8 
#PKG_FCFLAGS += $(FLAGS.MUMPS) $(FLAGS.MPI) -w #-fcray-pointer -fdefault-integer-8 

SOURCES.C = $(wildcard $(DMETIS)/libmetis/*.c $(DMETIS)/GKlib/*.c)

FUSELIBS=$(DMUMPS)/lib/libdmumps.a $(DMUMPS)/lib/libmumps_common.a $(DMUMPS)/libseq/libmpiseq.a $(DMUMPS)/PORD/lib/libpord.a $(DMETIS)/libmetis/libmetis.a

#$(foreach v, $(.VARIABLES), $(info $(v) = $($(v))))

.PHONY: all allmumps scotchlibs

all: allmumps

scotchlibs:
	$(MAKE) -C $(DSCOTCH)/src CC="$(CC)" PKG_CFLAGS="$(CFLAGS) $(PKG_CFLAGS)" INC_OS="$(INC_OS)" scotch esmumps
	mkdir -p $(DLIB)
	mv -f $(DSCOTCH)/lib/libscotch.a $(DSCOTCH)/lib/libscotcherr.a $(DSCOTCH)/src/esmumps/libesmumps.a $(DLIB)/
metislib: $(SOURCES.C:.c=.o) 
	$(AR) cvrs $(DMETIS)/libmetis/libmetis.a $^

#allmumps: scotchlibs $(FUSELIBS)
#	mkdir -p $(DLIB) && $(AR) cqT $(DLIB)/MUMPS_thin $(FUSELIBS)
#	@(cd $(DLIB) && \
#	 echo -e 'create liballmumps.a\naddlib MUMPS_thin\nsave\nend' | $(AR) -M)
#	$(RM) $(DLIB)/MUMPS_thin $^
#	$(RANLIB) $(DLIB)/liballmumps.a
#	echo $(CURDIR)
allmumps: scotchlibs metislib
	@(mkdir -p $(DMUMPS)/lib)
	#$(MAKE) -C $(DMUMPS)/libseq CC="$(CC)" FC="$(FC)" PKG_CFLAGS="$(PKG_CFLAGS)" PKG_FCFLAGS="$(PKG_FCFLAGS)" AR="$(AR) cvr " RANLIB="$(RANLIB)"
	$(MAKE) -C $(DMUMPS) mumps_lib ARITH=$(ARITH) CC="$(CC)" FC="$(FC)" PKG_CFLAGS="$(PKG_CFLAGS) -DMUMPS_ARITH=MUMPS_ARITH_$(ARITH) $(FLAGS.MUMPS) $(FLAGS.MPI) $(FLAGS.PORD) $(FLAGS.SCOTCH)" PKG_FCFLAGS="$(PKG_FCFLAGS) $(FLAGS.MUMPS) $(FLAGS.MPI)" AR="$(AR) cvr " RANLIB="$(RANLIB)"
	#$(AR) cqT $(DLIB)/MUMPS_thin $(FUSELIBS)
	#@(cd $(DLIB) && \
	# echo -e 'create liballmumps.a\naddlib MUMPS_thin\nsave\nend' | $(AR) -M)
	cp $(FUSELIBS) $(DLIB)/
	#$(RANLIB) $(DLIB)/liballmumps.a

clean:
	$(RM) $(DMUMPS)/PORD/lib/*.[oa] $(DMUMPS)/libseq/*.[oa] $(DMETIS)/GKlib/*.[oa] $(DMETIS)/libmetis/*.[oa]
	@(cd $(DSCOTCH)/src && $(MAKE) clean)
	@(cd $(DMUMPS)/src && $(MAKE) clean)
#	@(cd $(DMUMPS)/src && $(MAKE) -f Makevars.win distclean)
distclean: clean
	$(RM) lib*/lib*.*
.c.o:
	$(CC) $(PKG_CFLAGS) $(FLAGS.GKlib) $(FLAGS.metis) -c $< -o $@
