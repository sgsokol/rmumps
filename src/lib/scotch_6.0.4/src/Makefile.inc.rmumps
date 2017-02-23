EXE		=
LIB		= .a
OBJ		= .o

#MAKE		= make
AR		= $(AR)
#ARFLAGS		= -shared -o
ARFLAGS		= rv
CAT		= cat
CCS		= $(CC)
#CCP		= mpicc
CCD		= $(CC)
HCFLAGS 	= $(CFLAGS) -O3 -DCOMMON_FILE_COMPRESS_GZ -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -DSCOTCH_PTHREAD -Drestrict=__restrict -DIDXSIZE64
CLIBFLAGS	= -shared -fPIC
#LDFLAGS		= -lz -lm -pthread #-lrt
HLDFLAGS		= $(LDFLAGS)
CP		= cp
LEX		= flex -Pscotchyy -olex.yy.c
LN		= ln
MKDIR		= mkdir
MV		= mv
RANLIB		= echo
YACC		= bison -pscotchyy -y -b y
