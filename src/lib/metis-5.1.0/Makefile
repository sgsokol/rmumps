#R vars
#PKG_CFLAGS=-O3 -mtune=native -fpic
#AR=ar


DMETIS=.
FLAGS.GKlib=-I$(DMETIS)/GKlib -DUSE_GKREGEX
FLAGS.metis=-I. -I$(DMETIS)/libmetis -I$(DMETIS)/include -w
SOURCES.C = $(wildcard $(DMETIS)/libmetis/*.c $(DMETIS)/GKlib/*.c)

all: metislib

metislib: $(SOURCES.C:.c=.o) 
	$(AR) cvrs $(DMETIS)/libmetis/libmetis.a $^
.c.o:
	$(CC) $(PKG_CFLAGS) $(FLAGS.GKlib) $(FLAGS.metis) -c $< -o $@
