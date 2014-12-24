DEBUG = -g -O3
# DEBUG = -O3
CC = gcc

# HACK
# need gfortran for slarefro package
ifeq ($(shell which gfortran), )
SLAINC =
SLALIB =
else
VPATH = ../sofa:../slarefro:.
SLAINC = -I../slarefro -DSLAREFRO
SLALIB = -L../slarefro -lslarefro
endif

ifeq ($(HEALPIX), )
HPXINC =
HPXLIB = -lchealpix
else
HPXINC = -I$(HEALPIX)/include
HPXLIB = -L$(HEALPIX)/lib -lchealpix
endif

CFLAGS = $(DEBUG) -Wall -std=c99 -I../src -I../sofa $(SLAINC) -I/usr/local/include -I/opt/local/include $(HPXINC) -fPIC -fopenmp
LDFLAGS = -L../sofa -lsofa_c $(SLALIB) -L/opt/local/lib $(HPXLIB) -lgomp

default: all
all: test

test: test_qpoint test_math

test_qpoint: test_qpoint.o $(LIB)
	gcc $(CFLAGS) -o $@ $< $(LDFLAGS) -lgetdata -L../src -lqpoint

test_math: test_math.o ../src/sincos.o
	gcc $(DEBUG) -o $@ $< ../src/sincos.o

.PHONY: tidy clean

tidy:
	rm -f *~

clean: tidy
	rm -f *.o
