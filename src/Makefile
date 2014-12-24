VPATH=../sofa:.
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

CFLAGS = $(DEBUG) -Wall -std=c99 -I. -I../sofa $(SLAINC) -I/usr/local/include -I/opt/local/include $(HPXINC) -fPIC -fopenmp
LDFLAGS = -L../sofa -lsofa_c $(SLALIB) -L/opt/local/lib $(HPXLIB) -lgomp

SRCS = $(filter-out test_qpoint.c test_math.c,$(wildcard *.c))
OBJS = $(SRCS:.c=.o)
HEADERS = $(wildcard *.h)
HEADER = qpoint.h
LIB = libqpoint.a
PREFIX?=/usr/local

default: all
all: $(LIB)

iers_bulletin_a.c: iers_bulletin_a.dat
	python make_iers_bulletin_a_dot_c.py

$(LIB): $(OBJS)
	ar ru $@ $(OBJS)

install: $(LIB) $(HEADER)
	install -m 644 $(HEADER) $(PREFIX)/include
	install -m 644 $(LIB) $(PREFIX)/lib

uninstall:
	cd $(PREFIX)/include && rm -f $(HEADER)
	cd $(PREFIX)/lib && rm -f $(LIB)

.PHONY: tidy clean

tidy:
	rm -f *~

clean: tidy
	rm -f *.o
	rm -f $(LIB)
