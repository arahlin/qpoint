VPATH=../sofa:.
DEBUG = -g
# DEBUG =
CC = gcc

# HACK
# need gfortran for slarefro package
ifneq ($(shell which gfortran), )
ifeq ($(DISABLE_SLAREFRO), )
ifeq ($(ENABLE_LITE), )
VPATH = ../sofa:../slarefro:.
SLAINC = -I../slarefro -DENABLE_SLAREFRO
SLALIB = -L../slarefro -lslarefro
endif
endif
endif

ifeq ($(HEALPIX), )
HPXINC =
HPXLIB = -lchealpix
else
HPXINC = -I$(HEALPIX)/include
HPXLIB = -L$(HEALPIX)/lib -lchealpix
endif

ifeq ($(EXTERNAL_CFITSIO), yes)
ifeq ($(CFITSIO_EXT_PREFIX), )
CFITSINC = -I$(CFITSIO_EXT_INC)
CFITSLIB = -L$(dir $(CFITSIO_EXT_LIB)) -l$(notdir $(CFITSIO_EXT_LIB))
else
CFITSINC = -I$(CFITSIO_EXT_PREFIX)/include
CFITSLIB = -L$(CFITSIO_EXT_PREFIX)/lib -lcfitsio
endif
else
CFITSINC =
CFITSLIB = -lcfitsio
endif

LIBNAME_FULL = libqpoint
LIBNAME_LITE = libqpntlt

ifneq ($(ENABLE_LITE), )
CFITSINC =
CFITSLIB =
HPXINC =
HPXLIB =
PFLAGS = -DENABLE_LITE
LIBNAME = $(LIBNAME_LITE)
else
PFLAGS =
LIBNAME = $(LIBNAME_FULL)
endif

CFLAGS = $(DEBUG) -O3 -Wall -std=c99 -I. -I../sofa $(SLAINC) $(PFLAGS) -I/usr/local/include -I/opt/local/include $(HPXINC) $(CFITSINC) -fPIC -fopenmp
LDFLAGS = -L../sofa -lsofa_c $(SLALIB) -L/opt/local/lib $(HPXLIB) $(CFITSLIB) -lgomp

SRCS = $(filter-out %iers_bulletin_a.c, $(wildcard *.c)) qp_iers_bulletin_a.c

ifneq ($(PFLAGS), )
SRCS := $(filter-out qp_map.c qp_pixel.c, $(SRCS))
endif

OBJS = $(SRCS:.c=.o)
HEADERS = $(wildcard *.h)
HEADER = qpoint.h
ifneq ($(ENABLE_SHARED), )
LIB = $(LIBNAME).so
CFLAGS := -fPIC $(CFLAGS)
LDFLAGS := -shared $(LDFLAGS)
else
LIB = $(LIBNAME).a
endif
PREFIX?=/usr/local
INCPREFIX=$(PREFIX)/include
LIBPREFIX=$(PREFIX)/lib

default: all
all: $(LIB)

qp_iers_bulletin_a.c: iers_bulletin_a.dat make_iers_bulletin_a_dot_c.py
	python make_iers_bulletin_a_dot_c.py

$(LIB): $(OBJS)
ifneq ($(ENABLE_SHARED), )
	$(CC) $(LDFLAGS) $(OBJS) -o $@
else
	ar ru $@ $(OBJS)
endif

install: $(LIB) $(HEADER)
	install -m 644 $(HEADER) $(INCPREFIX)
	install -m 644 $(LIB) $(LIBPREFIX)

uninstall:
	cd $(INCPREFIX) && rm -f $(HEADER)
	cd $(LIBPREFIX) && rm -f $(LIBNAME).*

.PHONY: tidy clean clean-obj

tidy:
	rm -f *~

clean-obj:
	rm -f *.o

clean: tidy clean-obj
	rm -f *iers_bulletin_a.c
	rm -f $(LIBNAME_FULL).* $(LIBNAME_LITE).*
