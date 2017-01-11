VPATH=../sofa:../chealpix:.
DEBUG = -g
# DEBUG =
CC = gcc

# HACK
# need gfortran for slarefro package
ifneq ($(shell which gfortran), )
ifeq ($(DISABLE_SLAREFRO), )
ifeq ($(ENABLE_LITE), )
VPATH = ../sofa:../chealpix:../slarefro:.
SLAINC = -I../slarefro -DENABLE_SLAREFRO
SLALIB = -L../slarefro -lslarefro
endif
endif
endif

LIBNAME_FULL = libqpoint
LIBNAME_LITE = libqpntlt

QP_VERSION = $(shell git describe --abbrev=4 --dirty --always --tags)

ifneq ($(ENABLE_LITE), )
HPXINC =
HPXLIB =
PFLAGS = -DENABLE_LITE
LIBNAME = $(LIBNAME_LITE)
OBJNAME = -lt.o
else
HPXINC = -I../chealpix
HPXLIB = -L../chealpix -lchealpix_qp
PFLAGS =
LIBNAME = $(LIBNAME_FULL)
OBJNAME = .o
endif

ifeq ($(CC),icc)
OPENMP = -qopenmp
OPENMP_LD =
else
OPENMP = -fopenmp
OPENMP_LD = -lgomp
endif

CFLAGS = $(DEBUG) -O3 -Wall -std=c99 -DQP_VERSION=\"$(QP_VERSION)\"
# add some options even when CFLAGS overridden on command line
override CFLAGS += -std=c99 -I. -I../sofa $(SLAINC) $(HPXINC) $(PFLAGS) -fPIC $(OPENMP)
LDFLAGS = -L../sofa -lsofa_c $(SLALIB) $(HPXLIB) $(OPENMP_LD)

SRCS = $(filter-out %iers_bulletin_a.c, $(wildcard *.c)) qp_iers_bulletin_a.c

ifneq ($(PFLAGS), )
SRCS := $(filter-out qp_map.c qp_pixel.c qp_pixhash.c, $(SRCS))
endif

OBJS = $(SRCS:.c=$(OBJNAME))
HEADERS = $(wildcard *.h)
HEADER = qpoint.h
ifneq ($(ENABLE_SHARED), )
LIB = $(LIBNAME).so
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

$(OBJS): %$(OBJNAME): %.c
	$(CC) -c $(CFLAGS) $< -o $@

$(LIB): $(OBJS)
ifneq ($(ENABLE_SHARED), )
	$(CC) $(LDFLAGS) $(OBJS) -o $@
else
	ar ru $@ $(OBJS)
endif

install: $(LIB) $(HEADER)
	install $(HEADER) $(INCPREFIX)
	install $(LIB) $(LIBPREFIX)

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
