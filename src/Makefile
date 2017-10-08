VPATH=../sofa:../chealpix:.
# DEBUG = -g
DEBUG =
CC = gcc

LIBNAME_FULL = libqpoint
LIBNAME_LITE = libqpntlt

QP_VERSION = $(shell git describe --abbrev=4 --dirty --always --tags)

SOFAINC = -I../sofa
SOFALIB = # -L../sofa -lsofa_c
SOFAOBJS = $(wildcard ../sofa/*.o)
SOFAHDRS = $(wildcard ../sofa/*.h)

ifneq ($(ENABLE_LITE), )
HPXINC =
HPXLIB =
HPXOBJS =
HPXHDRS =
PFLAGS = -DENABLE_LITE
LIBNAME = $(LIBNAME_LITE)
OBJNAME = -lt.o
else
HPXINC = -I../chealpix
HPXLIB = # -L../chealpix -lchealpix_qp
HPXOBJS = $(wildcard ../chealpix/*.o)
HPXHDRS = $(wildcard ../chealpix/*.h)
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
override CFLAGS += -std=c99 -I. $(SOFAINC) $(HPXINC) $(PFLAGS) -fPIC $(OPENMP)
LDFLAGS = $(SOFALIB) $(HPXLIB) $(OPENMP_LD)

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
INCPREFIX=$(PREFIX)/include/qpoint
LIBPREFIX=$(PREFIX)/lib

default: all
all: $(LIB)

qp_iers_bulletin_a.c: make_iers_bulletin_a_dot_c.py
	python make_iers_bulletin_a_dot_c.py

$(OBJS): %$(OBJNAME): %.c
	$(CC) -c $(CFLAGS) $< -o $@

$(LIB): $(SOFAOBJS) $(HPXOBJS) $(OBJS)
ifneq ($(ENABLE_SHARED), )
	$(CC) $(LDFLAGS) $(SOFAOBJS) $(HPXOBJS) $(OBJS) -o $@
else
	ar ru $@ $(SOFAOBJS) $(HPXOBJS) $(OBJS)
endif

install: $(LIB) $(SOFAHDRS) $(HPXHDRS) $(HEADER)
	if [ ! -d $(INCPREFIX) ]; then install -d $(INCPREFIX); fi
	install $(SOFAHDRS) $(HPXHDRS) $(HEADER) $(INCPREFIX)
	if [ ! -d $(LIBPREFIX) ]; then install -d $(LIBPREFIX); fi
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
