VPATH=../erfa:../chealpix:.
# DEBUG = -g
DEBUG =

LIBNAME_FULL = libqpoint
LIBNAME_LITE = libqpntlt

PYTHON ?= python
QP_VERSION = $(shell git describe --abbrev=4 --dirty --always --tags)

ERFAINC = -I../erfa
ERFALIB = # -L../erfa -lerfa_qp
ERFAOBJS = $(wildcard ../erfa/*.o)
ERFAHDRS = $(wildcard ../erfa/*.h)

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

ifneq ($(ENABLE_OMP), )
ifeq ($(CC),icc)
OPENMP = -qopenmp
OPENMP_LD =
else
OPENMP = -fopenmp
OPENMP_LD = -lgomp
endif
endif

CFLAGS = $(DEBUG) -O3 -Wall -std=c99 -DQP_VERSION=\"$(QP_VERSION)\"
# add some options even when CFLAGS overridden on command line
override CFLAGS += -std=c99 -I. $(ERFAINC) $(HPXINC) $(PFLAGS) -fPIC $(OPENMP)
LDFLAGS = $(ERFALIB) $(HPXLIB) $(OPENMP_LD)

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
	$(PYTHON) make_iers_bulletin_a_dot_c.py

$(OBJS): %$(OBJNAME): %.c
	$(CC) -c $(CFLAGS) $< -o $@

$(LIB): $(ERFAOBJS) $(HPXOBJS) $(OBJS)
ifneq ($(ENABLE_SHARED), )
	$(CC) $(LDFLAGS) $(ERFAOBJS) $(HPXOBJS) $(OBJS) -o $@
else
	ar ru $@ $(ERFAOBJS) $(HPXOBJS) $(OBJS)
endif

install: $(LIB) $(ERFAHDRS) $(HPXHDRS) $(HEADER)
	if [ ! -d $(INCPREFIX) ]; then install -d $(INCPREFIX); fi
	install $(ERFAHDRS) $(HPXHDRS) $(HEADER) $(INCPREFIX)
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
