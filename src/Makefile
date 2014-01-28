VPATH=../sofa:../slarefro:.
DEBUG = -g -O3
# DEBUG = -O3
CC = gcc
CFLAGS = $(DEBUG) -Wall -std=c99 -I. -I../sofa -I../slarefro -I/usr/local/include -fPIC
LDFLAGS = -L../sofa -lsofa_c -L../slarefro -lslarefro
SRCS = $(filter-out test_qpoint.c test_math.c,$(wildcard *.c))
OBJS = $(SRCS:.c=.o)
HEADERS = $(wildcard *.h)
HEADER = qpoint.h
LIB = libqpoint.a
PREFIX?=/usr/local

default: all
all: $(LIB)

test: test_qpoint test_math

test_qpoint: test_qpoint.o $(LIB)
	gcc $(DEBUG) -o $@ $< $(LDFLAGS) -lgetdata -L. -lqpoint

test_math: test_math.o sincos.o
	gcc $(DEBUG) -o $@ $< sincos.o

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
