DEBUG= -g
TESTSRC=t_sofa_c.c
TESTOBJ=$(TESTSRC:.c=.o)
TESTTARG=$(TESTSRC:.c=)
SRCS=$(filter-out $(TESTSRC),$(wildcard *.c))
OBJS=$(SRCS:.c=.o)
HEADERS=$(wildcard *.h)
LIBNAME=sofa_c
TARGET=lib$(LIBNAME).a
PREFIX?=/usr/local

CC = gcc
CFLAGS = -pedantic -Wall -O3 -fPIC -I. $(DEBUG)

all: lib

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

lib: $(TARGET)

$(TARGET): $(OBJS) $(HEADERS)
	ar ru $(TARGET) $(OBJS)

test: $(TESTTARG)

$(TESTTARG): $(TESTOBJ) $(TARGET)
	gcc -o $@ $(TESTOBJ) -L. -l$(LIBNAME)

.PHONY: clean tidy
clean:
	rm -f *.o $(TARGET) $(TESTTARG)

install: $(TARGET) $(HEADERS)
	install -m 644 $(TARGET) $(PREFIX)/lib/$(TARGET)
	install -m 644 $(HEADERS) $(PREFIX)/include/

.PHONY: uninstall
uninstall:
	cd $(PREFIX)/lib && rm -f $(TARGET)
	cd $(PREFIX)/include && rm -rf $(HEADERS)
