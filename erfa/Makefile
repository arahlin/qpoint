DEBUG= -g
TESTSRC=test_erfa.c
TESTOBJ=$(TESTSRC:.c=.o)
TESTTARG=$(TESTSRC:.c=)
SRCS=$(filter-out $(TESTSRC),$(wildcard *.c))
OBJS=$(SRCS:.c=.o)
HEADERS=$(wildcard *.h)
LIBNAME=erfa_qp
PREFIX?=/usr/local

CC = gcc
CFLAGS = -pedantic -Wall -O3 $(DEBUG)
# add some options even when CFLAGS overridden on command line
override CFLAGS += -fPIC -I.

ifeq ($(ENABLE_SHARED), )
TARGET=lib$(LIBNAME).a
else
TARGET=lib$(LIBNAME).so
LDFLAGS := -shared $(LDFLAGS)
endif

all: lib

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

lib: $(TARGET)

$(TARGET): $(OBJS) $(HEADERS)
ifeq ($(ENABLE_SHARED), )
	$(AR) ru $(TARGET) $(OBJS)
else
	$(CC) $(LDFLAGS) $(OBJS) -o $@
endif

test: $(TESTTARG)

$(TESTTARG): $(TESTOBJ) $(TARGET)
	$(CC) -o $@ $(TESTOBJ) -L. -l$(LIBNAME)

.PHONY: clean tidy
clean:
	rm -f *.o $(TARGET) $(TESTTARG)

install: $(TARGET) $(HEADERS)
	if [ ! -d $(PREFIX)/lib ]; then install -d $(PREFIX)/lib; fi
	if [ ! -d $(PREFIX)/include/qpoint ]; then install -d $(PREFIX)/include/qpoint; fi
	install $(TARGET) $(PREFIX)/lib/
	install $(HEADERS) $(PREFIX)/include/qpoint/

.PHONY: uninstall
uninstall:
	cd $(PREFIX)/lib && rm -f $(TARGET)
	cd $(PREFIX)/include && rm -rf $(HEADERS)
