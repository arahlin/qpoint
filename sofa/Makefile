DEBUG= -g
TESTSRC=t_sofa_c.c
TESTOBJ=$(TESTSRC:.c=.o)
TESTTARG=$(TESTSRC:.c=)
SRCS=$(filter-out $(TESTSRC),$(wildcard *.c))
OBJS=$(SRCS:.c=.o)
HEADERS=$(wildcard *.h)
LIBNAME=sofa_c
PREFIX?=/usr/local

CC = gcc
CFLAGS = -pedantic -Wall -O3 -fPIC -I. $(DEBUG)

ifeq ($(ENABLE_SHARED), )
TARGET=lib$(LIBNAME).a
else
TARGET=lib$(LIBNAME).so
LDFLAGS := -shared $(LDFLAGS)
CFLAGS := -fPIC $(CFLAGS)
endif

all: lib

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

lib: $(TARGET)

$(TARGET): $(OBJS) $(HEADERS)
ifeq ($(ENABLE_SHARED), )
	ar ru $(TARGET) $(OBJS)
else
	gcc $(LDFLAGS) $(OBJS) -o $@
endif

test: $(TESTTARG)

$(TESTTARG): $(TESTOBJ) $(TARGET)
	gcc -o $@ $(TESTOBJ) -L. -l$(LIBNAME)

.PHONY: clean tidy
clean:
	rm -f *.o $(TARGET) $(TESTTARG)

install: $(TARGET) $(HEADERS)
	if [ ! -d $(PREFIX)/lib ]; then install -d $(PREFIX)/lib; fi
	if [ ! -d $(PREFIX)/include ]; then install -d $(PREFIX)/include; fi
	install $(TARGET) $(PREFIX)/lib/
	install $(HEADERS) $(PREFIX)/include/

.PHONY: uninstall
uninstall:
	cd $(PREFIX)/lib && rm -f $(TARGET)
	cd $(PREFIX)/include && rm -rf $(HEADERS)
