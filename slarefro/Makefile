# 
# Set PREFIX before installing.
#

FSOURCE = $(wildcard *.f)
FOBJECT = $(FSOURCE:%.f=%.o)

CSOURCE = wrap.c
COBJECT = wrap.o

OBJECT = $(FOBJECT) $(COBJECT)

NAME = slarefro
# LIBNAME = lib$(NAME).so
LIBNAME = lib$(NAME).a
HDRNAME = $(NAME).h

FC = gfortran

default: $(LIBNAME) $(HDRNAME)

$(LIBNAME): $(OBJECT)
#	$(CC) -shared -fPIC $^ -o $@ -fPIC -shared
	ar -ru $@ $^

$(HDRNAME): wrap.h
	cp wrap.h $(HDRNAME)

%.o : %.c
	$(CC) -fPIC -g -c -o $@ $<

%.o : %.f
	$(FC) -fPIC -g -c -o $@ $<

install:
	if [ ! -d $(PREFIX)/lib ]; then install -d $(PREFIX)/lib; fi
	if [ ! -d $(PREFIX)/include ]; then install -d $(PREFIX)/include; fi
	install $(LIBNAME) $(PREFIX)/lib
	install $(HDRNAME) $(PREFIX)/include

uninstall:
	cd $(PREFIX)/lib && rm -f $(LIBNAME)
	cd $(PREFIX)/include && rm -f $(HDRNAME)

clean:
	-rm *.o $(LIBNAME) $(HDRNAME)
