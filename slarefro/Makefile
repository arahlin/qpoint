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
	install -d $(PREFIX)/lib $(PREFIX)/include
	install $(LIBNAME) $(PREFIX)/lib
	install $(HDRNAME) $(PREFIX)/include


clean:
	-rm *.o $(LIBNAME) $(HDRNAME)

