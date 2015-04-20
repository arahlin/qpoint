CC=gcc

# HACK
# need gfortran to build slarefro package
ifeq ($(shell which gfortran), )
SLAR =
INSTALL_SLAR =
UNINSTALL_SLAR =
else
SLAR = slarefro
INSTALL_SLAR = install-slarefro
UNINSTALL_SLAR = uninstall-slarefro
endif

default: all

all: qpoint python

.PHONY: qpoint sofa slarefro python
qpoint: sofa $(SLAR)
	make -C src

sofa:
	make -C sofa

slarefro:
	make -C slarefro

install-sofa: sofa
	make -C sofa install

install-slarefro: slarefro
	make -C slarefro install

python: sofa $(SLAR)
	CC=$(CC) python setup.py build

install-python: python
	python setup.py install

install-qpoint: qpoint install-sofa $(INSTALL_SLAR)
	make -C src install

install: install-qpoint install-python

install-all: install-qpoint install-python

uninstall-sofa:
	make -C sofa uninstall

uninstall-slarefro:
	make -C slarefro uninstall

uninstall-qpoint:
	make -C src uninstall

uninstall: uninstall-qpoint

uninstall-all: uninstall-sofa $(UNINSTALL_SLAR) uninstall-slarefro uninstall-qpoint

clean-sofa:
	make -C sofa clean

clean-slarefro:
	make -C slarefro clean

clean-python:
	python setup.py clean --all

clean-qpoint:
	make -C src clean

clean: clean-qpoint clean-python

clean-all: clean clean-python clean-sofa clean-slarefro
