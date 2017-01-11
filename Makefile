CC=gcc

# HACK
# need gfortran to build slarefro package
ifeq ($(shell which gfortran), )
SLAR =
else
ifneq ($(DISABLE_SLAREFRO), )
SLAR =
else
SLAR = slarefro
endif
endif

ifneq ($(SLAR), )
INSTALL_SLAR = install-$(SLAR)
INSTALL_SLAR_USER = install-$(SLAR)-user
UNINSTALL_SLAR = uninstall-$(SLAR)
endif

LOCALPREFIX=$(HOME)/.local

default: all

all: qpoint python

.PHONY: qpoint sofa slarefro chealpix python
qpoint: sofa $(SLAR) chealpix
	make -C src

qpoint-shared: sofa-shared $(SLAR) chealpix-shared
	ENABLE_SHARED=yes make -C src

qpoint-shared-lite: sofa-shared $(SLAR)
	ENABLE_SHARED=yes ENABLE_LITE=yes make -C src

sofa:
	make -C sofa

sofa-shared:
	ENABLE_SHARED=yes make -C sofa

slarefro:
	make -C slarefro

chealpix:
	make -C chealpix

chealpix-shared:
	ENABLE_SHARED=yes make -C chealpix-shared

install-sofa: sofa
	make -C sofa install

install-sofa-shared: sofa-shared
	ENABLE_SHARED=yes make -C sofa install

install-slarefro: slarefro
	make -C slarefro install

install-chealpix: chealpix
	make -C chealpix install

install-chealpix-shared: chealpix-shared
	ENABLE_SHARED=yes make -C chealpix install

python: sofa $(SLAR) chealpix
	CC=$(CC) python setup.py build

install-python: python
	python setup.py install

install-qpoint: qpoint install-sofa $(INSTALL_SLAR) install-chealpix
	make -C src install

install-qpoint-shared: qpoint-shared install-sofa-shared install-chealpix-shared
	ENABLE_SHARED=yes make -C src install

install-qpoint-shared-lite: qpoint-shared install-sofa-shared
	ENABLE_SHARED=yes ENABLE_LITE=yes make -C src install

install: install-qpoint install-python

install-all: install-qpoint install-python

install-python-user: python
	python setup.py install --user

install-sofa-user: sofa
	PREFIX=$(LOCALPREFIX) make -C sofa install

install-slarefro-user: slarefro
	PREFIX=$(LOCALPREFIX) make -C slarefro install

install-chealpix-user: chealpix
	PREFIX=$(LOCALPREFIX) make -C chealpix install

install-qpoint-user: install-sofa-user $(INSTALL_SLAR_USER) install-chealpix-user
	PREFIX=$(LOCALPREFIX) make -C src install

install-user: install-qpoint-user install-python-user

uninstall-sofa:
	make -C sofa uninstall

uninstall-slarefro:
	make -C slarefro uninstall

uninstall-chealpix:
	make -C chealpix uninstall

uninstall-qpoint:
	make -C src uninstall

uninstall: uninstall-qpoint

uninstall-all: uninstall-sofa $(UNINSTALL_SLAR) uninstall-slarefro uninstall-chealpix uninstall-qpoint

clean-sofa:
	make -C sofa clean

clean-slarefro:
	make -C slarefro clean

clean-chealpix:
	make -C chealpix clean

clean-python:
	python setup.py clean --all

clean-qpoint:
	make -C src clean

clean: clean-qpoint clean-python

clean-all: clean clean-python clean-sofa clean-slarefro clean-chealpix
