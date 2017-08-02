CC=gcc

LOCALPREFIX = $(HOME)/.local

ifeq ($(PREFIX), )
PYTHONPREFIX = 
else
PYTHONPREFIX = --prefix=$(PREFIX)
endif

default: all

all: qpoint python

.PHONY: qpoint sofa chealpix python
qpoint: sofa chealpix
	make -C src

qpoint-lite: sofa
	ENABLE_LITE=yes make -C src

qpoint-shared: sofa-shared chealpix-shared
	ENABLE_SHARED=yes make -C src

qpoint-shared-lite: sofa-shared
	ENABLE_SHARED=yes ENABLE_LITE=yes make -C src

sofa:
	make -C sofa

sofa-shared:
	ENABLE_SHARED=yes make -C sofa

chealpix:
	make -C chealpix

chealpix-shared:
	ENABLE_SHARED=yes make -C chealpix

install-sofa: sofa
	make -C sofa install

install-sofa-shared: sofa-shared
	ENABLE_SHARED=yes make -C sofa install

install-chealpix: chealpix
	make -C chealpix install

install-chealpix-shared: chealpix-shared
	ENABLE_SHARED=yes make -C chealpix install

python: sofa chealpix
	CC=$(CC) python setup.py build

install-python: python
	python setup.py install $(PYTHONPREFIX)

install-qpoint: qpoint
	make -C src install

install-qpoint-lite: qpoint-lite

install-qpoint-shared: qpoint-shared
	ENABLE_SHARED=yes make -C src install

install-qpoint-shared-lite: qpoint-shared
	ENABLE_SHARED=yes ENABLE_LITE=yes make -C src install

install: install-qpoint install-python

install-all: install-sofa install-chealpix install-qpoint install-python

install-python-user: python
	python setup.py install --prefix=$(LOCALPREFIX)

install-sofa-user: sofa
	PREFIX=$(LOCALPREFIX) make -C sofa install

install-chealpix-user: chealpix
	PREFIX=$(LOCALPREFIX) make -C chealpix install

install-qpoint-user: qpoint
	PREFIX=$(LOCALPREFIX) make -C src install

install-user: install-qpoint-user install-python-user

uninstall-sofa:
	make -C sofa uninstall

uninstall-chealpix:
	make -C chealpix uninstall

uninstall-qpoint:
	make -C src uninstall

uninstall: uninstall-qpoint

uninstall-all: uninstall-sofa uninstall-chealpix uninstall-qpoint

clean-sofa:
	make -C sofa clean

clean-chealpix:
	make -C chealpix clean

clean-python:
	python setup.py clean --all

clean-qpoint:
	make -C src clean

clean: clean-qpoint clean-python

clean-all: clean clean-python clean-sofa clean-chealpix
