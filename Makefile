LOCALPREFIX = $(HOME)/.local

PYTHON ?= python

ifeq ($(PREFIX), )
PYTHONPREFIX = 
else
PYTHONPREFIX = --prefix=$(PREFIX)
endif

default: all

all: qpoint

.PHONY: qpoint erfa chealpix python docs
qpoint: erfa chealpix
	make -C src

qpoint-debug: erfa chealpix
	CFLAGS="-g -DDEBUG" make -C src

qpoint-lite: erfa
	ENABLE_LITE=yes make -C src

qpoint-shared: erfa-shared chealpix-shared
	ENABLE_SHARED=yes make -C src

qpoint-shared-lite: erfa-shared
	ENABLE_SHARED=yes ENABLE_LITE=yes make -C src

erfa:
	make -C erfa

erfa-shared:
	ENABLE_SHARED=yes make -C erfa

chealpix:
	make -C chealpix

chealpix-shared:
	ENABLE_SHARED=yes make -C chealpix

install-erfa: erfa
	make -C erfa install

install-erfa-shared: erfa-shared
	ENABLE_SHARED=yes make -C erfa install

install-chealpix: chealpix
	make -C chealpix install

install-chealpix-shared: chealpix-shared
	ENABLE_SHARED=yes make -C chealpix install

python:
	$(PYTHON) setup.py build

python-debug:
	CFLAGS="-g -DDEBUG" $(PYTHON) setup.py build

docs:
	$(PYTHON) setup.py build_sphinx -a

install-python: python
	$(PYTHON) setup.py install $(PYTHONPREFIX)

install-qpoint: qpoint
	make -C src install

install-qpoint-lite: qpoint-lite

install-qpoint-shared: qpoint-shared
	ENABLE_SHARED=yes make -C src install

install-qpoint-shared-lite: qpoint-shared
	ENABLE_SHARED=yes ENABLE_LITE=yes make -C src install

install: install-qpoint install-python

install-all: install-erfa install-chealpix install-qpoint install-python

install-python-user: python
	$(PYTHON) setup.py install --prefix=$(LOCALPREFIX)

install-erfa-user: erfa
	PREFIX=$(LOCALPREFIX) make -C erfa install

install-chealpix-user: chealpix
	PREFIX=$(LOCALPREFIX) make -C chealpix install

install-qpoint-user: qpoint
	PREFIX=$(LOCALPREFIX) make -C src install

install-user: install-qpoint-user install-python-user

uninstall-erfa:
	make -C erfa uninstall

uninstall-chealpix:
	make -C chealpix uninstall

uninstall-qpoint:
	make -C src uninstall

uninstall: uninstall-qpoint

uninstall-all: uninstall-erfa uninstall-chealpix uninstall-qpoint

clean-erfa:
	make -C erfa clean

clean-chealpix:
	make -C chealpix clean

clean-python:
	$(PYTHON) setup.py clean --all

clean-qpoint:
	make -C src clean

clean: clean-qpoint clean-python

clean-all: clean clean-python clean-erfa clean-chealpix
