CC=gcc

default: all

all: qpoint python

.PHONY: qpoint sofa slarefro python
qpoint: sofa slarefro
	make -C src

sofa:
	make -C sofa

slarefro:
	make -C slarefro

install-sofa: sofa
	make -C sofa install

install-slarefro: slarefro
	make -C slarefro install

python: sofa slarefro
	CC=$(CC) python setup.py build

install-python: python
	python setup.py install

install: qpoint
	make -C src install

install-all: install-sofa install install-python

uninstall-sofa:
	make -C sofa uninstall

uninstall-slarefro:
	make -C slarefro uninstall

uninstall:
	make -C src uninstall

uninstall-all: uninstall-sofa uninstall-slarefro uninstall

clean-sofa:
	make -C sofa clean

clean-slarefro:
	make -C slarefro clean

clean-python:
	python setup.py clean --all

clean:
	make -C src clean

clean-all: clean clean-python clean-sofa clean-slarefro
