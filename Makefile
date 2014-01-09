default: all

all: qpoint python

.PHONY: qpoint sofa python
qpoint: sofa
	make -C src

sofa:
	make -C sofa

install-sofa: sofa
	make -C sofa install

python: sofa
	python setup.py build

install-python: python
	python setup.py install

install: qpoint
	make -C src install

install-all: install-sofa install install-python

uninstall-sofa:
	make -C sofa uninstall

uninstall:
	make -C src uninstall

uninstall-all: uninstall-sofa uninstall

clean-sofa:
	make -C sofa clean

clean-python:
	python setup.py clean --all

clean:
	make -C src clean

clean-all: clean clean-python clean-sofa
