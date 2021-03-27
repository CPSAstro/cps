Requirements
============

The cpsAStro package requires the following packages:

Python 3.8 or later
`Numpy 1.16  <https://numpy.org/devdocs/release/1.16.0-notes.html>`_ or later

`Astropy 2.0  <https://www.astropy.org/announcements/release-2.0.html>`_  or later

`astrquery 0.4.2 <https://astroquery.readthedocs.io/en/latest/>`_  or later

`photutils 1.1.0 <https://photutils.readthedocs.io/en/stable/>`_  or later

`datackasses 3.8 <https://docs.python.org/3/library/dataclasses.html>`_  or later

`requests 2.25.1 <https://pypi.org/project/requests/>`_  or later

`logging 3.9.2 <https://docs.python.org/3/howto/logging.html>`_  or later

The following packahe is optional: 

`spectral-cube 0.5 <https://spectral-cube.readthedocs.io/en/latest/#>`_  or later


Installation
============
To install the tool
::

	git clone https://github.com/CPSAstro/cps
	cd cps
	pip install .


For development, run the following comand to install the dependencies.
::

	poetry install




For testing, run 
::

	poetry run pytest -s test/






