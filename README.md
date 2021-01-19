# CPS (Combined Photometric-Spectra python analysis tool)

## About
This package aims to allow a effective way to query intensity maps from archival astronomical surveys (e.g. Herschel, WISE) and measure intensity as a function of wavelengths. An important product is the spectral energy distribution (SED), which allow use uncover loads of physical properties of the astronomical object (e.g. evolutionary stage, mass, etc). The package can also obtain spectral information from position-position-velocity cubes (still in progress...). 

## Installation
To install the tool
```
$ git clone https://github.com/CPSAstro/cps
$ cd cps
$ pip install .
```

For development, run
```
$ poetry install
```
to install the dependencies

For testing, run 
```
$ poetry run pytest -s test/
```

## Usage
Currently, you have to assign the interested `DataSet` and `Survey` supported by <span style="color:red">*`astroquery`*</span> when you using the package.

