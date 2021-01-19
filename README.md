# CPSASTRO (Combined Photometric-Spectra python analysis tool)

## Synoptic introduction
The aim of this package is core of a more full version of a to be astronomical data analysis package with a clear and as transparaent as possible techqniue in obtaining various crucial physical parameters that store in various astronomical distributions (e.g. SED; Turbulence spectral power spectrum, probability density function, etc). The goal is to allow the data being analysed through these routines in a relatively uniform way, but still retain the flexibility for custom analyses on top the fundamental routines. This current version provides a relatively effective way to query intensity maps from archival astronomical surveys (e.g. Herschel, WISE) and measure intensity as a function of wavelengths, and form an important product: spectral energy distribution (SED), which allow use uncover loads of physical properties of the astronomical object (e.g. evolutionary stage, mass, etc). The package can also obtain spectral information from position-position-velocity cubes, with more deticated gas velocity analyses routine to be added in later.  

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

