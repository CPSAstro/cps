Usage
==================
For the continuum intensity measurement, we have implmented the Skyview, Herschel, and Spitzer astroquery services as these are the most widely used images for SED analyses. 
::

	import matplotlib.pyplot as plt
	import numpy as np
	from astropy.coordinates import SkyCoord
	from cpsastro.data_source import DataImage, DataCube
	
	ds = DataImage('SkyView', 'R CrA Cloud', 'AKARI N60')
	data = ds.data
	plt.imshow(data)
	plt.savefig("test/Figures/TestImage.png")

When we measure the intensity, it is normal to have comtaimination from the background or neighboring sources, hence we apply photutils to help user definineing an annulus region from the aperture with manuual selected radius and the central coordiante of the aperture. Here we assume the aperture radius is 15 arcsec and the annulus has a radius such that the area of the annulus is the same as the aperture. 
::

	coord = SkyCoord(359.901, -17.853, unit='deg', frame='galactic') 
	ds.set_mask(inner=15, outer=15*1.4, position=coord, method='annulus')


For the spectral cube anlyses, many of the these cubes are not in public archive, and normally require manual download. Hence, here we require the user to use the url link to the file and download to disk before performing the anlyses. 
::

	url = ("http://jvo.nao.ac.jp/skynode/do/download/nobeyama/coming/coming_meta/CMG00000000")
	dp = data_source.DataCube(url)
	dp.query()
	print(dp.spectrum)
