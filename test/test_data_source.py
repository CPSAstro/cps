
import matplotlib.pyplot as plt
import numpy as np
from cps.data_source import DataSource,DataSpectral#get_cube,get_spectral

from astroquery.skyview import SkyView
from astropy.coordinates import SkyCoord

def test_init():
    ds = DataSource(SkyView, list_survey=True)

def test_query():
    position='Eta Carinae'
    survey=['Fermi 5', 'HRI', 'DSS']
    # ds = DataSource(SkyView)
    # position = SkyCoord(ra, dec, unit='deg',frame='fk5')
    # survey = 'WISE 3.4'
    ds = DataSource(SkyView, position=position, survey=survey)
    ds.query()

def test_images_and_headers():
    ds = DataSource(SkyView, position='Eta Carinae', survey='DSS')
    ds.query()
    print(ds.headers())

    data = ds.images()
    plt.imshow(data[0])
    plt.savefig('test/Figures/TestImage.png')

def test_masks_and_background():

    coord = SkyCoord(280.7156971, -4.0573715, unit='deg', frame='fk5')
    ds = DataSource(SkyView, position=coord, survey='WISE 3.4')
    ds.images()
    ds.set_mask(inner=15, outer=15*np.sqrt(2), pos=coord, method='annulus')
    print(ds.flux_info())


def test_DataSpectral():
    url=('http://jvo.nao.ac.jp/skynode/do/download/nobeyama/coming/coming_meta/CMG00000000')
    filename=('test_new.fits')
    dp= DataSpectral(url, filename)
    dp.query()
    line_spectrum=dp.study_spectral()
#    spectrum=get_spectral('test/input_test_file/test.fits')
    print(line_spectrum)
    
    
