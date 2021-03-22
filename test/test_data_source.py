import matplotlib.pyplot as plt
import numpy as np
from astroquery.skyview import SkyView
from astropy.coordinates import SkyCoord

from cpsastro.data_source import DataSource, DataCube


def test_init():
    ds = DataSource(SkyView, list_survey=True)


def test_query():
    position = "Eta Carinae"
    survey = ["Fermi 5", "HRI", "DSS"]
    # ds = DataSource(SkyView)
    # position = SkyCoord(ra, dec, unit='deg',frame='fk5')
    # survey = 'WISE 3.4'
    ds = DataSource(SkyView, position=position, survey=survey)
    ds.query()


def test_images_and_headers():
    ds = DataSource(SkyView, position="Eta Carinae", survey="DSS")
    ds.query()
    print(ds.headers())

    data = ds.images()
    plt.imshow(data[0])
    plt.savefig("test/Figures/TestImage.png")


def test_masks_and_background():

    coord = SkyCoord(280.7156971, -4.0573715, unit="deg", frame="fk5")
    ds = DataSource(SkyView, position=coord, survey="WISE 3.4")
    ds.images()
    ds.set_mask(inner=15, outer=15 * np.sqrt(2), pos=coord, method="annulus")
    print(ds.flux_info())


def test_datacube_spectrum():
    url = "http://jvo.nao.ac.jp/skynode/do/download/nobeyama/coming/coming_meta/CMG00000000"
    dp = DataCube(url)
    dp.query()
    print(dp.spectrum)
