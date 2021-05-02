import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord

from cpsastro.data_source import DataImage, DataCube


def test_dataimage_query():
    position = "Eta Carinae"
    survey = ["Fermi 5", "HRI", "DSS"]
    # ds = DataSource(SkyView)
    # position = SkyCoord(ra, dec, unit='deg',frame='fk5')
    # survey = 'WISE 3.4'
    dslist = [DataImage("skyview", position=position, survey=s) for s in survey]
    for ds in dslist:
        ds.query()
    assert ds.data.any()

def test_dataimage_data_and_header():
    ds = DataImage("skyview", position="Eta Carinae", survey="DSS")
    print(ds.header)
    assert ds.header.naxis==2
    data = ds.data
    plt.imshow(data)
    plt.savefig("test/Figures/TestImage.png")


def test_dataimage_masks_and_flux():

    coord = SkyCoord(280.7156971, -4.0573715, unit="deg", frame="fk5")
    ds = DataImage("skyview", position=coord, survey="WISE 3.4")
    ds.set_mask(inner=15, outer=15 * np.sqrt(2), position=coord, method="annulus")
    print(type(ds._mask))
    print(ds.center_flux())
    print(ds.background_flux())
    assert ds.center_flux().any()

def test_datacube_spectrum():
    url = "http://jvo.nao.ac.jp/skynode/do/download/nobeyama/coming/coming_meta/CMG00000000"
    dp = DataCube(url)
    dp.query()
    print(dp.spectrum)
    assert dp.data.any()

    
