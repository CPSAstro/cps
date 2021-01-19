import astropy.coordinates as coord
import astropy.units as u
from astroquery.irsa import Irsa
from astroquery.vizier import Vizier

def cross_match():
    """
    Nothing but test query data and cross match.
    """
    twomass = Irsa.query_region(coord.SkyCoord(28.2,-0.049, unit=(u.deg,u.deg),frame='galactic'),catalog='fp_psc', radius='1d0m0s')

    v = Vizier(columns = ["**","RAJ2000","DEJ2000"])
    v.ROW_LIMIT = 9000000
    result = v.query_region(coord.SkyCoord(ra=280.7421273, dec=-4.2326516,\
    unit=(u.deg, u.deg), frame='icrs'),radius=0.05*u.deg, catalog=["GAIA DR2"])

    gaia_data = result[4]

    coo_wise = coord.SkyCoord(twomass['ra'], twomass['dec'])
    coo_twomass = coord.SkyCoord(gaia_data['RAJ2000'], gaia_data['DEJ2000'])
    idx_twomass, d2d_twomass, d3d_twomass = coo_wise.match_to_catalog_sky(coo_twomass)
    YSO = gaia_data[:0].copy()
    for k in idx_twomass:
        YSO.add_row(gaia_data[k])

    return YSO
