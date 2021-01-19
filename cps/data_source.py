# from cps import annulus_class
import requests
import logging
import urllib
import numpy as np
import astropy.units as u
from io import BytesIO
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.stats import SigmaClip, sigma_clipped_stats
from photutils import SkyCircularAperture, SkyCircularAnnulus, ApertureMask, \
                      aperture_photometry, make_source_mask

class DataSource:
    def __init__(self, database: object, *args, list_survey=False, **kwargs) -> None:
        """
        Astronomical data source

        Args:
            database (object): one of database provide by astroquery
        """
        self.database = database
        self.args = args
        self.kwargs = kwargs
        self.imgs = None

        # TODO: reduce the usage of aper or mask. The function is repeated right now
        self.aper = None
        self.mask = None
        self.mask_member = None

        if list_survey:
            database.list_surveys()

    # TODO: Use generator?
    def query(self, *args, **kwargs):
        if args or kwargs:
            self.args = args
            self.kwargs = kwargs
        self.imgs = list(self.database.get_images(*self.args, **self.kwargs))

    def headers(self):
        """
        Get general information of images

        Returns:
            header (str): the header of images
        """
        return [WCS(x[0].header) for x in self.imgs]

    def images(self):
        if not self.imgs:
            self.query()
        return [x[0].data for x in self.imgs]

    def set_mask(self, inner: float = None, outer: float = None, pos: SkyCoord = None, wcs: WCS = None, source: str = 'all', method: str = 'source_mask', *args, **kwargs):

        # TODO: rename the variable
        self.mask_member = source

        if source == 'each':
            data_list = self.images()

            if method == 'source_mask':
                self.mask = [
                    make_source_mask(data, args, kwargs)
                    for data in data_list
                ]

            if method == 'annulus':

                raise NotImplementedError

        elif source == 'all':
            if method == 'source_mask':
                logging.info("Use the first image source to set source mask. Could not compatible with other sources")

                data = self.images()[0]
                self.mask = make_source_mask(data, args, kwargs)

            elif method == 'annulus':

                if not (inner and outer and pos):
                    # TODO: support saving coordinate when initializing
                    raise TypeError

                if not wcs:
                    logging.info("Not assign the wcs, use the wcs from the header of the first image source")
                    wcs = self.headers()[0]

                # TODO: should allow input different parameters
                self.aper = SkyCircularAnnulus(pos, inner*u.arcsec, outer*u.arcsec).to_pixel(wcs)
                self.mask = self.aper.to_mask('center')

        else:
            # TODO: Update error type
            raise NameError

    @staticmethod
    def _calc_flux(data, *args, method: str = 'rms', **kwargs):

        if method == 'rms':
            return np.sqrt(np.average(data**2))
        
        elif method == 'median':
            _, median, _ = sigma_clipped_stats(data, args, kwargs)
            return median

        else:
            raise NameError

    @classmethod
    def _center_flux(cls, image: object, mask: object, *args, **kwargs):

        if isinstance(mask, SkyCircularAnnulus):
            pos = mask.positions
            inner = mask.r_in
            center_mask = SkyCircularAperture(inner*u.arcsec).to_pixel(pos).to_mask('center')

            data = image[center_mask]

            return cls._calc_flux(data, args, kwargs)

        elif isinstance(mask, SkyCircularAperture):

            data = image[mask]

            return cls._calc_flux(data, args, kwargs)

    @classmethod
    def _background_flux(cls, image: object, mask: object, *args, **kwargs):
        
        annulus_data = mask.multiply(image)
        data = annulus_data[mask.data>0]

        return cls._calc_flux(data, args, kwargs)

    def backgrounds(self, *args, **kwargs):

        if not self.mask_member:
            raise RuntimeError

        elif self.mask_member == 'all':
            return [self._background_flux(image, self.mask, args, kwargs) for image in self.images()]

        elif self.mask_member == 'each':
            return [self._background_flux(image, mask, args, kwargs) for image, mask in zip(self.images, self.mask)]

        else:
            raise RuntimeError

    def flux_info(self, *args, **kwargs):

        if self.aper == None:

            raise TypeError

        else:

            if self.mask_member == 'each':
                
                raise NotImplementedError

            elif self.mask_member == 'all':
                infos = [aperture_photometry(image, self.aper) for image in self.images()]

            if not infos:
                raise RuntimeError

            bkgs = self.backgrounds(args, kwargs)
            for idx, info in enumerate(infos):
                info['background_flux'] = bkgs[idx]

        return infos

def get_cube(url: str):
    urllib.request.urlopen(url)
    r = requests.get(url)
    open('test.fits', 'wb').write(r.content)
                
    
def get_spectral(filename: str,vindex:int='3'):
    file_fits = fits.open(filename)
    spectral_data = file_fits[0].data
    line_header=file_fits[0].header
    vaxis=f'CRVAL{vindex}'
    pix_axis=f'CDELT{vindex}'
    averaged_spectrum= np.nanmean(np.nanmean(spectral_data,2),1)
    velocity_axis_end=line_header[vaxis]+len(averaged_spectrum)*line_header[pix_axis]
    v_axes =np.linspace(line_header[vaxis],velocity_axis_end,len(averaged_spectrum))

    return averaged_spectrum
        
    
        
              
              
               
                    
