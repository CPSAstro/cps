import requests
import logging
import numpy as np
import astropy.units as u
from os import path
from dataclasses import dataclass
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.stats import SigmaClip, sigma_clipped_stats
from photutils import (
    SkyCircularAperture,
    SkyCircularAnnulus,
    aperture_photometry,
    make_source_mask,
)

from astroquery.irsa import Irsa
from astroquery.esasky import ESASky
from astroquery.skyview import SkyView


class DataImage:
    """
    This class handles the file in some archived databases and allow
    the analysis on those object. One file should be handle by one
    instance.

    Attributes:
        database: the target database
        position: the sky coordinate of the object or its alias
        survey: the target survey

    """

    @dataclass
    class ImageMask:
        """
        The dataclass stores the information of mask.
        """

        shape: str = None
        position: object = None
        rin: float = None
        rout: float = None
        mask: object = None

    supported_archives = {
        "irsa": Irsa,
        "esasky": ESASky,
        "skyview": SkyView,
    }

    def __init__(
        self, database: str, position: SkyCoord = None, survey: str = None
    ) -> None:
        """
        Astronomical image map (2D)

        :param database: the source archives (allowed: skyview, irsa, esasky)
        :type database: str
        :param position: the sky coordinate of the object (or alias), defaults to None
        :type position: SkyCoord, optional
        :param survey: the target survey, defaults to None
        :type survey: str, optional
        """
        self.database = self.supported_archives.get(database.lower())
        self.position = position
        self.survey = survey
        self._image = None
        self._mask = None

    @property
    def image(self):
        """
        The image requested from the database. Ideally should be a HDUList

        :return: image map content
        :rtype: HDUList
        """
        if self._image:
            return self._image
        else:
            self.query()
            return self._image

    @property
    def surveylist(self):
        """
        The list of survey in the archived database

        :return: the list of survey
        :rtype: string
        """
        return self.database.list_surveys()

    def query(self, position: SkyCoord = None, survey: str = None):
        """
        request data from the database

        :param position: the sky coordinate of the object or alias, the value passed into the constructor will be try if no value is assigned here, defaults to None
        :type position: SkyCoord, optional
        :param survey: the target survey, the value passed into the constructor will be try if no value is assigned here, defaults to None
        :type survey: str, optional
        :raises ValueError: [description]
        """

        if position and survey:
            logging.info("Update the position and survey")
            self.position = position
            self.survey = survey

        if not self.position or not self.survey:
            raise ValueError("Missing position or survey")

        # Only get one image
        self._image = self.database.get_images(
            position=self.position, survey=self.survey
        )[0]

    @property
    def header(self):
        """
        Get the header (WCS) of the image

        :return: the header of the image
        :rtype: WCS
        """
        return WCS(self.image[0].header)

    @property
    def data(self):
        """
        Get the raw data of the image

        :return: raw data
        :rtype: HDU
        """
        return self.image[0].data

    def set_mask(
        self,
        *args,
        method: str = "source",
        inner: float = None,
        outer: float = None,
        position: SkyCoord = None,
        wcs: WCS = None,
        **kwargs,
    ):
        """
        Set a mask to the data image

        :param method: the method to set mask, "source", or "annulus"defaults to "source"
        :type method: str, optional
        :param inner: the inner radius of the annulus, defaults to None
        :type inner: float, optional
        :param outer: the outer radius of the annulus, defaults to None
        :type outer: float, optional
        :param position: the center of the annulus, defaults to None
        :type position: SkyCoord, optional
        :param wcs: the wcs used to convert mask to pixel, defaults to None
        :type wcs: WCS, optional
        :raises ValueError: "annulus" must have position and inner/outer radius
        """

        if method == "source":
            self._mask = make_source_mask(self.data, *args, **kwargs)

        elif method == "annulus":
            if not position or not inner or not outer:
                raise ValueError(
                    "Use annulus mask but missing position or inner/outer radius"
                )

            else:
                if not wcs:
                    logging.warning("Not assign the wcs, try the wcs from the header")
                    wcs = self.header

                mask = (
                    SkyCircularAnnulus(position, inner * u.arcsec, outer * u.arcsec)
                    .to_pixel(wcs)
                    .to_mask("center")
                )
                self._mask = ImageMask(
                    shape="annulus", position=position, rin=inner, rout=outer, mask=mask
                )
                # self.mask = self.aper.to_mask("center")

        else:
            raise ValueError(f"No supported method: {method}")

    @staticmethod
    def _calc_flux(data, *args, method: str = "rms", **kwargs):
        """
        calculate the flux of the input data

        :param data: the input data
        :type data: array
        :param method: the method to calculate flux, "rms" or "median", defaults to "rms"
        :type method: str, optional
        :return: the flux of the input data
        :rtype: [type]
        """

        if method == "rms":
            return np.sqrt(np.average(data ** 2))

        elif method == "median":
            _, median, _ = sigma_clipped_stats(data, *args, **kwargs)
            return median

        else:
            return NotImplemented

    def center_flux(self, *args, **kwargs):
        """
        calculate the flux in the center if the mask is "annulus" or "circle"

        :return: the flux in the center
        :rtype: float
        """

        if self._mask.shape == "annulus":

            pos = self._mask.position
            wcs = self.header
            inner = self._mask.rin
            center_mask = (
                SkyCircularAperture(pos, inner * u.arcsec)
                .to_pixel(wcs)
                .to_mask("center")
            )

            data = center_mask.multiply(self.data)
            data = data[center_mask.data > 0]

            return self._calc_flux(data, *args, **kwargs)

        elif self._mask.shape == "circle":

            data = self._mask.mask.multiply(self.data)
            data = data[self._mask.mask.data > 0]

            return self._calc_flux(data, *args, **kwargs)

    def background_flux(self, *args, **kwargs):
        """
        return the background flux if the mask is "annulus"

        :return: the background flux
        :rtype: float
        """

        if self._mask.shape != "annulus":
            raise ValueError("The mask is not an annulus")

        data = self._mask.mask.multiply(self.data)
        data = data[self._mask.mask.data > 0]

        return self._calc_flux(data, *args, **kwargs)

    # def backgrounds(self, *args, **kwargs):

    #     if not self.mask_member:
    #         raise RuntimeError

    #     elif self.mask_member == "all":
    #         return [
    #             self._background_flux(image, self.mask, args, kwargs)
    #             for image in self.images()
    #         ]

    #     elif self.mask_member == "each":
    #         return [
    #             self._background_flux(image, mask, args, kwargs)
    #             for image, mask in zip(self.images, self.mask)
    #         ]

    #     else:
    #         raise RuntimeError

    # def flux_info(self, *args, **kwargs):

    #     if self.aper == None:

    #         raise TypeError

    #     else:

    #         if self.mask_member == "each":

    #             raise NotImplementedError

    #         elif self.mask_member == "all":
    #             infos = [
    #                 aperture_photometry(image, self.aper) for image in self.images()
    #             ]

    #         if not infos:
    #             raise RuntimeError

    #         bkgs = self.backgrounds(args, kwargs)
    #         for idx, info in enumerate(infos):
    #             info["background_flux"] = bkgs[idx]

    #     return infos


class DataCube:
    """
    This class handles cube data. It opens a local file or download a
    file from remote server by url.

    Attributes:
        name: the filename or url
    """

    def __init__(self, name: str) -> None:
        """
        Astronomical data cube (3D)

        :param name: the target filename or url
        :type name: str
        """
        if path.exists(name):
            self.name = name
        else:
            try:
                response = requests.get(name)
                self.name = name
            except requests.ConnectionError as exception:
                logging.ERROR("Could not find the file or URL")

        self._hdul = None
        # self.data = None
        # self.header = None

    @property
    def data(self):
        """
        the raw data of the target file

        :return: raw data
        :rtype: array
        """
        return self.hdul[0].data

    @property
    def header(self):
        """
        the header of the target file

        :return: header
        :rtype: HDU header
        """
        return self.hdul[0].header

    @property
    def hdul(self):
        """
        The HDUList of the file content

        :return: file content in HDUList
        :rtype: HDUList
        """
        if self._hdul:
            return self._hdul
        else:
            self.query()
            return self._hdul

    def query(self) -> object:
        """
        open the file from local or remote database
        """
        self._hdul = fits.open(self.name)
        # self.data = hdul[0].data
        # self.header = hdul[0].header

    def save(self, filename: str) -> None:
        """
        save file to the assign filename

        :param filename: the target filename
        :type filename: str
        """
        if not path.exist(filename):
            self.hdul.writeto(filename)
        else:
            print("Failed! File exists.")

    @property
    def spectrum(self):
        """
        get the mean spectrum of the cube

        :return: line spectrum
        :rtype: array
        """
        spectrum = np.nanmean(self.data, axis=(1, 2))

        # cube = SpectralCube(data=self.data, wcs=WCS(self.header))
        # cube_spectrum = cube.mean(axis=(1, 2))
        return spectrum

        # name = self.filename
        # cube = SpectralCube.read(self.name)
        # cube_spectrum = cube.mean(axis=(1, 2))
        # moment_0_filename = name + "_moment_0.fits"
        # moment_1_filename = name + "_moment_1.fits"
        # spec_file = self.filename + "_spectrum.jpg"
        # file_fits = fits.open(name)
        # spectral_data = file_fits[0].data
        # line_header = file_fits[0].header
        # moment_0_test = cube.moment(order=0)
        # moment_1_test = cube.moment(order=1)
        # moment_0_test.write(moment_0_filename)
        # moment_1_test.write(moment_1_filename)
        # vaxis = "CRVAL" + str(self.vindex)
        # pix_axis = "CDELT" + str(self.vindex)
        # velocity_axis_end = (
        #     line_header[vaxis] + len(cube_spectrum) * line_header[pix_axis]
        # ) / 1000
        # v_axes = np.linspace(
        #     line_header[vaxis] / 1000, velocity_axis_end, len(cube_spectrum)
        # )
        # plt.ioff()
        # plt.figure()
        # plt.step(v_axes, cube_spectrum)
        # plt.xlabel("Velocity [km/s]", fontsize=15)
        # plt.ylabel("Intensity [K]", fontsize=15)
        # plt.savefig(spec_file)
        # return cube_spectrum
