import numpy as np

def gaussian_profile(x, band, strength, fwhm):
    stdev = fwhm / (2 * np.sqrt(2 * np.log(2)))
    bandshape = 9.2369174e7 * strength / stdev * np.exp(-((x - band)**2 / (2 * stdev**2)))
    return bandshape

