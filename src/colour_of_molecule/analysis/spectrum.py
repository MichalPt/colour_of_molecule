import numpy as np
import matplotlib.pyplot as plt
import colour
from colour_of_molecule.classes.classes import EnergyRange, Energy, EnergyAxis
from colour_of_molecule.utils.energy_units import get_current_energy_units
from colour_of_molecule.plotting.lineshapes import gaussian_profile

def abslines_to_molar_abs(file, title="Plot1", show_plot=False):
    abslines = file.abs_lines
    plot_range = file._plot_range
    normalize = file.normalize_absorption_spectrum
    npoints = file.npoints
    bands = list(map(lambda k: k.energy, abslines))
    fs = list(map(lambda l: l.oscillator_strength, abslines))

    # Basic check that we have the same number of bands and oscillator strengths
    if len(bands) != len(fs):
        raise Exception('ERROR:   Number of bands does not match the number of oscillator strengths.')

    x = plot_range.energy_axis
    converted_x = x.in_current_units()
    composite = np.zeros_like(x.value)

    for i, band in enumerate(bands):
        # lineshape always calculated in linear energy units (kcal/mol - default unit)
        peak = gaussian_profile(x.value,
                                float(band.value), 
                                float(fs[i]), 
                                file._fwhm.value)
        composite += peak

    if show_plot == True:
        figg, axx = plt.subplots()
        axx.plot(x.value, composite)
        plt.xlabel('$E$ / {}'.format(get_current_energy_units()))
        plt.ylabel('$\epsilon$ / L mol$^{-1}$ cm$^{-1}$')
        plt.show()

    if normalize is True:
        from colour_of_molecule.analysis.common_tools import normalize_list
        composite = np.array(normalize_list(composite))

    x_sorted = np.sort(converted_x)
    y_sorted = composite[np.argsort(converted_x)]

    data = colour.SpectralDistribution(data=y_sorted, domain=x_sorted, name=title)

    return data


def molar_abs_to_transmittance(spectrum, OD=0.15):
    import colour
    import numpy as np

    val = spectrum.values
    wav = spectrum.wavelengths
    title = spectrum.name
    export = list()

    for j in range(0, len(val), 1):
        A = val[j]
        T = 10 ** (-A * OD)
        export.append(T)
    out = colour.SpectralDistribution(data=export, domain=wav, name=title)
    return out


def find_colour(spectrum):
    import colour

    col_map_f='CIE 1931 2 Degree Standard Observer'
    cmfs = colour.MSDS_CMFS[col_map_f]
    illuminant = colour.SDS_ILLUMINANTS['D65']

    energies = Energy(spectrum.wavelengths)
    wavelengths = np.array(energies.in_units("nm"), dtype=int)

    x_sorted = np.sort(wavelengths)
    y_sorted = spectrum.values[np.argsort(wavelengths)]

    wl_spectrum =colour.SpectralDistribution(data=y_sorted, domain=x_sorted, name=spectrum.name)
    XYZ = colour.sd_to_XYZ(wl_spectrum, cmfs, illuminant)

    # try:
    #     XYZ = colour.sd_to_XYZ(wl_spectrum, cmfs, illuminant)
    # except:
    #     XYZ = colour.sd_to_XYZ(spectrum.interpolate(colour.SpectralShape(min(wavelengths), max(wavelengths), 1)), cmfs, illuminant)

    RGB = colour.XYZ_to_sRGB(XYZ / 100)

    for i in range(0, 3, 1):
        if RGB[i] < 0:
            RGB[i] = 0
        if RGB[i] > 1:
            RGB[i] = 1

    return RGB


def find_colour_single(wl):
    import colour
    from colour.colorimetry import wavelength_to_XYZ

    if wl < 360 or wl > 830:
        RGB = (0, 0, 0)
    else:
        XYZ = wavelength_to_XYZ(wl)
        RGB = colour.XYZ_to_sRGB(XYZ)
        for i in range(0, 3, 1):
            if RGB[i] < 0:
                RGB[i] = 0
            if RGB[i] > 1:
                RGB[i] = 1
    return (RGB)


