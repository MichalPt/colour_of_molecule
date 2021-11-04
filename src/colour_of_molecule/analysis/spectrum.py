
def abslines_to_molar_abs(inputf, title="Plot1", show_plot=False, stdev=3099.6, wav_range=(200, 1000), normalize=False):
    import numpy as np
    import matplotlib.pyplot as plt
    import colour

    abslines = inputf
    wav = list(map(lambda k: k.wavelength, abslines))
    fs = list(map(lambda l: l.oscillator_strength, abslines))

    start = wav_range[0]
    finish = wav_range[1]
    points = finish - start  # i.e. 1 nm

    # A sqrt(2) * standard deviation of 0.4 eV is 3099.6 nm. 0.1 eV is 12398.4 nm. 0.2 eV is 6199.2 nm, 0.33 eV = 3723.01 nm
    bands = wav
    f = fs

    # Basic check that we have the same number of bands and oscillator strengths
    if len(bands) != len(f):
        raise Exception('ERROR:   Number of bands does not match the number of oscillator strengths.')

    def gauss_band(x, band, strength, stdev):
        bandshape = 1.3062974e8 * (strength / (1e7 / stdev)) * np.exp(
            -(((1.0 / x) - (1.0 / band)) / (1.0 / stdev)) ** 2)
        return bandshape

    x = np.linspace(start, finish, points)

    composite = 0
    for i in range(0, len(bands), 1):
        peak = gauss_band(x, float(bands[i]), float(f[i]), stdev)
        composite += peak

    if show_plot == True:
        figg, axx = plt.subplots()
        axx.plot(x, composite)
        plt.xlabel('$\lambda$ / nm')
        plt.ylabel('$\epsilon$ / L mol$^{-1}$ cm$^{-1}$')
        plt.show()

    if normalize is True:
        from colour_of_molecule.analysis.common_tools import normalize_list
        composite = normalize_list(composite)

    data = colour.SpectralDistribution(data=composite, domain=x, name=title)

    return data


def molar_abs_to_complement_abs(spectrum, OD=0.15, normalize=True):
    import colour
    import numpy as np

    val = spectrum.values
    wav = spectrum.wavelengths
    tit = spectrum.name
    export = list()
    top = max(spectrum.values)

    norm = OD / top if normalize is True else 1

    for j in range(0, len(val), 1):
        expo = val[j] * norm
        if expo < 1e-10:
            expo = 1e-10
        exval = -np.log10(1 - 10 ** (-expo))
        export.append(exval)
    out = colour.SpectralDistribution(data=export, domain=wav, name=tit)
    return out


def find_colour(spectrum, col_map_f='CIE 1931 2 Degree Standard Observer'):

    import colour

    cmfs = colour.MSDS_CMFS[col_map_f]
    illuminant = colour.SDS_ILLUMINANTS['D65']

    wavs = spectrum.wavelengths

    try:
        XYZ = colour.sd_to_XYZ(spectrum, cmfs, illuminant)
    except:
        XYZ = colour.sd_to_XYZ(spectrum.interpolate(colour.SpectralShape(min(wavs), max(wavs), 1)), cmfs, illuminant)

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


