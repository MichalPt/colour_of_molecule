from ..classes.classes import FontSettings
import numpy as np

def plot_single_spectrum(file, save="", size=(6,3), dpi=400, rainbow=True,
                         title="", fonts=FontSettings(),
                         xaxis_label="wavelength [nm]", yaxis_label="relative absorbance",
                         lines_show=True, lines_ratio=(14, 1), lines_colours=True, lines_width=1.2, lines_lim=0.0001):
    from matplotlib import rcParams
    from matplotlib import pyplot as plt
    from colour_of_molecule.classes.classes import File
    from colour_of_molecule.analysis.spectrum import find_colour_single

    #sanity check:
    if not isinstance(file, File):
        raise Exception("ERROR:\tUnrecognized input. First argument has to be derived from class \"File\".")

    if title != "":
        file.plot_title = title

    wavelengths = file.molar_abs_spectrum.wavelengths
    values = file.molar_abs_spectrum.values
    wav_range = file.wavelength_range

    rcParams['font.size'] = fonts.sizedict['all']
    rcParams['font.family'] = fonts.fontdict['all']

    nrows = 1 if lines_show is False else 2
    gskw = dict(dict(height_ratios=list(lines_ratio), hspace=0.0)) if lines_show is True else None
    fig, ax = plt.subplots(nrows=nrows, figsize=size, dpi=dpi, facecolor='w', edgecolor='k',
                           gridspec_kw=gskw, sharex=True)

    boo = type(ax) == np.ndarray
    ax0 = ax[0] if boo is True else ax
    ax1 = ax[1] if boo is True else ax

    ax1.set_xlabel(xaxis_label)
    ax1.set_xlim(*wav_range)
    ax0.set_ylabel(yaxis_label)
    ax0.set_ylim(0, max(values))
    ax0.ticklabel_format(axis="y", style="sci", scilimits=(-1, 1))
    ax0.locator_params(axis='y', nbins=5)

    ax0.plot(wavelengths, values, linewidth=1.2, color='k')

    if lines_show is True:
        abslines = file.abs_lines
        for absline in abslines:
            wv = absline.wavelength
            if lines_colours is False:
                col = 'k'
            else:
                if absline.oscillator_strength > lines_lim:
                    col = find_colour_single(wv)
                else:
                    col = 'k'
            ax1.plot([wv, wv],[0, 1], color=col, linewidth=lines_width)
        plt.setp(ax1.get_yticklabels(), visible=False)
        plt.setp(ax1.get_yticklines(), visible=False)

    if rainbow is True:
        add_rainbow(ax0, wavelengths, values)

    for label in [*ax1.get_xticklabels(), *ax0.get_yticklabels()]:
        label.set_fontproperties(fonts.fonts['axis_tick_labels'])

    plot_title = file.plot_title if file.plot_title != "" else "File: {0}".format(file.filename)
    ax0.set_title(plot_title)

    plt.tight_layout()

    if save != "":
        if isinstance(save, str):
            from colour_of_molecule.analysis.common_tools import file_saver
            fpath = file_saver(save)
            if fpath != save:
                print("INFO:\tThe provided filename already exists. The file will be saved in \"{new}\" instead.".format(new=fpath))
            fig.savefig(fpath, dpi=dpi)
            print("INFO:\tImage saved successfully in \"{path}\".".format(path=fpath))

    plt.show()


def add_rainbow(axis, wavelengths, values, opacity=100):
    # sanity check:
    if not hasattr(axis, 'plot') and not hasattr(axis, 'add_patch'):
        raise Exception("ERROR:\tFirst argument needs to have method \"plot\".")

    from colour.plotting import XYZ_to_plotting_colourspace, filter_cmfs, CONSTANTS_COLOUR_STYLE
    from colour.colorimetry import CCS_ILLUMINANTS, wavelength_to_XYZ
    from colour.utilities import first_item, normalise_maximum
    from matplotlib.patches import Polygon

    col_map_f = "CIE 1931 2 Degree Standard Observer"

    cmfs = first_item(filter_cmfs(col_map_f).values())
    wlen_cmfs = [n for n in wavelengths if n > cmfs.shape.start and n < cmfs.shape.end]

    clr = XYZ_to_plotting_colourspace(
        wavelength_to_XYZ(wlen_cmfs, cmfs),
        CCS_ILLUMINANTS['CIE 1931 2 Degree Standard Observer']['E'],
        apply_cctf_encoding=False)

    clr = normalise_maximum(clr)
    clr = CONSTANTS_COLOUR_STYLE.colour.colourspace.cctf_encoding(clr)

    polygon = Polygon(
        np.vstack([
            [min(wavelengths), 0],
            np.array([wavelengths, values]).T.tolist(),
            [max(wavelengths), 0],
        ]),
        facecolor='none',
        edgecolor='none')
    axis.add_patch(polygon)

    if opacity < 100:
        padding = 0
    else:
        padding = 0.1

    for dom, col in [(wavelengths - padding, 'black'), (wlen_cmfs, clr)]:
        axis.bar(
            x=dom,
            height=max(values),
            width=1 + padding,
            color=col,
            align='edge',
            alpha=opacity/100,
            clip_path=polygon
        )

    pass


def create_label(list1, list2):
    if type(list1) != list:
        list1 = [list1]
    if type(list2) != list:
        list2 = [list2]

    def equal_length(list1, list2):
        len1 = len(list1)
        len2 = len(list2)
        dif = len1 - len2
        if dif < 0:
            list1.extend([""] * abs(dif))
        elif dif > 0:
            list2.extend([""] * abs(dif))
        return list1, list2

    list1, list2 = equal_length(list1, list2)

    len1 = max([len(i) for i in list1])
    len2 = max([len(i) for i in list2])

    part = '{0: {align}{length}}'

    lines = str()
    for i, (m, n) in enumerate(zip(list1, list2)):
        if i == 0:
            line = part.format(" $\\rightarrow$ ", align='^', length=4).join(
                (part.format(m, align='>', length=len1), part.format(n, align='<', length=len2)))
            lines = line
        else:
            line = part.format("", align='^', length=4).join(
                (part.format(m, align='>', length=len1), part.format(n, align='<', length=len2)))
            lines = lines + "\n" + line

    return lines


def plot_abs_lines(file, save="", size=(7,7), dpi=200, fonts=FontSettings(),
                   title="",
                   xaxis_label="wavelength [nm]", yaxis_label="$\epsilon$", yaxis_label_right="oscillator strength"):
    from colour_of_molecule.classes.classes import File
    from matplotlib import pyplot as plt
    from matplotlib import rcParams

    rcParams['font.size'] = fonts.sizedict['all']
    rcParams['font.family'] = fonts.fontdict['all']

    rcParams['xtick.bottom'] = rcParams['xtick.labelbottom'] = True
    rcParams['xtick.top'] = rcParams['xtick.labeltop'] = False

    # sanity check:
    if not isinstance(file, File):
        raise Exception("ERROR:\tFunction argument has to be of \"File\" class.")

    if title != "":
        file.plot_title = title

    ab_spectrum = file.molar_abs_spectrum
    wavelengths = ab_spectrum.wavelengths
    values = ab_spectrum.values
    wav_range = file.wavelength_range
    abs_lines = file.abs_lines

    fig2, axis1 = plt.subplots(dpi=dpi, figsize=size, facecolor='w', edgecolor='k')
    axis1.set_xlim(wav_range)
    axis2 = axis1.twinx()
    axis2.set_xlim(wav_range)

    axis1.locator_params(axis='y', nbins=5)
    axis1.ticklabel_format(axis="y", style="sci", scilimits=(-1, 1))

    axis1.plot(wavelengths, values, color='black', alpha=0.3, linewidth=0.6)

    def check_transition_amplitudes(absl):
        ls = [ab[2] for ab in absl if len(ab) >= 3]
        maximal = max(ls) if len(ls) > 0 else None
        return maximal

    for ab in abs_lines:
        text_label = ""
        maximal = check_transition_amplitudes(ab.transitions[0])
        #print("   ABL: ", ab.transitions[0])
        for args in ab.transitions[0]:

            length = len(args)
            if length >= 3:
                amplitude = args[2]
                c = " ({:.2f})".format(amplitude)
            else:
                amplitude = None
                c = ""
            a, b = args[:2]
            lab = create_label(a, b)
            #print("Lab: ", lab, c, amplitude, maximal, amplitude >= file.transition_minimal_amplitude or amplitude == maximal)
            if amplitude is not None:
                text_label += lab + c + "\n" if amplitude >= file.transition_minimal_amplitude or amplitude == maximal else ""
            else:
                text_label += lab + "\n"
        text_label = text_label[:-1]

        if ab.wavelength > wav_range[0] and ab.wavelength < wav_range[1]:
            axis2.plot([ab.wavelength, ab.wavelength], [0, ab.oscillator_strength], label=text_label, linewidth=2)

    axis1.set_ylabel(yaxis_label)
    axis1.set_xlabel(xaxis_label)
    axis2.set_ylabel(yaxis_label_right)

    for ax in [axis1, axis2]:
        for label in [*ax.get_xticklabels(), *ax.get_yticklabels()]:
            label.set_fontproperties(fonts.fonts['axis_tick_labels'])

    plot_title = file.plot_title if file.plot_title == "" else "File: {0}".format(file.filename)
    legend_title = file.legend_title if file.legend_title != "" \
        else "Electron transitions with amplitude greater than {:.2f} contributing to absorption lines:".format(file.transition_minimal_amplitude)

    plt.title(plot_title)
    axis2.legend(bbox_to_anchor=(0.5, -0.18), loc='upper center', ncol=3,
                 title=legend_title,
                 frameon=False, columnspacing=2, labelspacing=1.5, prop=fonts.fonts['legend'])

    plt.tight_layout()

    if save != "":
        if isinstance(save, str):
            from colour_of_molecule.analysis.common_tools import file_saver
            fpath = file_saver(save)
            if fpath != save:
                print("INFO:\tThe provided filename already exists. The file will be saved in \"{new}\" instead.".format(new=fpath))
            fig2.savefig(fpath, dpi=dpi)
            print("INFO:\tImage saved successfully in \"{path}\".".format(path=fpath))
    plt.show()


def get_colour(file, save="", size=(3,3), dpi=200, fonts=FontSettings(), title="RGB: {RGB}",
               col_map_f='CIE 1931 2 Degree Standard Observer'):
    from colour_of_molecule.analysis.spectrum import find_colour, molar_abs_to_complement_abs
    from matplotlib import rcParams
    from matplotlib import pyplot as plt
    from colour_of_molecule.classes.classes import File

    # sanity check:
    if not isinstance(file, File):
        raise Exception("ERROR:\tUnrecognized input. First argument has to be derived from class \"File\".")

    if title != "":
        file.plot_title = title

    RGB = find_colour(file.complementary_abs_spectrum, col_map_f)

    print("RGB:\t{0}".format(RGB))

    rcParams['font.size'] = fonts.sizedict['all']
    rcParams['font.family'] = fonts.fontdict['all']
    fig, axes = plt.subplots(num=None, figsize=size, dpi=dpi, facecolor='w', edgecolor='k')

    axes.set_facecolor(RGB)
    axes.set_xlabel(title.format(RGB=", ".join(["%.3f" % i for i in RGB])), fontsize=fonts.sizedict['title'])
    plt.setp((axes.get_yticklabels(), axes.get_yticklines(), axes.get_xticklabels(), axes.get_xticklines()), visible=False)

    if save != "":
        if isinstance(save, str):
            from colour_of_molecule.analysis.common_tools import file_saver
            fpath = file_saver(save)
            if fpath != save:
                print("INFO:\tThe provided filename already exists. The file will be saved in \"{new}\" instead.".format(new=fpath))
            fig.savefig(fpath, dpi=dpi)
            print("INFO:\tImage saved successfully in \"{path}\".".format(path=fpath))

    plt.show()










