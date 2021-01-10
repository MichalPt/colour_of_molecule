def Gauss_log_to_abslines(log_file, title={},show_abs_lines_data=False,show_par=True, wav_shift=0.0):
    import re
    import os

    file = open(log_file, "r")
    
    if title=={}:
        title = os.path.basename(file.name)
    title = title.replace("_","-")
    
    list_file = list()
    fs = list()
    wav = list()
    
    parameters = list()

    reg_exc = re.compile("Excited State")
    reg_wav = re.compile("\s[\d.]*\snm")
    reg_f = re.compile("f\=[\d.]*\s")
    reg_par = re.compile("^(\s?#)p\s")
    reg_line = re.compile("^(\s?-)-*")
    reg_charge = re.compile("^\s?Charge.*\sMultiplicity")
    
    ki=0
    for x in file:
        #print(x)
        if reg_exc.search(x) != None:
            list_file.append(x)
        if ki==1:
            if reg_line.search(x) == None:
                parameters.append(x)
            else:
                ki = 0
        if reg_par.search(x) != None:
            parameters.append(x)
            ki = 1
        if reg_charge.search(x) != None:
            ch_mult = x[1:]
    file.close()
    
    header=""
    for q in parameters:
        header = header + q[1:len(q)-1]
    
    if show_par:
        print("Parameters:\n  ",header,"\n  ",ch_mult)

    for y in list_file:
        loc_wav = reg_wav.search(y)
        loc_f = reg_f.search(y)
        wav.append(float(y[loc_wav.start()+1 : loc_wav.end()-3])+wav_shift)
        fs.append(float(y[loc_f.start()+2 : loc_f.end()-1]))
        
    if wav==[]:
        print("Error in file import. Check the encoding of .txt file and eventually change it to ANSI.")
        
    if show_abs_lines_data:
        if wav_shift!=0: 
            print("Shift applied to wavelengths:\n  ",wav_shift," nm")
        print("Wavelenths:\n  ",wav)
        print("f:\n  ",fs)
    
    output = (title,wav,fs)
    
    return output


def abslines_to_molar_abs(input, show_plot=False, stdev=3099.6, wav_range=(200,1000)):
    import numpy as np
    import matplotlib.pyplot as plt
    import colour
    
    title,wav,fs = input
    
    start=wav_range[0]
    finish=wav_range[1]
    points=finish-start #i.e. 1 nm

    # A sqrt(2) * standard deviation of 0.4 eV is 3099.6 nm. 0.1 eV is 12398.4 nm. 0.2 eV is 6199.2 nm, 0.33 eV = 3723.01 nm
    bands = wav
    f = fs

    # Basic check that we have the same number of bands and oscillator strengths
    if len(bands) != len(f):
        print('Number of bands does not match the number of oscillator strengths.')
        sys.exit()

    def gaussBand(x, band, strength, stdev):
        bandshape = 1.3062974e8 * (strength / (1e7/stdev))  * np.exp(-(((1.0/x)-(1.0/band))/(1.0/stdev))**2)
        return bandshape

    x = np.linspace(start,finish,points)

    composite = 0
    for i in range(0,len(bands),1):
        peak = gaussBand(x, float(bands[i]), float(f[i]), stdev)
        composite += peak
    
    if show_plot==True:
        figg, axx = plt.subplots()
        axx.plot(x,composite)
        plt.xlabel('$\lambda$ / nm')
        plt.ylabel('$\epsilon$ / L mol$^{-1}$ cm$^{-1}$')
        plt.show()

    data = colour.SpectralDistribution(data=composite,domain=x,name=title)
    
    return data


def txt_spectrum_to_abs(text_file,title={},wav_range=(200,1000),from_Gaussian=False,
                        header=True, col_wav=1, col_val=2, show_abs_data=False):
    import re
    import os
    import colour
    import numpy as np

    txt_file = open(text_file, "r")
    
    if title=={}:
        title = os.path.basename(txt_file.name)
    title = title.replace("_","-")
    
    val = list()
    wav = list()

    reg_spectra = re.compile("# Spectra")
    reg_num = re.compile("[\d.]+")
    
    if from_Gaussian:
        ti=0
        ri=0
    elif header:
        ti=1
        ri=-1
    else:
        ti=1
        ri=-2
    
    for k,x in enumerate(txt_file):
        if reg_spectra.search(x) != None and from_Gaussian:
            ti=1
            ri=k
        if ti==1 and k>ri+1:
            found = reg_num.findall(x)
            #print(found)
            v = float(found[col_val-1])
            w = float(found[col_wav-1])
            #if w >= wav_range[0] and w <= wav_range[1]:
            val.append(v)
            wav.append(w)
    txt_file.close()

    if show_abs_data:
        print("Wavelenths:\n  ",wav)
        print("Abs:\n  ",val)
    #print(wav)
    try:
        data0 = colour.SpectralDistribution(data=val,domain=wav,name=title)
    except:
        print("Error while creating SD. Imported data might not be well formated or ordered.")
    data = data0.interpolate(colour.SpectralShape(np.floor(wav_range[0]),np.floor(wav_range[1]),1))
    #print(data)
    return data


def molar_abs_to_complement_abs(spectrum,OD=0.15,renormalize=True):
    import colour
    import numpy as np
    #numpy.seterr(divide='ignore') 
    
    val = spectrum.values
    wav = spectrum.wavelengths
    tit = spectrum.name
    export = list()
    top = max(spectrum.values)
    
    norm = OD/top if renormalize==True else 1
    
    for j in range(0,len(val),1):
        expo = val[j]*norm
        if expo < 1e-10:
            expo = 1e-10
        exval = -np.log10(1-10**(-expo))
        export.append(exval)
        #print((wav[j],val[j],exval))
    out = colour.SpectralDistribution(data=export,domain=wav,name=tit)
    #print(out)
    return out


def find_colour(spectrum, col_map_f='CIE 1931 2 Degree Standard Observer'):
    from matplotlib import rc
    import numpy as np
    import colour
    import matplotlib.patches as mpatches
    from matplotlib import pyplot as plt
    from matplotlib.pyplot import figure
    
    cmfs = colour.MSDS_CMFS[col_map_f]
    illuminant = colour.SDS_ILLUMINANTS['D65']
    
    val = spectrum.values
    wavs = spectrum.wavelengths
    
    try:
        XYZ = colour.sd_to_XYZ(spectrum, cmfs, illuminant)
    except:
        XYZ = colour.sd_to_XYZ(spectrum.interpolate(colour.SpectralShape(min(wavs), max(wavs), 1)), cmfs, illuminant)
        
    RGB = colour.XYZ_to_sRGB(XYZ / 100)
    for i in range(0,3,1):
        if RGB[i]<0:
            RGB[i] = 0
        if RGB[i]>1:
            RGB[i]=1

    return RGB


def spectrum_colour_analysis(file,title={}, stdev = 3096.01,renormalize = True,OD = 0.15, show_plot = False,
                         col_map_f="CIE 1931 2 Degree Standard Observer", wav_range=(200,1000), fancy_cols=True,
                         show_abs_lines_data=False, show_par=True, from_Gaussian=False, header=True, 
                         col_wav=1, col_val=2, show_abs_data=False, give_raw_data=True, wav_shift=0.0):
    # available kwargs:
    #    title                   - set the title for the spectrum (imported file name is set by default)
    #    wav_range = (min,max)   - set the limit wavelength to be plotted
    #    stdev = 3096.01         - Gaussian width
    #    renormalize = True      - the molar abs. sprectrum is converted into an absorption spectrum of defined peak OD
    #    OD = 0.15               - target OD to convert the molar abs. spectrum into (i.e. target sample dilution)
    #    show_plot = False       - show raw molar abs. spectra
    #    show_par = True         - show calculation parameters imported from .log
    #    show_abs_lines_data = False  - show excitation wavelength and abs. dipole strengths data imported from .log
    #    col_map_f               - color mapping function (see colour-science package documentation), "CIE 1931 2 Degree Standard Observer" set by default
    #    fancy_cols = True       - display colours
    #    from_Gaussian = False    - imported .txt file has GaussView spectrum format
    #    header = True           - imported .txt file has a single row header
    #    col_wav = 1             - index (1,2,...) of column with wavelength data in imported .txt file
    #    col_val = 2             - index (1,2,...) of column with abs. spectrum data in imported .txt file
    #    give_raw_data = True    - return raw data (colour.SpectralDistribution) and rgb code as [sd1,sd2,RGB]
    #    wav_shift = 0           - artificially move with the calculated spectral values along the wavelengths to see a colour change
    
    
    import colour
    import numpy as np
    from matplotlib import pyplot as plt
    from matplotlib import rc
    from matplotlib.pyplot import figure
    from matplotlib import gridspec
    from matplotlib.patches import Polygon
    from colour.plotting import (XYZ_to_plotting_colourspace,  filter_cmfs, CONSTANTS_COLOUR_STYLE )
    from colour.colorimetry import (CCS_ILLUMINANTS, wavelength_to_XYZ)
    from colour.utilities import (first_item, normalise_maximum)
    
    ###
    import re
    import os
    
    reg_exten = re.compile("\.\w{2,4}$")
    
    xfile = open(file, "r")
    file_name = os.path.basename(xfile.name)
    xfile.close()
    extension = reg_exten.findall(file_name)[0]
    
    if extension==".log":
        lines = Gauss_log_to_abslines(file,title=title,show_abs_lines_data=show_abs_lines_data,
                                      wav_shift=wav_shift,show_par=show_par)
        spec1 = abslines_to_molar_abs(lines, stdev=stdev, show_plot=show_plot,wav_range=wav_range)
    elif extension==".txt":
        spec1 = txt_spectrum_to_abs(file,title=title,wav_range=wav_range,from_Gaussian=from_Gaussian,
                        header=header, col_wav=col_wav, col_val=col_val, show_abs_data=show_abs_data)
    
    if title=={}:
        title=spec1.name
    
    spec2 = molar_abs_to_complement_abs(spec1, renormalize=renormalize,OD=OD)
    RGB = find_colour(spec2, col_map_f=col_map_f)
    
    output=[spec1,spec2,RGB]

    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size': 13})
    rc('text', usetex=True)
    rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{mhchem, physics} \usepackage[utf8]{inputenc} \usepackage{textcomp}' )

    plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = True
    plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = False

    fig, axes = plt.subplots(num=None, figsize=(11, 4), dpi=400, facecolor='w', edgecolor='k')

    gs = gridspec.GridSpec(ncols=3, nrows=1, width_ratios=[5, 5, 1], figure=fig, wspace=0.25) 

    for i in [0,1]:
        ax = plt.subplot(gs[i])
        wavelengths = output[i].wavelengths
        values = output[i].values
        ax.plot(wavelengths,values, linewidth=1.2, color='k')
        ax.set_xlabel(r'$\mathrm{vlnov\acute{a} \: d\acute{e}lka} \: \lambda \: [\mathrm{nm}]$')
        ax.set_xlim(min(wavelengths),max(wavelengths))
        ax.set_ylim(0, max(values))
        plt.locator_params(axis='y', nbins=5)
        ax.ticklabel_format(axis="y",style="sci",scilimits=(-1,1))
        if i==0:
            #ax.set_title(label=f'$\\mathrm{{{output[i].name}}}$')
            ax.set_title(label=output[i].name)
            ax.set_ylabel(r'$ \epsilon \: [\mathrm{dm}^3 \cdot \mathrm{mol}^{-1} \cdot \mathrm{cm}^{-1}]$')
        elif i==1: 
            #ax.set_title(label=f'$\\mathrm{{{output[i].name} \: (K)}}$')
            ax.set_title(label=output[i].name + "\: $\mathrm{(K)}$")
            ax.set_ylabel(r'$ "1-A" \: [\mathrm{a.u.}]$')
        
        if fancy_cols==True:
            cmfs = first_item(filter_cmfs(col_map_f).values())
            wlen_cmfs = [n for n in wavelengths if n > cmfs.shape.start and n < cmfs.shape.end]

            #global clr
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
            ax.add_patch(polygon)

            padding = 0.1
            
            for dom,col in [(wavelengths-padding,'black'),(wlen_cmfs,clr)]:
                ax.bar(
                x=dom,
                height=max(values),
                width=1 + padding,
                color=col,
                align='edge',
                clip_path=polygon)

    ax2 = plt.subplot(gs[2])

    ax2.set_facecolor(output[2])
    ax2.set_ylabel(f'$\\mathrm{{RGB}}({["%.3f" % i for i in output[2]]})$')
    plt.setp((ax2.get_yticklabels(), ax2.get_yticklines(),ax2.get_xticklabels(),ax2.get_xticklines()), visible=False)
    
    plt.show()
    
    if give_raw_data:
        return output

    
def plot_multiple(*inputs,wav_range=(300,1000),plot_kwargs={'linewidth':1},plot_kwargs_of_nth={}):
    import colour
    import numpy as np
    from matplotlib import pyplot as plt
    from matplotlib import rc
    from matplotlib.pyplot import figure
    from matplotlib import gridspec
    from matplotlib.patches import Polygon
    
    # Deactivating the MatplotlibDeprecationWarning ...
    import warnings
    import matplotlib.cbook
    warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
    #
    
    fig = plt.figure(figsize=(11,4), dpi=400, facecolor='w', edgecolor='k')
    gs0 = fig.add_gridspec(nrows=1, ncols=3, width_ratios=[6,6,1])    #, height_ratios=[1/len(inputs)]*len(inputs))
    
    gs00 = gs0[2].subgridspec(nrows=len(inputs), ncols=1,wspace=0.25)
    
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size': 13})
    rc('text', usetex=True)
    rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{mhchem, physics} \usepackage[utf8]{inputenc} \usepackage{textcomp}' )
    
    w_range = list(wav_range)
    legend = list()
    ax2 = [None]*len(inputs)
    
    ax0 = plt.subplot(gs0[0])
    ax1 = plt.subplot(gs0[1])
    
    for i, data in enumerate(inputs):
        w_min = min(data[0].wavelengths)
        w_max = max(data[0].wavelengths)
        if w_min > w_range[0]:
            w_range[0] = w_min
        if w_max < w_range[1]:
            w_range[1] = w_max
        legend.append("$\mathrm{("+str(i+1)+")}$\:"+data[0].name)
    for i, data in enumerate(inputs):
        wavelengths0 = data[0].wavelengths
        norm0=max([data[0].values[n] for n,m in enumerate(wavelengths0) if m<=w_range[1] and m>=w_range[0]])        
        values0 = data[0].values/norm0
        
        wavelengths1 = data[1].wavelengths          
        norm1=max([data[1].values[n] for n,m in enumerate(wavelengths1) if m<=w_range[1] and m>=w_range[0]]) 
        values1 = data[1].values/norm1          
        
        try:
            kw = plot_kwargs_of_nth[i+1]
            ax0.plot(wavelengths0, values0, **plot_kwargs, **kw)
            ax1.plot(wavelengths1, values1, **plot_kwargs, **kw)
        except:
            ax0.plot(wavelengths0, values0, **plot_kwargs)
            ax1.plot(wavelengths1, values1, **plot_kwargs)
        
        ax2[i] = plt.subplot(gs00[i])
        ax2[i].set_xlim(0,1)
        ax2[i].set_ylim(0,1)
        plt.setp((ax2[i].get_yticklabels(), ax2[i].get_yticklines(),ax2[i].get_xticklabels(),
                  ax2[i].get_xticklines()), visible=False)
        ax2[i].set_facecolor(data[2])
        ax2[i].text(0.52, 0.5, f'$({i+1})$', va='center', ha='center')
    
    ax0.set_xlim(w_range)
    ax0.set_ylim(0,1)
    ax0.set_yticks([0.0,0.5,1.0])
    ax0.set_xlabel(r'$\mathrm{vlnov\acute{a} \: d\acute{e}lka} \: \lambda \: [\mathrm{nm}]$')
    ax1.set_xlim(w_range)
    ax1.set_ylim(0,1)
    ax1.set_yticks([0.0,0.5,1.0])
    ax1.set_xlabel(r'$\mathrm{vlnov\acute{a} \: d\acute{e}lka} \: \lambda \: [\mathrm{nm}]$')
    
    ax0.set_title("$\\boldsymbol{\\mathrm{A}}$")
    ax1.set_title("$\\boldsymbol{\\mathrm{B}}$")
    ax2[0].set_title("$\\boldsymbol{\\mathrm{C}}$")
    
    numcol = 2 if len(inputs)==2 or len(inputs)==4 else 3
    ax1.legend(labels=legend, bbox_to_anchor=(0.05, -0.18), loc='upper center', ncol=numcol, 
           frameon=False, columnspacing=2)
    
    plt.subplots_adjust(hspace=.0)

    plt.show()
