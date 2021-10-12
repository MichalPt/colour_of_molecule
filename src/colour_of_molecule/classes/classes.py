import os

class AbsLine:
    def __init__(self, wavelength, strength, *transitions):
        self.wavelength = wavelength
        self.oscillator_strength = strength
        self.transitions = list(transitions)


class MolarAbsSpectrum:
    def __get__(self, instance, owner):
        from colour_of_molecule.analysis.spectrum import abslines_to_molar_abs

        data = abslines_to_molar_abs(instance.abs_lines,
                                     stdev=instance.standard_deviation,
                                     wav_range=instance.wavelength_range,
                                     normalize=instance.normalize_absorption_spectrum)
        return data


class ComplementaryAbsSpectrum:
    def __get__(self, instance, owner):
        from colour_of_molecule.analysis.spectrum import molar_abs_to_complement_abs

        data = molar_abs_to_complement_abs(instance.molar_abs_spectrum,
                                           OD=instance.optical_density,
                                           normalize=instance.normalize_complementary_spectrum)
        return data


class ColourRGB:
    def __get__(self, instance, owner):
        from colour_of_molecule.analysis.spectrum import find_colour
        rgb = find_colour(instance.complementary_abs_spectrum)
        return rgb


class Energy:
    __default_unit = "kcal/mol"

    def convert_nm_to_rcm(self, x):
        return 1e7 / x

    supported_units = {  # 1 unit = ... kcal/mol
        "kcal/mol": 1,
        "Hartree": 630,
        "au": 630,
        "eV": 23.060541945329334,
        "kJ/mol": 0.2390057361376673,
        "Ry": 315,
        "cm-1": 0.0028591441663694465,
        "nm": convert_nm_to_rcm,
    }

    def __init__(self, value, units=__default_unit):
        if units is self.__default_unit:
            self.value = float(value)
            self.units = units
        elif units in self.supported_units:
            if callable(self.supported_units.get(units)) and units == "nm":
                func = self.supported_units.get(units)
                coeff = self.supported_units.get("cm-1") / self.supported_units.get(self.__default_unit)
                self.value = self.func(value) * coeff
            else:
                coeff = self.supported_units.get(units) / self.supported_units.get(self.__default_unit)
                self.value = float(value) * coeff
                self.units = self.__default_unit
        else:
            raise Exception("ERROR:\tUnsupported unit was encountered. Only these units are currently supported: " +
                            "".join([str(i) + "  " for i in self.supported_units.keys()]))

    def __sub__(self, other):
        result = self.value - other.value
        return Energy(result)

    def __add__(self, other):
        result = self.value + other.value
        return Energy(result)

    def in_units(self, new_units):
        if new_units == "nm":
            out = self.convert_nm_to_rcm(self.value / self.supported_units.get("cm-1") * self.supported_units.get(self.__default_unit))
            return round(out, 3)
        elif new_units in self.supported_units:
            out = self.value / self.supported_units.get(new_units) * self.supported_units.get(self.__default_unit)
            return out
        else:
            raise Exception("ERROR:\tUnsupported unit was encountered. Only these units are currently supported: " +
                            "".join([str(i) + "  " for i in self.supported_units.keys()]))


class OpList(list):
    def __add__(self, other):
        return OpList([x + y for (x, y) in zip(self, other)])

    def __sub__(self, other):
        return OpList([x - y for (x, y) in zip(self, other)])

    def __neg__(self):
        return OpList([-x for x in self])

    def get_abs_length(self):
        return int(sum([abs(x) for x in self]))

    def as_lengths(self):
        return OpList([len(x) for x in self])


class FontSettings():

    def __init__(self, newfonts=list(), newsizes=list(), use_all=False):
        self.fontdict = {'title': 'Calibri',
                    'axis': 'Calibri',
                    'axis_tick_labels': 'Calibri',
                    'legend':'Calibri',
                    'all':'Calibri'}
        self.sizedict = {'title': 14, 'axis': 12, 'axis_tick_labels': 12, 'legend': 12, 'all': 12}
        self.fontdict.update(newfonts)
        self.sizedict.update(newsizes)

        if use_all is True:
            self.fontdict = {key:self.fontdict['all'] for key in self.fontdict.keys()}

        import matplotlib.font_manager as font_manager
        self.fonts = {key:font_manager.FontProperties(family=self.fontdict[key], weight='normal', style='normal',
                                                      size=self.sizedict[key])
                 for key in self.fontdict.keys()}


class File:
    supported_formats = {"gaussian": (0, "Entering Link 1"),
                         "orca": (3, "* O   R   C   A *"),
                         "mndo": (0, "PROGRAM MNDO"),
                         "molpro": (0, "***  PROGRAM SYSTEM MOLPRO  ***")
                         }

    def __init__(self, path):
        def check_line(dict, line):
            ans = [k for k in list((i, j[0]) if j[1] in line else False for i, j in dict.items()) if k]
            if ans:
                return ans[0]

        def sanity_check_type_module_class(self):
            if not hasattr(self, "type"):
                raise Exception("ERROR:\tAttribute \"type\" not found.")
            else:
                from importlib.util import find_spec

                mod_path = "colour_of_molecule.input." + self.type
                class_name = "File_" + self.type

                check_module = find_spec(mod_path)
                if check_module is None:
                    raise Exception("ERROR:\tModule \"" + mod_path + "\" not found. Is this file type implemented yet?")
                else:
                    class_check = hasattr(import_module(mod_path), class_name)
                    if class_check is not True:
                        raise Exception("ERROR:\tClass \"" + class_name + "\" not found in module \"" + mod_path + "\"")
            #print("INFO:\tAll sanity checks passed successfully.")
            pass

        self.path = path
        self.filename = os.path.basename(path)
        self.ranges_of_comps = dict()
        self.number_of_comps = 0
        self.more_than_one_comp = False

        self.standard_deviation = 3096.01
        self.optical_density = 0.15
        self.wavelength_range = (200, 800)
        self.transition_minimal_amplitude = 0.5
        self.normalize_absorption_spectrum = False
        self.normalize_complementary_spectrum = True

        self.plot_title = ""
        self.legend_title = ""

        with open(path, "r") as file:

            self.name = os.path.basename(file.name).replace("_", "-")

            for index, line in enumerate(file):
                out = check_line(self.supported_formats, line)
                if out:
                    self.type = out[0]
                    self.ranges_of_comps.update({self.number_of_comps: (index - out[1], None)})
                    if self.number_of_comps > 0:
                        self.more_than_one_comp = True
                        self.ranges_of_comps.update({self.number_of_comps -1 : (self.ranges_of_comps.get(self.number_of_comps-1)[0], index - out[1] - 1)})
                    self.number_of_comps += 1

            if self.type is None:
                raise Exception("ERROR:\tInput file type is not implemented.\n\tOnly these types are currently supported:  " +
                                "".join([str(i).upper()+"  " for i in self.supported_formats.keys()]))
            print("INFO:\tNumber of recognised computations in "+self.type.capitalize()+" file \""+self.name+"\" is:   "+str(self.number_of_comps))

            ### Switch for multiple files should be here!!!  ###

            from importlib import import_module

            sanity_check_type_module_class(self)

            print("INFO:\tChanging class to \"File_"+self.type+"\"")
            clss = getattr(import_module("colour_of_molecule.input."+self.type), "File_"+self.type)
            self.__class__ = clss

    molar_abs_spectrum = MolarAbsSpectrum()
    complementary_abs_spectrum = ComplementaryAbsSpectrum()
    colour_rgb = ColourRGB()




