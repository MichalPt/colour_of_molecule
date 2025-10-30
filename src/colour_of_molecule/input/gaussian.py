import re
from colour_of_molecule.classes import AbsLine, File, Energy
from colour_of_molecule.analysis.common_tools import homo_lumo
from colour_of_molecule.utils.energy_units import get_current_energy_units


class File_gaussian(File):
    class AbsLines_getter:
        def __get__(self, instance, owner):
            data = Gauss_log_to_abslines(instance)
            return data

    abs_lines = AbsLines_getter()


def Gauss_log_to_abslines(log_file):
    file = open(log_file.path, "r")
    shift = log_file._shift

    list_file = list()
    fs = list()
    wav = list()

    parameters = list()
    # electrons = list()
    orbitals = list()
    all_orbitals = list()

    reg_exc = re.compile("Excited State")
    reg_wav = re.compile("\s[\d.]*\snm")
    reg_f = re.compile("f\=[\d.]*\s")
    reg_par = re.compile("^(\s?#)p\s")
    reg_line = re.compile("^(\s?-)-*")
    reg_charge = re.compile("^\s?Charge.*\sMultiplicity")
    reg_MOs = re.compile("\d+\s[(alpha)(beta)]+\selectrons")
    reg_num = re.compile("[\d.]+")

    ki = False
    ori = False
    HOMO = 0

    multiplicity = {0: "singlet",
                    1: "doublet",
                    2: "triplet",
                    3: "quartet",
                    4: "quintuplet"}

    for x in file:
        if reg_MOs.search(x) is not None and HOMO == 0:
            electrons = reg_num.findall(x)
            diff = abs(int(electrons[0]) - int(electrons[1]))

            if diff in multiplicity:
                pass
                #print("INFO:  Spin multiplicity is " + str(diff + 1) + ". It's a " + multiplicity[diff] + ".")
            else:
                pass
            HOMO = max(electrons)

        if ori is True:
            if len(reg_num.findall(x)) == 3:
                numerals = reg_num.findall(x)
                orb_nos = numerals[:2]
                amplitude = float(numerals[2])
                ltt = list(map(lambda z: homo_lumo(HOMO, z), orb_nos))
                ltt.append(amplitude)
                orbitals.append(ltt)
            else:
                all_orbitals.append(orbitals)
                orbitals = []
                ori = False

        if reg_exc.search(x) is not None and ori is False:
            list_file.append(x)
            ori = True

        if ki is True:
            if reg_line.search(x) is None:
                parameters.append(x)
            else:
                ki = False

        if reg_par.search(x) is not None:
            parameters.append(x)
            ki = True

        if reg_charge.search(x) is not None:
            ch_mult = x[1:]

        if reg_MOs.search(x) is not None:
            reg_num.findall(x)
    file.close()

    header = ""
    for q in parameters:
        header = header + q[1:len(q) - 1]

    for y in list_file:
        loc_wav = reg_wav.search(y)
        loc_f = reg_f.search(y)
        wavelength = Energy(float(y[loc_wav.start() + 1: loc_wav.end() - 3]), 'nm')
        # energy shift needs to be applied in linear energy units, not wavelength units
        shifted_wavelength = wavelength + shift
        wav.append(shifted_wavelength.in_units("nm"))
        fs.append(float(y[loc_f.start() + 2: loc_f.end() - 1]))

    if wav is []:
        raise Exception("Error in file import. Check the encoding of .txt file and eventually change it to ANSI.")

    if shift.value != 0:
        print(f"Shift of {shift.in_current_units()} {get_current_energy_units()} was applied to the energy axis:\n")

    output = list(map(lambda i, j, k: AbsLine(i, j, k), wav, fs, all_orbitals))

    return output
