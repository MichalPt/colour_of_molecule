import re
import numpy as np
from colour_of_molecule.classes import AbsLine, File, Energy
from colour_of_molecule.analysis.common_tools import homo_lumo
from colour_of_molecule.utils.energy_units import get_current_energy_units


class File_orca(File):
    class AbsLines_getter:
        def __get__(self, instance, owner):
            data = Orca_file_to_abslines(instance)
            return data

    abs_lines = AbsLines_getter()


def Orca_file_to_abslines(log_file):
    shift = log_file._shift

    with open(log_file.path, "r") as file:

        fs = list()
        wav = list()

        # parameters = list()
        state_no = 0
        states = list()
        state_numbers = list()
        orbitals = list()
        all_orbitals = list()
        num_of_el = list()

        reg_state = re.compile(r"^STATE\s*\d+\:")
        reg_line = re.compile(r"^-+")
        reg_parantheses = re.compile(r"\((.*?)\)")
        # reg_charge = re.compile("^\s?Charge.*\sMultiplicity")

        reg_MOs = re.compile(r"^N(\(Alpha\)|\(Beta\))\s*\:")
        reg_num = re.compile(r"[\d.]+")
        reg_num_float = re.compile(r"[\d\.]{3,}")
        
        edm = False
        edmf = False
        li = 0
        HOMO = int()

        for line in file:
            if reg_MOs.search(line) is not None:
                num_of_el.append(reg_num.findall(line)[0])
                if len(num_of_el) == 2:
                    HOMO = int(round(float(max(num_of_el))))

            if state_no != 0:
                if "->" in line:
                    numerals = reg_num.findall(line)
                    orb_nos = numerals[:2]
                    amplitude = float(numerals[3])
                    ltt = list(map(lambda z: homo_lumo(HOMO, z), orb_nos))
                    ltt.append(amplitude)
                    orbitals.append(ltt)
                else:
                    all_orbitals.append(orbitals)
                    state_numbers.append(state_no)
                    orbitals = []
                    state_no = 0

            if reg_state.search(line) is not None and state_no == 0:
                state_no = int(reg_num.findall(line)[0])
                states.append(state_no)

            if edm is True:
                if li == 1:
                    parantheses = np.array(reg_parantheses.findall(line))
                    if parantheses.size > 2:
                        nm_index = np.argwhere(parantheses == 'nm').item()
                if li == 2:
                    edmf = True
                    edm = False
                elif reg_line.search(line) is not None:
                    li += 1

            if edmf is True:
                if reg_num.search(line) is not None:
                    nums = reg_num_float.findall(line)
                    val = float(nums[nm_index])
                    energy = Energy(val, "nm") + shift
                    wavelength = energy.in_units("nm")
                    wav.append(wavelength)
                    fs.append(float(nums[3]))
                else:
                    edm = False
                    edmf = False

            if "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in line:
                edm = True

    if wav is []:
        raise Exception("Error in file import. Check the encoding of .txt file and eventually change it to ANSI.")

    if shift.value != 0:
        print("Shift applied to wavelengths:\n  {} {}".format(shift.in_current_units(), get_current_energy_units()))

    output = list(map(lambda i, j, k: AbsLine(i, j, k), wav, fs, all_orbitals))

    return output

