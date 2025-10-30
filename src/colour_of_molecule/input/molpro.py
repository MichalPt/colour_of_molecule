import re
from colour_of_molecule.classes import AbsLine, File, Energy
from colour_of_molecule.analysis.common_tools import homo_lumo, wrap_from_both_sides


class DipoleTransMoment:
    def __init__(self):
        self.dmx = 0
        self.dmy = 0
        self.dmz = 0
        self.energy_difference = float()
        self.oscillator_strength = float()

    def get_oscillator_strength(self, dE):
        if isinstance(dE, float):
            self.energy_difference = dE
        else:
            raise Exception("ERROR:\tArgument has to be either a float (in atomic units) or a tuple of two \"Energy\" class instances.")

        self.oscillator_strength = 2/3 * self.energy_difference * (self.dmx**2 + self.dmy**2 + self.dmz**2)

        return self.oscillator_strength


class SingleState:
    def __init__(self, number, energy):
        self.number = float(number)
        self.energy = Energy(energy, "au")


class StateTransition:
    def __init__(self, num_final, energy, num_initial=1.1):
        self.num_final = float(num_final)
        self.num_initial = float(num_initial)
        self.symmetry = int(str(num_final).split(".")[1])
        self.energy = Energy(energy, "au")
        self.wavelength = self.energy.in_units("nm")
        self.osc_strength = float()
        self.MO_indices = list()

        self.dipole_trans_moment = DipoleTransMoment()

    def add_MO_transition(self, mo1, mo2, coeff):
        self.MO_indices.append([float(mo1), float(mo2), float(coeff)])

    def add_osc_strength(self, fosc):
        self.osc_strength = float(fosc)


class MainReader:
    reg_num = re.compile("-?\d+[.\d]*")
    reg_alphabet = re.compile("[a-zA-Z]+")
    reg_res_for_states = re.compile("^\s*Results for state\s*\d+\.\d+")
    reg_oscillator = re.compile("^\s*Oscillator strength")
    reg_orbitals = re.compile("^\sELECTRON ORBITALS")


class EomReader(MainReader):
    reg_eom_final = re.compile("^\s*Final Results for EOM")

    def __call__(self, log_file):
        shift = log_file._shift

        with open(log_file.path, "r") as file:

            state_transitions = dict()
            index = float()
            id_homo = float()
            #id_lumo = float()
            listlong = list()
            longest = int()
            c_empty_line = 0
            c_results_for_states = False
            eom_final = False
            c_oscill = False
            c_orbitals = False

            for line in file:
                # first part - states and transitions
                if self.reg_res_for_states.search(line) is not None:
                    index, eng = [float(x) for x in self.reg_num.findall(line)[:2]]
                    state_transitions.update({
                        index: StateTransition(index, eng)
                    })
                    c_results_for_states = True
                    c_empty_line = 0

                elif c_results_for_states is True and index > 0:
                    if self.reg_num.search(line) is not None and "->" in line:
                        mos = self.reg_num.findall(line)
                        state_transitions.get(index).add_MO_transition(mos[1], mos[2], mos[0])
                        c_empty_line = 0

                    else:
                        c_empty_line += 1

                if c_empty_line > 3:
                    c_results_for_states = False

                # second part - transition dipole moments
                if self.reg_eom_final.search(line) is not None:
                    eom_final = True

                if eom_final is True:
                    if self.reg_oscillator.search(line) is not None:
                        c_oscill = True

                    elif c_oscill is True and self.reg_num.search(line) is not None:
                        nums = [float(x) for x in self.reg_num.findall(line)]
                        state_transitions.get(nums[0]).add_osc_strength(nums[2])
                        c_oscill = False

                # third part - MO indexes, names
                if self.reg_orbitals.search(line) is not None:
                    c_orbitals = True

                if c_orbitals is True and self.reg_alphabet.search(line) is None:
                    length = len(self.reg_num.findall(line))
                    if length > longest and length > 0:
                        longest = length
                        listlong = []
                        listlong.append(self.reg_num.findall(line))
                    elif length == longest:
                        listlong.append(self.reg_num.findall(line))

                if " HOMO " in line:
                    c_orbitals = False
                    id_homo = float(self.reg_num.findall(line)[0])
                #
                # if " LUMO " in line:
                #     c_orbitals = False
                #     id_lumo = float(self.reg_num.findall(line)[0])

        def format_listlong(listlong):
                if not isinstance(listlong, list):
                    raise Exception("ERROR:\tArgument has to be a list.")

                sortlist = list()

                for x in listlong:
                    sortlist.append([float(x[2]), float(x[0]), int(x[1])])

                out = list(enumerate(sorted(sortlist)))

                return {y[1]: x for (x, y) in out}

        def homo_lumo_for_MO(dct, homo, mos):
            if not isinstance(dct, dict):
                raise Exception("ERROR:\tFirst argument has to be a dictionary of floats in both key and value.")
            if not isinstance(id_homo, float):
                raise Exception("ERROR:\tSecond argument has to be a float (orbital.symmetry).")
            if not isinstance(mos, list) and len(mos) != 2:
                raise Exception("ERROR:\tThird argument has to be a list of two floats (initial and final orbital.symmetry).")

            def get_and_check(dct, homo_id, id):
                homo = dct.get(homo_id)

                if dct.get(id) is None:
                    return wrap_from_both_sides(id, symbol="**")
                    #return "*"+str(id)+"*"
                else:
                    return homo_lumo(homo, dct.get(id))

            out = [
                get_and_check(dct, homo, mos[0]),
                get_and_check(dct, homo, mos[1]),
                abs(mos[2])
                   ]
            return out

        dictlong = format_listlong(listlong)
        output = list()

        for transition in state_transitions.values():
            energy = Energy(transition.wavelength, "nm") + shift
            output.append(AbsLine(energy.in_units("nm"),
                                  transition.osc_strength,
                                  [homo_lumo_for_MO(dictlong, id_homo, x) for x in transition.MO_indices]
                                  )
                          )
        return output


class MrciReader(MainReader):
    reg_mrci_energy = re.compile("^\s*!MRCI STATE \d+\.\d+ Energy")
    reg_mrci_trans = re.compile("^\s*!MRCI trans")
    reg_mrci_DMX = re.compile("\|DMX\|")
    reg_mrci_DMY = re.compile("\|DMY\|")
    reg_mrci_DMZ = re.compile("\|DMZ\|")

    def __call__(self, log_file):
        shift = log_file._shift
        states = dict()
        state_transitions = dict()

        with open(log_file.path, "r") as file:
            for line in file:
                if self.reg_mrci_energy.search(line) is not None:
                    nums = [float(x) for x in self.reg_num.findall(line)]
                    states.update({nums[0]: SingleState(nums[0], nums[1])})

                if self.reg_mrci_trans.search(line) is not None:
                    nums_in_string = [float(x) for x in self.reg_num.findall(line)[:2]]

                    if len([x for x in list(states.keys()) if x in nums_in_string]) == 2:
                        state1, state2 = sorted([float(x) for x in nums_in_string])[:2]

                        trans_id = (state1, state2)
                        eng1 = states.get(state1).energy
                        eng2 = states.get(state2).energy
                        delta_eng = eng2 - eng1

                        state_transitions.update({
                            trans_id: StateTransition(state2, delta_eng.in_units("au"), state1)
                        })

                        if self.reg_mrci_DMX.search(line) is not None:
                            dmx = self.reg_num.findall(line)[2]
                            state_transitions.get(trans_id).dipole_trans_moment.dmx = float(dmx)

                        if self.reg_mrci_DMY.search(line) is not None:
                            dmy = self.reg_num.findall(line)[2]
                            state_transitions.get(trans_id).dipole_trans_moment.dmy = float(dmy)

                        if self.reg_mrci_DMZ.search(line) is not None:
                            dmz = self.reg_num.findall(line)[2]
                            state_transitions.get(trans_id).dipole_trans_moment.dmz = float(dmz)

        output = list()

        for transition in state_transitions.values():
            energy = transition.energy + shift
            output.append(AbsLine(energy.in_units("nm"),
                                  transition.dipole_trans_moment.get_oscillator_strength(energy.in_units("au")),
                                  [[wrap_from_both_sides(transition.num_initial), wrap_from_both_sides(transition.num_final)]]
                                  )
                          )
        return output


class Rs3Reader(MainReader):
    def __call__(self, log_file_path):
        print("ERROR:   Rs3 type is not implemented yet.")
        pass
    pass


class File_molpro(File):
    class AbsLines_getter:
        def __get__(self, instance, owner):
            typ, reader = instance.molpro_subtype
            the_reader = reader()
            data = the_reader(instance)
            return data

    class MolproSubtype:
        def __get__(self, instance, owner):
            def check_tuple_re_func(tup):
                from re import Pattern

                if not isinstance(tup, tuple) and len(tup) != 2:
                    raise Exception("ERROR:\tDictionary has to have format \"a: (re.Pattern, function)\".")

                reg, func = tup

                if not isinstance(reg, Pattern):
                    raise Exception("ERROR:\tDictionary has to have format \"a: (re.Pattern, function)\".")

                if not callable(func):
                    raise Exception("ERROR:\tDictionary has to have format \"a: (re.Pattern, function)\".")

            def determine_subtype(lin, supported):
                if not isinstance(supported, dict):
                    raise Exception("ERROR:\tWrong second argument type. Has to be \"dict\".")

                out = tuple()
                for typ, tools in supported.items():
                    check_tuple_re_func(tools)

                    reg_type, func = tools

                    if reg_type.search(lin) is not None:
                        out = (typ, func)
                if out == ():
                    pass
                else:
                    return out

            with open(instance.path, 'r') as file:
                for line in file:
                    found = determine_subtype(line, instance.supported_subtypes)
                    if found is not None:
                        typ, func = found
                        break
            return typ, func

    abs_lines = AbsLines_getter()
    molpro_subtype = MolproSubtype()

    supported_subtypes = {
        "mrci": (re.compile("^\s?\{?ci"), MrciReader),
        "eom": (re.compile("^\s?\{?eom"), EomReader),
        "rs3": (re.compile("^\s?\{?rs3"), Rs3Reader)
    }



