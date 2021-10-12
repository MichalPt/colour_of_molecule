import re
from colour_of_molecule.classes import AbsLine, File, Energy, OpList
from colour_of_molecule.analysis.common_tools import homo_lumo, remove_dashes, wrap_from_both_sides


class Configuration:
    def __init__(self, indices, config_a, config_b, amplitude):
        self.indices = indices
        self.alpha = OpList(config_a)
        self.beta = OpList(config_b)
        self.amplitude = amplitude
        self.length = self.get_length()

    def get_length(self):
        la = len(self.alpha)
        lb = len(self.beta)
        if la != lb:
            raise Exception("ERROR:\tConfigurations for alpha an beta electrons differ in length for a given state.")
        return la


class State:
    def __init__(self, number, multiplicity, energy):
        self.number = number
        self.multiplicity = multiplicity
        self.energy = Energy(energy, "eV")
        #self.conf = list()
        #
        self.configurations = list()

        self.MO_indices = list()

    def add_config_line(self, MO_indices, config_a, config_b, amplitude):
        self.configurations.append(Configuration(MO_indices, config_a, config_b, amplitude))
        self.MO_indices = MO_indices


def create_abs_lines(state1, state2, fosc, wavelength, homo):
    initial = state1.configurations[0]
    finals = state2.configurations

    def determine_change_in_occ(initial, final):
        change_in_occ = list(map(lambda i, j: len(i.replace("-", "")) - len(j.replace("-", "")), final, initial))
        return change_in_occ

    def determine_MOs_of_transition(change_in_occ, MO_ids, homo_id):
        trans_from = [homo_lumo(homo_id, y) for x, y in zip(change_in_occ, MO_ids) if x < 0]
        trans_to = [homo_lumo(homo_id, y) for x, y in zip(change_in_occ, MO_ids) if x > 0]
        return trans_from, trans_to

    def ids_to_MO_names(initial, final, MO_ids, homo_id):
        dif_a = determine_change_in_occ(initial.alpha, final.alpha)
        dif_b = determine_change_in_occ(initial.beta, final.beta)
        output = list()

        for x in dif_a, dif_b:
            if sum([abs(y) for y in x]) != 0:
                output.append(list(determine_MOs_of_transition(x, MO_ids, homo_id)))
                output[0].append(final.amplitude)
        return output

    transitions = list()
    for cf in finals:
        data = ids_to_MO_names(initial, cf, state1.MO_indices, homo)
        if len(data) != 0:
            transitions.append(data[0])

    output = AbsLine(wavelength, fosc, transitions)
    return output


class File_mndo(File):
    class AbsLines_getter:
        def __get__(self, instance, owner):
            data = Mndo_file_to_abslines(instance.path)
            return data

    abs_lines = AbsLines_getter()


def Mndo_file_to_abslines(log_file, wav_shift=0.0):
    with open(log_file, "r") as file:

        states = dict()
        config = list()
        output = list()

        reg_state = re.compile("^\s?State\s*\d+,")
        reg_line = re.compile("^\s?-+")
        reg_num = re.compile("[\d.]+")
        reg_ab = re.compile("ab|[ab-]{1}")

        def filter_ab(lst, ab):
            out = [x.replace("-", "").replace(ab, "") for x in lst]
            return out

        read = False
        transitions = False

        state_no = 0
        empty_line = 0
        empty_line2 = 0
        init_state = int()

        for line in file:
            ### first part - first tables with states descriptions
            if reg_state.search(line) is not None:
                nums = reg_num.findall(line)
                state_no = int(nums[0])
                states.update({state_no: State(state_no, nums[2], nums[4])})
                empty_line = 0

            if read is True:
                if reg_num.search(line) is not None:
                    numerals = reg_num.findall(line)
                    amplitude = float(numerals[1])
                    indices = numerals[3:]

                elif reg_ab.search(line) is not None:
                    config = reg_ab.findall(line)
                    config_a = filter_ab(config, "a")
                    config_b = filter_ab(config, "b")
                    states.get(state_no).add_config_line(indices, config_a, config_b, amplitude)

                    config = []
                    config_a = []
                    config_b = []
                    indices = []
                    amplitude = 0
                    empty_line = 0

                else:
                    empty_line += 1

            if empty_line == 2:
                read = False

            if state_no != 0 and reg_line.search(line) is not None and empty_line != 2:
                read = True


            ### second part - second table
            if empty_line2 >2:
                transitions = False

            if "Properties of transitions" in line:
                transitions = True
                init_state = int(reg_num.findall(line)[0])

            elif transitions is True:
                if "Symmetry" not in line and reg_num.search(line) is not None:
                    numerics = [float(x) for x in reg_num.findall(line)]
                    homo_id = max([i for i, j in zip(states.get(1).MO_indices, states.get(1).configurations[0].alpha) if len(j) > 0])
                    abline = create_abs_lines(states.get(init_state), states.get(int(numerics[0])), numerics[5], numerics[4], homo_id)
                    output.append(abline)

                elif "Symmetry" not in line:
                    empty_line2 += 1

    if wav_shift != 0:
        print("Shift applied to wavelengths:\n  ", wav_shift, " nm")

    return output

