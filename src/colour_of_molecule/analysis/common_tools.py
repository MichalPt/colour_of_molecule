
def homo_lumo(HOMO_no, orbital_no):
    if int(orbital_no) <= int(HOMO_no):
        output = "HOMO"
        d = int(orbital_no) - int(HOMO_no)
    else:
        output = "LUMO"
        d = int(orbital_no) - int(HOMO_no) - 1
    if d > 0:
        output += "+"+str(d)
    elif d < 0:
        output += str(d)
    return output


def remove_dashes(lst):
    return list(map(lambda x: x.replace("-",""), lst))


def determine_number_of_calculations(file):
    type = file.type
    switch = {
        "gaussian": "",
        "orca": ""
    }
    pass


def normalize_list(lst):
    norm = max(lst)
    return [x/norm for x in lst]


def wrap_from_both_sides(middle, symbol="*"):
    return str(symbol) + str(middle) + str(symbol)[::-1]


def file_saver(path):
    import os
    number = 0

    while os.path.isfile(path) is True:
        pth, ext = path.rsplit('.', 1)

        if number == 0:
            path = "_{:d}.".format(number).join([pth, ext])
        else:
            path = "_{:d}.".format(number).join([pth[:-(len(str(number - 1)) + 1)], ext])
        number += 1

    return path


def get_output_path(input_path, name="IMG_col-of-mol.png"):
    import os
    dir = os.path.dirname(input_path)
    output = os.path.join(dir, name)
    return output



# def check_save():
#     import sys, os
#     from shutil import copy
#
#     if sys.argv[1] == 'save':
#         com_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#         templates_dir = os.path.join(com_dir, 'templates')
#         script_name = os.path.basename(sys.argv[0])
#         script_path = os.path.join(os.getcwd(), script_name)
#         target_path = os.path.join(templates_dir, script_name)
#
#         print(templates_dir)
#         print(script_path)
#
#         rename = input("\nDo you want to archive the script \"{}\" in package templates? [Y/n] \n".format(script_name))
#         if rename == "Y":
#             copy_path = file_saver(target_path)
#         else:
#             print("\nINFO:\tScript archiving canceled.")
#             sys.exit(0)
#
#         print("\nINFO:\tFile \"{}\" will be copied \n\t> from \"{}\" \n\t> to \"{}\"".format(script_name, script_path, copy_path))
#
#         confirmation = input("\nPress Enter to proceed. \n")
#
#         if confirmation == "":
#             copy(script_path, copy_path)
#             print("INFO:\tScript successfully copied to \"{}\"\n".format(copy_path))
#             sys.exit(0)
#
#         else:
#             print("INFO:\tScript archiving canceled.")
#             sys.exit(0)


