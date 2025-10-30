from colour_of_molecule.classes import File
from colour_of_molecule.utils.energy_units import energy_units, get_current_energy_units
from colour_of_molecule.analysis.common_tools import wrap_from_both_sides
from colour_of_molecule.plotting.plot_spectra import *
import sys, os


def file_in(path):
    return File(path)


if len(sys.argv) > 1:
    script_name = os.path.basename(sys.argv[0])
    script_file_path = os.path.abspath(script_name)
    args = sys.argv[1:]

    if "--save" in args or "-s" in args:
        from shutil import copy

        path0 = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(path0, "templates")

        all_templates = [f for f in os.listdir(path) if ".py" in f and f != "__init__.py"]

        def check_template_format(path):
            with open(path, 'r') as script:
                for line in script:
                    if "#" in line and "category:" in line and len(line.split(":")[-1].replace(" ","")) > 1:
                        return True
                    break

        if check_template_format(script_file_path) is not True:
            raise Exception("ERROR:\tScript couldn't be archived. The first line of the script has to contain \"### category: type_your_category_name_here\"")

        if ".py" not in script_name:
            raise Exception("ERROR:\tScript couldn't be archived. It doesn't seem to have the \".py\" extension.")

        while script_name in all_templates:

            print("WARNING:\tScript called"+wrap_from_both_sides(script_name, " \"")+"is already present in templates.")
            overwrite = input("\nDo you want to overwrite it? [Y/n] \n")

            if overwrite == "Y":
                break

            elif overwrite == "n":
                rename = input("\nDo you want to save it with different name? [Y/n] \n")

                if rename == "Y":
                    new_script_name = input("\nEnter new script name: \n") or script_name.replace(".py", "_new.py")

                    if ".py" in new_script_name:
                        script_name = new_script_name
                    else:
                        print("ERROR:\tInvalid file name. It doesn't seem to have the \".py\" extension. Reinitializing ...\n")
                        continue

                else:
                    print("INFO:\tScript archiving canceled.")
                    sys.exit(0)

            else:
                print("INFO:\tArgument hasn't been recognised. Script archiving canceled.")
                sys.exit(0)


        print("INFO:\tFile" + wrap_from_both_sides(script_name, " \"") + "will be copied"+
              "\n\t> from" + wrap_from_both_sides(script_file_path, " \"") +
              "\n\t> to" + wrap_from_both_sides(os.path.join(path, script_name), " \"")
              )
        confirmation = input("\nPress Enter to proceed. \n")

        if confirmation == "":
            target_path = os.path.join(path, script_name)
            copy(script_file_path, target_path)
            print("INFO:\tScript successfully archived in"+wrap_from_both_sides(target_path, " \""))
            sys.exit(0)

        else:
            print("INFO:\tScript archiving canceled.")
            sys.exit(0)



