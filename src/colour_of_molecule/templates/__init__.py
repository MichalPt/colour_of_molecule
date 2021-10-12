from os import listdir
import re, os, sys
from shutil import copy
from colour_of_molecule.analysis.common_tools import wrap_from_both_sides


def update_dict(dictionary, key, value):
    if dictionary.get(key):
        dictionary[key].append(value)
    else:
        dictionary[key] = [value]


def interactive_cmd(opt1, *args):
    from InquirerPy import prompt, inquirer

    def add_points(lst):
        for j in lst:
            key = list(j.keys())[0]
            val = "> " + list(j.values())[0]
            j.update({key: val})
        return lst

    opt1 = add_points(opt1)
    back = True

    while back is True:
        exit_option = [{"name": "-- exit", "value": "exit"}]
        question1 = [
            {
                "type": "list",
                "message": "What category of templates are you interested in?  (use arrows to navigate)",
                "choices": opt1 + exit_option,
            }, ]
        result1 = prompt(questions=question1)

        if result1[0] == "exit":
            raise IOError

        options2 = args[result1[0]]

        back_option = [{"name": "-- back", "value": True}]
        question2 = [
            {
                "type": "list",
                "message": "Choose template you wish to create.  (use arrows to navigate)",
                "choices": options2 + back_option
            }, ]

        result2 = prompt(questions=question2)

        if result2[0] is not True:
            back = False

    return result1[0], result2[0]


def static_cmd(opt1, *args):
    def print_static(level_1, *args):
        print("Select a template you wish to import:\n")
        counter = 0
        output = dict()
        for i in level_1:
            label = str(i.get("name"))
            index = i.get("value")
            print(label)
            for k in args[index]:
                print("    " + str(counter) + "  >  " + k.get("name"))
                output.update({counter: k})
                counter += 1

        print("\nThen run a function \"colour_of_molecule.templates.create(#)\" where # is the index of selected file to copy it into current directory.")
        return output

    global indexed_templates
    indexed_templates = print_static(opt1, *args)

    global create
    def create(index, path, list=indexed_templates):
        answ = list.get(index).get("name")
        template_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), answ)
        print("path: ", path)
        if "." in path and ".py" in path:
            print("INFO:\tCopying template" + wrap_from_both_sides(answ, " \"") +
                  "\n\t> from" + wrap_from_both_sides(template_path, " \"") +
                  "\n\t> to" + wrap_from_both_sides(path, " \"")
                  )
            copy(template_path, path)
        else:
            target_path = os.path.join(os.path.dirname(path), answ)
            print("INFO:\tCopying template"+wrap_from_both_sides(answ, " \"")+
                  "\n\t> from"+wrap_from_both_sides(template_path, " \"")+
                  "\n\t> to"+wrap_from_both_sides(target_path, " \"")
                  )
            copy(template_path, target_path)

        return answ


### always running:

path = os.path.dirname(os.path.abspath(__file__))

reg_extension = re.compile("\.py$")
files = [f for f in listdir(path) if reg_extension.search(f) != None and f != "__init__.py"]

files_dict = dict()

for filename in files:
    with open(path+"/"+filename, "r") as file:
        for line in file:
            if "category:" in line and "#" in line:
                cat = line.split(":")[1].replace("\n","").replace(" ","")
                update_dict(files_dict, cat, filename)
            break

options1 = [{"name": x, "value": i} for (i, x) in enumerate(files_dict.keys())]
options2 = [[{"name": y, "value": j} for (j, y) in enumerate(files_dict.get(z))] for z in files_dict.keys()]
name_is_ok = False

try:
    answer = interactive_cmd(options1, *options2)
    file_dir_path = os.path.abspath(sys.argv[0])
    selected_file = options2[answer[0]][answer[1]].get("name")
    template_path = os.path.join(path, selected_file)

    files_in_target_dir = [f for f in os.listdir(file_dir_path) if not f.startswith("~$")]

    while selected_file in files_in_target_dir:
        print("WARNING:\tScript called" + wrap_from_both_sides(selected_file, " \"") + "already exists in current directory.")
        overwrite = input("\nDo you want to overwrite it? [Y/n] \n")

        if overwrite == "Y":
            break

        elif overwrite == "n":
            rename = input("\nDo you want to save it with different name? [Y/n] \n")

            if rename == "Y":
                while name_is_ok is False:
                    new_file_name = input("\nEnter new script name: \n") or selected_file.replace(".py", "_new.py")

                    if ".py" in new_file_name:
                        name_is_ok = True
                        selected_file = new_file_name
                    else:
                        print("ERROR:\tInvalid file name. It doesn't seem to have the \".py\" extension. Reinitializing ...\n")
                        continue

            else:
                print("\nINFO:\tScript archiving canceled.")
                sys.exit(0)

        else:
            print("INFO:\tArgument hasn't been recognised. Script archiving canceled.")
            sys.exit(0)

    target_path = os.path.join(file_dir_path, selected_file)
    print("INFO:\tFile" + wrap_from_both_sides(selected_file, " \"") + "will be copied" +
          "\n\t> from"+wrap_from_both_sides(template_path, " \"") +
          "\n\t> to" + wrap_from_both_sides(target_path, " \"")
          )

    confirmation = input("\nPress Enter to proceed. \n")

    if confirmation == "":
        copy(template_path, target_path)
        print("INFO:\tScript successfully copied to" + wrap_from_both_sides(target_path, " \""))
        sys.exit(0)

    else:
        print("INFO:\tScript archiving canceled.")
        sys.exit(0)

except IOError:
    pass

except SystemExit:
    pass

except:
    static_cmd(options1, *options2)



