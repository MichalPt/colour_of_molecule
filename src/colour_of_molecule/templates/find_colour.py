### category: Colours_and_plotting

import sys
import colour_of_molecule as com

# Path to the file can be entered here:
input_file = "C:/"
save_path = ""


# Accepting command line arguments:
if len(sys.argv) > 1:
    path = sys.argv[1]
else:
    path = input_file


# Loading the file:
file = com.file_in(path)


# Font settings
from colour_of_molecule.classes.classes import FontSettings
fonts = FontSettings(newfonts={'all': 'Arial'}, newsizes={'title': 11}, use_all=True)


# Generate output filenames
if save_path == "":
    from colour_of_molecule.analysis.common_tools import get_output_path, file_saver
    save_path = get_output_path(path, name="{0}_IMG_colmol.png".format(file.filename.split(".")[0]))


# Generating output:
com.get_colour(file, save=file_saver(save_path), fonts=fonts)


