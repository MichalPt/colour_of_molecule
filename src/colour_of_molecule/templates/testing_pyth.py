### category: My_first_category

import os
import importlib.util
import sys

MODULE_PATH = "C:/Users/Michal Ptáček/PycharmProjects/colour_of_molecule_v2/colour_of_molecule/__init__.py"
MODULE_NAME = "colour_of_molecule"

spec = importlib.util.spec_from_file_location(MODULE_NAME, MODULE_PATH)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module 
spec.loader.exec_module(module)

#print(os.path.dirname(os.path.abspath(__file__)))

from colour_of_molecule import templates


