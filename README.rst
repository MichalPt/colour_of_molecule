"""""""""""""
Documentation
"""""""""""""

.. contents:: Contents:

============
Installation
============

++++++++++++++
Prerequisities
++++++++++++++

To be able to install and use the package ``colour_of_molecule``, Python language (version >= 3.6) needs to be installed on your computer.

You can check the installation by passing ``python`` or ``python3`` to the console. If Python distribution is ready to be used, you will enter the Python console:

.. code-block:: console

 Python 3.9.7 (tags/v3.9.7:1016ef3, Aug 30 2021, 20:19:38) [MSC v.1929 64 bit (AMD64)] on win32
 Type "help", "copyright", "credits" or "license" for more information.
 >>> |

Any time you would like to close the Python console and return to the standard console, call:

.. code-block:: console

 quit()

Besides, you will need to have installed a ``pip`` package manager.
It can be installed by the following command:

.. code-block:: console

 python -m install pip

++++++++++++++++++++
Package installation
++++++++++++++++++++

The package, along with all necessary dependencies, can be installed by the following command:

.. code-block:: console

 pip install colour-of-molecule

The installation can be checked by calling:

.. code-block:: console

 python

and then trying to import the package:

.. code-block:: console

 import colour_of_molecule

You can again leave the Python console by calling ``quit()`` if no error has occurred.

++++++++++++++++++++
Updating the package
++++++++++++++++++++

You can upgrade the package to the current version by:

.. code-block:: console

 pip install colour-of-molecule --upgrade

or alternatively by a shorter command:

.. code-block:: console

 pip install colour-of-molecule -U

If you wish to install a specific version, the command might look like this:

.. code-block:: console

 pip install colour-of-molecule==0.0.2.dev3

++++++++++++++
Uninstallation
++++++++++++++

The package can be completely removed from your machine by following command:

.. code-block:: console

 pip uninstall colour-of-molecule

=====
Usage
=====
++++++++++++++++++++++++++
Importing template scripts
++++++++++++++++++++++++++

The package contains several preset template scripts which can be copied to current folder at any time by following commands.

Initialize Python console:

.. code-block:: console

 python

Then import the templates:

.. code-block:: console

 import colour_of_molecule.templates

An interactive menu should appear:

.. code-block:: console

 >>> import colour_of_molecule.templates
 ? What category of templates are you interested in? (use arrows to navigate)
  > Colours_and_plotting
  > Multiple_files_manipulation
  ... custom folders ...
  --exit

Follow the instructions and select the desired .py script by using arrows and enter keys. You will be asked to confirm the creation of the selected .py script in the directory the console was navigated into. For example if the Python console was invoked in ``C:\Users\Joe`` folder and the script ``plot_spectrum.py`` was selected, the confirmation might look like this:

.. code-block:: console

 INFO:   File "plot_spectrum.py" will be copied
         > from "C:\Users\Joe\miniconda3\envs\env-01\lib\site-packages\colour_of_molecule\templates\plot_spectrum.py"
         > to "C:\Users\Joe\plot_spectrum.py"

 Press Enter to proceed.
 |

The saving process contains failsafe against possible file overwrite. You will be asked to enter a new script filename or to confirm the ovewrite if any filename collision was found.

+++++++++++++++++++++++++++++++++++++++++++++
Alternative way of importing template scripts
+++++++++++++++++++++++++++++++++++++++++++++

If your console doesn't support interactive prompt provided by ``InquirerPy`` Python package (section `Importing template scripts`_), an alternative menu might be displayed:

.. code-block:: console

 >>> import colour_of_molecule.templates
 Select a template you wish to import:
 > Multiple_files_manipulation
     0  >  analyze_multiple_files.py
 > Colours_and_plotting
     1  >  find_colour.py
     2  >  plot_spectrum.py
 Then run a function "colour_of_molecule.templates.create(#)" where # is the index of selected file to copy it into current directory.

 >>> |

Follow the listed instructions and create the desired script by calling, for example (#=1):

.. code-block:: console

 colour_of_molecule.templates.create(1)

++++++++++++++++++++++++++++++++++++++++++++++
Archive a new script inside the package folder
++++++++++++++++++++++++++++++++++++++++++++++

If you want to make your script easily accessible by the template script importing mechanism listed above, you can archive your own custom script inside the package folder along with the template ones. Please **keep in mind that the** ``colour_of_molecule`` **package update might remove or overwrite these archived custom scripts** so please store them somewhere else as well to keep them safe in a longterm perspective.

The scripts can be distributed into "virtual" folders or categories that will be displayed during the script import process to keep it organized. The folder assignment is done in the first line of the script itself by a following text:

.. code-block:: python

 ### category: folder-name

To add the custom scipt to the templates folder, within the package, navigate to the folder your custom script is currently stored. Then use the following command similar to the one normally used to run the script itself but with the ``--save`` keyword added to the command instead of the input file path. For example it might look like followlingly:

.. code-block:: console

 python plot_spectrum2.py --save

You will be asked to confirm the archiving or to enter a new filename if the current is already used in the templates folder.

+++++++++++++++++++
Running the scripts
+++++++++++++++++++

There are multiple ways how the scripts can be used. You can either specify the input file path in the script itself, for example into a variable called ``input_file``:

.. code-block:: python

 import colour_of_molecule as com      # importing the package
 input_file = "C:/..."                 # specifying the input file path
 file = com.file_in(input_file)        # loading the input file

The script would be then run by a simple command:

.. code-block:: console

 python my_script.py

Alternatively, you can also pass the input file path straight from the command line, for example:

.. code-block:: console

 python my_script.py ./gaussian/asp-B3LYP-pVDZ.log

It's important to mention that this way of passing the input file path as an in-line argument is possible if and only if the script contains appropriate piece of code which enables it:

.. code-block:: python

 import sys                 # importing python built-in package
 if len(sys.argv) > 1:      # checks the number of arguments passed to python (0 - script, 1 - input file path)
     path = sys.argv[1]     # sets the variable path to the input file path (argument with index 1)
 else:
     path = input_file      # if only one argument was passed (i.e. only the script), use the in-file specified path (see the beginning of this docs section)


===========================
Code structure and commands
===========================

All settings related to numerical parameters or analysis enters the process via the class ``File``. Setting related to fonts are managed by class ``FontSettings``.

++++++++++
class File
++++++++++

The first step every script has to contain is the command to load the input file. This is done by ``file_in()`` function directly accessible directly from the package directly. It takes a single argument - path to the input file. For example:

.. code-block:: python

 import colour_of_molecule as com
 file = com.file_in(PATH)

Currently, output formats of four QCh programs are supported: **Gaussian**, **ORCA**, **MNDO**, and **MOLPRO**. The format will be recognised automatically during the loading process.

Any settings are now passed to the ``file`` object (an instance of ``File`` class) as attributes: ``file.X`` where ``X`` can be:

o ``.wavelength_range``
 wavelength range to be plotted

 e.g.: ``file.wavelength_range = (250,850)``

o ``.standard_deviation``
 sets the width of gaussian peaks used to create absorption spectrum

 e.g.: ``file.standard_deviation = 3096.01`` (default value)

o ``.optical_density``
 sets the optical density used to calculate the complementary absorption spectrum needed to determine the actual colour

 e.g.: ``file.optical_density = 0.15`` (default value)

o ``.transition_minimal_amplitude``
 sets the minimal transition amplitude which will be included in the plot of absorption lines

 e.g.: ``file.transition_minimal_amplitude = 0.5`` (default value)

o ``.normalize_absorption_spectrum``
 determine if the absorption spectrum should be normalized to 1 at maximum value

 e.g.: ``file.normalize_absorption_spectrum = False`` (default value)

o ``.normalize_complementary_spectrum``
 determine if the complementary absorption spectrum should be normalized

 e.g.: ``file.normalize_complementary_spectrum = True`` (default value)

Setting related to plotting:

o ``.plot_title``
 sets custom title to the plots, string needs to be enquoted

 e.g.: ``file.plot_title = ""`` (default value)

o ``.legend_title``
 sets custom title to the legend, string needs to be enquoted

 e.g.: ``file.legend_title = ""`` (default value)

++++++++++++++++++
class FontSettings
++++++++++++++++++

All settings related to fonts used and displayed in the plots are managed by the ``FontSettings`` class. To begin with, the class needs to be imported:

.. code-block:: python

 from colour_of_molecule.classes.classes import FontSettings

After that, the class can be instatiated while taking up to two keyword arguments: ``newfonts``, ``newsizes``; and a single boolean keyword argument ``use_all``.
Both keyword arguments has to be dictionaries and the can specify font or font size for these keys:

o ``all``
 it is used for all text if ``use_all = True``

o ``title``
 title of the plot

o ``axis``
 x and y axis labels

o ``axis_tick_labels``
 x and y axis tick labels (i.e. numbers adjacent to axis ticks)

o ``legend``
 title of the legend and the whole legend itself

The default font is *Calibri* and the default font size is *14* for plot title and *12* for everything else.

The final usage might look like this:

.. code-block:: python

 font_settings = FontSettings(newfonts={'all': 'Consolas'}, newsizes={'title': 11, 'legend': 8}, use_all=True)

The instance can be then passed to any of the plotting functions, for example:

.. code-block:: python

 com.plot_single_spectrum(file, fonts=font_settings)

++++++++++++++++++
Plotting functions
++++++++++++++++++

There are currently three functions capable of returning an image of a plot:

o ``plot_single_spectrum()``

o ``plot_abs_lines()``

o ``get_colour()``

Each of these functions takes a single positional argument - an instance of class ``File`` - and various keyword arguments.

The keyword arguments can be categorised into two groups - **general** and **function-specific**.

--------------------------
General keyword arguments
--------------------------

o ``save``
 sets the path where to save the output image

 e.g.: ``com.plot_single_spectrum(file, save="C:/...")``

o ``title``
 title of the plot displayed in its header

 e.g.: ``com.plot_single_spectrum(file, title="Example1")``
 
 to hide the title use expression ``title=None`` (please note that in case of axis labels the preferred way is to use null string ``""`` instead)

o ``xaxis_label``
 label for the x-axis displayed below the plot

 e.g.: ``com.plot_single_spectrum(file, xaxis_label="wavelength [nm]")`` (default value)

o ``yaxis_label``
 label for the y-axis displayed on the left side of the plot
 
 e.g.: ``com.plot_single_spectrum(file, xaxis_label="wavelength [nm]")`` (default value)
 
o ``yaxis_label``
 label for the y-axis displayed on the left side of the plot
 
 e.g.: ``com.plot_single_spectrum(file, yaxis_label="relative absorbance")``

o ``yaxis_label_right``
 label for the right y-axis displayed on the right side of the plot

 e.g.: ``com.plot_single_spectrum(file, yaxis_label_right="oscillator strength")``

o ``size``
 tuple, diameters of the plot expressed by a tuple of values, i.e. ``(width, height)``

 e.g.: ``com.plot_single_spectrum(file, size=(6,4), )``

o ``dpi``
 resolution of the generated image (dots per inch)

 e.g.: ``com.plot_single_spectrum(file, dpi=400)`` (default value)

o ``fonts``
 ... already mentioned above

-----------------------------------
Function-specific keyword arguments
-----------------------------------

**com.plot_single_spectrum()**:

o ``lines_show``
 boolean, True if absorption lines should be plotted below the spectrum, False if not

 e.g.: ``com.plot_single_spectrum(file, lines_show=True)`` (default value)

o ``lines_ratio``
 tuple, sets the relative height of the main plot area (where spectrum is plotted) and the supportive stripe with positions of abs. lines

 e.g.: ``com.plot_single_spectrum(file, lines_ratio=(14,1), )`` (default value)

o ``lines_colours``
 boolean, True if absorption lines with oscillator strength larger than ``lines_lim`` should be coloured according to their corresponding wavelength, otherwise they will be coloured black

 e.g.: ``com.plot_single_spectrum(file, lines_colours=True)`` (default value)

o ``lines_lim``
 float, limiting value of oscillator strength separating so-called dark and bright transitions

 e.g.: ``com.plot_single_spectrum(file, lines_lim=0.0001)`` (default value)

o ``lines_width``
 float, width of plotted abs. lines in pts.

 e.g.: ``com.plot_single_spectrum(file, lines_width=1.2)`` (default value)

o ``rainbow``
 boolean, True if a colour spectrum should be displayed below the line of plotted abs. spectrum line

 e.g.: ``com.plot_single_spectrum(file, rainbow=True)`` (default value)

===============
Example outputs
===============

+++++++++
Example 1
+++++++++

.. code-block:: python

 file = com.file_in("C:/Users/xyz/carotenoid.out")
 file.wavelength_range = (100, 1000)
 com.plot_single_spectrum(file, save="C:/Users/xyz/exp1.png", dpi=200, size=(10, 3),
                         title="Carotenoid", xaxis_label="wavelength [nm]", yaxis_label="relative absorbance",
                         lines_show=False,
                        )

.. image:: https://github.com/MichalPt/colour_of_molecule/blob/6855ea3d8a149b7eb3b4c72048ecf5a42d50af85/exp1_0.png


+++++++++
Example 2
+++++++++

.. code-block:: python
 
 file = com.file_in("C:/Users/xyz/phenolphtalein.log")
 file.wavelength_range = (200, 700)
 com.plot_single_spectrum(file, save="C:/Users/xyz/exp2.png", dpi=200, size=(10, 3),
                         title=None, xaxis_label="wavelength [nanometers]", yaxis_label="rel. abs.", 
                         lines_show=True, lines_colours=True, lines_lim=0.001, lines_ratio=(12,2), lines_width=1.8,
                        )

.. image:: https://github.com/MichalPt/colour_of_molecule/blob/6855ea3d8a149b7eb3b4c72048ecf5a42d50af85/exp1.png

