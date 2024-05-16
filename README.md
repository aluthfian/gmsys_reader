# gmsys_reader
A legacy GM-SYS 2D model file reader coded in Python

<img src="https://github.com/aluthfian/gmsys_reader/blob/main/example_rangitoto_luthfian23.png" height="500">

# Aim
This Python file is used to read GM-SYS legacy model files:
1. The `sur` file contains model bodies' geometries.
2. The `blk` file contains model bodies' physical properties (density, magnetic susceptibility, and remanent magnetisation).
3. The `ecs` file contains model origin, azimuth, and coordinate system information.
4. The `wel` file contains wells used in the model.
5. The `grv` file contains the gravity data used to make the model.
6. The `mag` file contains magnetic observation data.

# Warning
⚠️ The `ecs` file reader is optimised for New Zealand use (it reads [NZTM](https://epsg.io/2193) and [NZMG](https://epsg.io/27200) coordinate systems). You may want to tweak the `gmsys_reader.py` to suit your local coordinate system.

# How to Use
Download the `gmsys_reader.py` into a folder containing the legacy GM-SYS files, and import it onto your project by typing `from gmsys_reader import *`. BLK and SUR files will be returned as a Pandas DataFrame, while the other files will be read as a mix of NumPy array and Pandas DataFrame. Please see `sample_rangitoto_plot.ipynb` for example.

# About the Model Used Here
The model used in the example comes from [Fig. 5C](https://www.sciencedirect.com/science/article/pii/S0377027323000811#f0025) of my publication [Luthfian et al. (2023)](https://www.sciencedirect.com/science/article/pii/S0377027323000811). The article and figures are licensed under [CC BY 4.0 DEED](https://creativecommons.org/licenses/by/4.0/) Creative Commons license. You are free to share, copy, redistribute, adapt, and build upon the material for any commercial or non-commercial purpose, with **attribution** and **no additional restrictions** applicable.

# Reference
* [GM-SYS User's Guide Version 4.9](https://docplayer.net/19843267-Gm-sys-user-s-guide-version-4-9.html)
