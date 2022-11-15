# Polhemus "DIGIGUI"
MATLAB GUI for the Polhemus PATRIOT Digitiser, most recently tested with MATLAB R2016b.

The GUI can be launched by running `DIGIGUI.m` using MATLAB. Ensure that all other `.m`, `.mat` and `.fig` files in this folder are included in the matlab directory.
Once the GUI is running, hover the curser over buttons in the user interface for more information about their function.

A list of points to digitise is imported from a text file where each point is on a new line. The `headpoints test.txt` and `headpoints test long.txt` files serve as a reference and can be imported as examples.

Points are digitised by pressing the stylus button.

For more information see the comments in `DIGIGUI.m`.

This has been designed for use with MATLAB compiler using the `DIGIGUI.prj` file.

Copyright (C) 2022 Reuben W. Nixon-Hill (formerly Reuben W. Hill)

This package is licenced under the MIT Licence (see LICENCE.txt).