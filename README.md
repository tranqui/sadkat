# sadkat: single-aerosol drying kinetics and trajectories
**A python package to model aerosol drying, brought to you by this cat's tears.**

<img alt="The droplets are brought to you by a very sad kitty" src="https://thypix.com/wp-content/uploads/sad-cat-2.jpg" width="300" />

## About

This is a software package designed to model the kinetics of aerosols under various conditions.

This is not intended as a black box: we tried to keep the python code as simple as possible so the user can see how the model is implemented. No knowledge of python is needed to actually use the model through the GUI, however with a little python experience it _should_ be straightforward to extend the model (by e.g. introducing new solutions) or automate running simulations over a range of conditions.

This package implements a widely used theoretical model for aerosol drying described in the following literature: 
* Xie, X. _et al_, "How far droplets can move in indoor environments - revisiting the Wells evaporation-falling curve", Indoor Air 17(3), 211-225 (2007). DOI: [10.1111/j.1600-0668.2007.00469.x](https://doi.org/10.1111/j.1600-0668.2007.00469.x).
* Liu, L. _et al_, "Evaporation and dispersion of respiratory droplets from coughing", Indoor Air 27(1), 179-180 (2017). DOI: [10.1111/ina.12297](https://doi.org/10.1111/ina.12297).
* Walker, J. S. _et al_, "Accurate Representations of the Microphysical Processes Occurring
during the Transport of Exhaled Aerosols and Droplets", ACS Cent. Sci.(7), 200-209 (2021). DOI: [10.1021/acscentsci.0c01522](https://doi.org/10.1021/acscentsci.0c01522).

Authors:
* Joshua Robinson (joshua.robinson@bristol.ac.uk)
* Dan Hardy (dan.hardy@bristol.ac.uk)

Elements of the code were inspired by an earlier numerical model built by Jim Walker (jim.walker@bristol.ac.uk).

## Installation instructions

### Windows

TBD.

### Linux/Mac

Make sure [jupytext](https://jupytext.readthedocs.io/en/latest/install.html) is installed, then run:

    python make.py

inside the source directory. This should generate a notebook `sadkat.ipynb` which can be opened from your browser after launching the jupyter notebook kernel. That is, run:

    jupyter notebook

which should open a new window in your browser. You should then be able to select the new notebook.


## How to use

TBD.
