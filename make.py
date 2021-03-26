#!/usr/bin/env python3

# Merge all the source files into one single document for the notebook.

files = ['src/preamble.py',
         'src/solvents.py',
         'src/solutes.py',
         'src/gas.py',
         'src/droplet.py',
         'src/gui.py',
         'src/equations.py']
source = []
for path in files:
    with open(path) as f:
        source += f.readlines()
source = ''.join(source)

# Convert the raw text into a jupyter notebook.

import jupytext
notebook = jupytext.reads(source, fmt='py')
jupytext.write(notebook, 'sadkat.ipynb')
