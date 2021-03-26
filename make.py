#!/usr/bin/env python3

import jupytext
notebook = jupytext.read('src/sadkat.py')
jupytext.write(notebook, 'sadkat.ipynb')
