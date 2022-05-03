#!/usr/bin/env python3

# Merge all the source files into one single document for the notebook.

files = ['src/preamble.py',
         'src/solvents.py',
         'src/solutes.py',
         'src/gas.py',
         'src/droplet.py',
         'src/gui.py',
         'src/equations.py',
	 'src/benchmarking.py']
source = []
for path in files:
    with open(path) as f:
        source += f.readlines()
source = ''.join(source)

# Convert the raw text into a jupyter notebook.

import jupytext
notebook = jupytext.reads(source, fmt='py')

# Remove any cells we have been instructed to ignore in the source files.

notebook.cells = list(filter(lambda cell: 'ignore' not in cell.metadata or
                             cell.metadata['ignore'] == False, notebook.cells))

# Remove "if __name__ == '__main__':" boilerplate from any cells that have it:

for cell in notebook.cells:
    metadata = cell.metadata
    code = cell.source

    split_by_lines = code.split('\n')
    first_line = split_by_lines[0]

    nspaces = 4
    if "if __name__" in first_line and "'__main__':" in first_line:
        # Remove first line with the boilerplate and any indentation that follows.
        split_by_lines = split_by_lines[1:]
        remove_indent = lambda line: line[nspaces:] if len(line) > nspaces else line
        split_by_lines = [remove_indent(line) for line in split_by_lines]
        cell.source = '\n'.join(split_by_lines)

# Save the notebook.

jupytext.write(notebook, 'sadkat.ipynb')
