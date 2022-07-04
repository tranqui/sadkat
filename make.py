#!/usr/bin/env python3

import jupytext
import sys

def read_source(paths):
    """Parse source file(s) into strings for compilation.

    Args:
        paths: file paths to process.
    Returns:
        String containing contents of source files read from given paths (in order).
    """
    source = []
    for path in paths:
        with open(path) as f:
            source += f.readlines()
    source = ''.join(source)
    return source

class Notebook:
    """Generate a jupyter notebook for the end-user by merging python source files together."""

    def __init__(self, input_paths, output):
        """Generate the notebook.

        Args:
           input_paths (container): locations of the source files in order to combine.
           output: name of generated notebook.
        """
        sys.stderr.write('building %s...' % output)
        sys.stderr.flush()

        self.source = read_source(input_paths)

        # Convert the raw text into a jupyter notebook.
        self.notebook = jupytext.reads(self.source, fmt='py')

        # Remove any cells we have been instructed to ignore in the source files.
        self.notebook.cells = list(filter(lambda cell: 'ignore' not in cell.metadata or
                                          cell.metadata['ignore'] == False, self.notebook.cells))

        self.markdown = []

        # Remove "if __name__ == '__main__':" boilerplate from any cells that have it:
        for cell in self.notebook.cells:
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

            cell.source = cell.source.strip()

            if cell.cell_type == 'markdown':
                self.markdown += cell.source.split('\n')

        self.sections = [line for line in self.markdown if len(line.strip()) and line.strip()[0] == '#']

        # Save the notebook.
        jupytext.write(self.notebook, output)

        sys.stderr.write(' done.\n')
        sys.stderr.write('    contents:%s\n' % '\n        '.join(['']+self.sections))
        sys.stderr.flush()

Notebook(['src/preamble.py', 'src/solvents.py'], '01_solvents.ipynb')
Notebook(['src/preamble.py', 'src/gas.py'],      '02_gas.ipynb')
Notebook(['src/preamble.py', 'src/droplet.py'],  '03_droplet.ipynb')
Notebook(['src/preamble.py', 'src/gui.py'],      '04_gui.ipynb')

Notebook(['src/preamble.py',
          'src/solvents.py',
          'src/gas.py',
          'src/droplet.py',
          'src/gui.py'], '05_everything.ipynb')

Notebook(['src/preamble.py', 'src/benchmarking.py'], '06_benchmarking.ipynb')
