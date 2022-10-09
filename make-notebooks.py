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
            source += f.readlines() + ['\n']

    source = ''.join(source)
    return source

class Notebook:
    """Generate a jupyter notebook for the end-user by merging python source files together."""

    def __init__(self, input_paths, output, include_ignored=True):
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
        if not include_ignored:
            self.notebook.cells = list(filter(lambda cell: 'ignore' not in cell.metadata or
                                              cell.metadata['ignore'] == 'False', self.notebook.cells))

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

        # Show table of contents for generated notebook (from section headers).

        # Indent section headings based on level of section for prettier contents.
        indented_sections = ['    ' * (len(line) - len(line.lstrip('#')) - 1) + line for line in self.sections]
        sys.stderr.write('    contents:%s\n' % '\n        '.join([''] + indented_sections))
        sys.stderr.flush()

source_dir = 'src/sadkat'

Notebook(['%s/solvents.py' % source_dir],     '01_solvents.ipynb')
Notebook(['%s/solutes.py' % source_dir],      '02_solutes.ipynb')
Notebook(['%s/gas.py' % source_dir],          '03_gas.ipynb')
Notebook(['%s/droplet.py' % source_dir],      '04_droplet.ipynb')
Notebook(['%s/gui.py' % source_dir],          '05_gui.ipynb')
Notebook(['%s/benchmarking.py' % source_dir], '06_benchmarking.ipynb')
Notebook(['%s/equations.py' % source_dir],    '07_equations.ipynb')

Notebook(['%s/common.py' % source_dir,
          '%s/solvents.py' % source_dir,
          '%s/solutes.py' % source_dir,
          '%s/gas.py' % source_dir,
          '%s/droplet.py' % source_dir,
          '%s/gui.py' % source_dir,
          '%s/benchmarking.py' % source_dir,
          '%s/equations.py' % source_dir], '08_everything.ipynb', include_ignored=False)
