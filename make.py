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

        # Show table of contents for generated notebook (from section headers).

        # Indent section headings based on level of section for prettier contents.
        indented_sections = ['    ' * (len(line) - len(line.lstrip('#')) - 1) + line for line in self.sections]
        sys.stderr.write('    contents:%s\n' % '\n        '.join([''] + indented_sections))
        sys.stderr.flush()

Notebook(['src/solvents.py'],     '01_solvents.ipynb')
Notebook(['src/solutes.py'],      '02_solutes.ipynb')
Notebook(['src/gas.py'],          '03_gas.ipynb')
Notebook(['src/droplet.py'],      '04_droplet.ipynb')
Notebook(['src/gui.py'],          '05_gui.ipynb')
Notebook(['src/benchmarking.py'], '06_benchmarking.ipynb')
Notebook(['src/equations.py'],    '07_equations.ipynb')

Notebook(['src/preamble.py',
          'src/solvents.py',
          'src/solutes.py',
          'src/gas.py',
          'src/droplet.py',
          'src/gui.py',
          'src/benchmarking.py',
          'src/equations.py'], '08_everything.ipynb', include_ignored=False)
