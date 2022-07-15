#!/usr/bin/env python3

from setuptools import setup, Extension, find_packages
import subprocess

try:
    version = subprocess.check_output(["git", "describe", "--tags"]).decode('ascii').strip()
    if '-' in version:
        version, ncommits, current_commit = version.split('-')
        version = '%s.dev%s' % (version, ncommits)
except: version = 'unknown-version'

with open("README.md", "r") as f:
    long_description = f.read()

required_modules = ['numpy', 'scipy', 'matplotlib', 'pandas', 'jupytext', 'chemicals']
requirements = []
for module in required_modules:
    try: exec("import %s" % module)
    except: requirements += [module]

setup(
    name='sadkat',
    version=version,
    license='MIT License',

    author='Joshua F. Robinson and Dan Hardy',
    author_email='joshuarrr@protonmail.com',

    url='https://github.com/tranqui/sadkat',
    description='SADKAT: single-aerosol drying and trajectories',
    long_description=long_description,
    long_description_content_type="text/markdown",

    python_requires='>=3',
    ext_modules=[],
    install_requires=requirements,
    package_dir={'': 'src'},
    packages=find_packages('src'),

    classifiers=[
        "Programming Language :: Python :: 3",
    ],
 )
