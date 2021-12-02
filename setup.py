#!/usr/bin/env python

from distutils.core import setup
from setuptools import find_packages

# parse requirements.txt
required = []

with open('requirements.txt') as f:
    lines = f.read().splitlines()
    for item in lines:
        if item.startswith('#'):
            continue
        elif item.startswith('-e'):
            continue
        else:
            required.append(item)

setup(name='sbml_wrapper',
      version='1.0',
      description='SBML Wrapper',
      author='Janosch Brandhorst',
      package_dir={'sbml_wrapper': 'src/sbml_wrapper'},
      packages=["sbml_wrapper"],
      package_data={
          '': ['./requirements.txt'],
      },
    include_package_data=True,
    python_requires='>=3.7',
     )