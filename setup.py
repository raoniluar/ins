#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 15:25:06 2022

@author: raoni
"""
from setuptools import setup, find_packages

with open('README.md', encoding='utf-8') as fh:
    long_description = fh.read()

setup(
      name='ins',
      version='1.0',
      description='Calculate Index of Nonstationarity',
      author='Raoni Alcantara',
      author_email='raoni@ime.eb.br',
      url='https://github.com/raoniluar/ins',
      long_description=long_description,
      long_description_content_type='text/markdown',
      license='',
      install_requires=['numpy', 'scipy', 'matplotlib'],
      classifiers=[
          'Intended Audience :: Science/Research'
          ],
      packages=find_packages()
      )