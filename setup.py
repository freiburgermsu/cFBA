# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

with open('README.rst') as file:
    readme = file.read()

setup(
  name = 'dFBApy',      
  package_dir = {'dfba':'dfbapy'},
  packages = find_packages(),
  package_data = {
	'dfbapy':['*'],
    'test': ['*']
  },
  version = '0.0.3',
  license = 'MIT',
  description = "Simulates dynamic Flux Balance Analysis (dFBA) of a BiGG-formatted XML model and with kinetics data in the format of the BiGG-SABIO module.", 
  long_description = readme,
  author = 'Andrew Freiburger',               
  author_email = 'andrewfreiburger@gmail.com',
  url = 'https://github.com/freiburgermsu/dfbapy',   
  keywords = ['chemistry', 'biology', 'metabolism', 'metabolome', 'flux', 'balance', 'analysis', 'biochemistry', 'FBA', 'dynamic'],
  install_requires = ['scipy', 'cobra', 'sigfig', 'numpy', 'pandas', 'matplotlib']
)
