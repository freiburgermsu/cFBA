# dFBApy

A module foracquiring kinetic parameters and executing GEM dFBA is defined. 

### scraping
The scraping section scrapes the SABIO kinetic database for the reactions that are defined in a GEM model. These scraping efforts are organized into a JSON file, which can parameterized in dFBA simulations.

### dynamic-FBA
The aforementioned kinetic parameters are parameterized into a dynamic FBA simulation of the corresponding GEM. The core script is formulated around Cobrapy, and generates a descriptive dataframe of metabolic concentrations with each timestep that can be externally investigated and visualized. 
