# dFBApy

A module foracquiring kinetic parameters and executing GEM dFBA is defined. 

### scraping
The scraping portion of the module determines the set of reactions that are defiend in a GEM from the BiGG repository, and accordingly scrapes the SABIO kinetic database for applicable expressions. These scraping efforts are organized into a JSON file, which can parameterize dFBA simulations.

### dynamic-FBA
The aforementioned kinetic parameters are used to execute dynamic FBA of a BiGG GEM. The core script is formulated around the Cobrapy module, and generates a descriptive dataframe of metabolic concentrations with each timestep that can be subsequently investigated and visualized. 
