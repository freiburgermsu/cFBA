# -*- coding: utf-8 -*-

# import statements
from scipy.constants import micro, milli, nano, hour, minute, femto
from datetime import date
from numpy import nan
from math import inf 
import cobra
import pandas
import json, re, os

   
# add the units of logarithm to the Magnesium concentration
def isnumber(string):
    try:
        string = str(string)
        string = string.strip()
        if re.sub('([0-9\.])', '',string) == '':
            return True
    except:
        return False
    
def average(num_1, num_2 = None):
    if isnumber(num_1): 
        if isnumber(num_2):
            numbers = [num_1, num_2]
            return sum(numbers) / len(numbers)
        else:
            return num_1
    elif type(num_1) is list:
        summation = total = 0
        for num in num_1:
            if num is not None:
                summation += num
                total += 1
        if total > 0:
            return summation/total
        return None
    else:
        return None
    
#    export_directory: Optional[str] = None, export_name: Optional[str] = None
#    ) -> None
            
# define chemical concentrations
class dFBA():
    def __init__(self, 
                 bigg_model_path,    # BiGG model            
                 kinetics_data: dict, # unsure about this structure 
                 verbose = False,
                 printing = False,
                 jupyter = False):
        
        # define object content
        self.bigg_metabolites = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'scraping', 'SABIO', 'BiGG_metabolites, parsed.json')))
        self.model = cobra.io.read_sbml_model(bigg_model_path)
        self.kinetics_data = kinetics_data
        self.verbose = verbose
        self.printing = printing
        self.jupyter = jupyter
        
        # define the parameter and variable dictionaries
        self.parameters = {}
        self.parameters['bigg_model_name'] = os.path.basename(bigg_model_path)
        
        self.variables = {}
        self.variables['concentrations'] = {}
        self.variables['time_series'] = {}
        
        # define a time-series value for each metabolite in the model
        for metabolite in self.model.metabolites:
            self.variables['time_series'][metabolite.name] = []
            
    def bigg_metabolite_name(self, bigg_id):
        if 'bigg_name' in self.bigg_metabolites[bigg_id]:
            return self.bigg_metabolites[bigg_id]['bigg_name']
        return self.bigg_metabolites[bigg_id]['name']
        
    def _find_data_match(self,enzyme, source):        
        # define the closest match of the data to the parameterized conditions
        if isnumber(self.kinetics_data[enzyme][source]["Temperature"]):
            temperature_deviation = abs(self.parameters['temperature'] - float(self.kinetics_data[enzyme][source]["Temperature"]))/self.parameters['temperature']
        else:
            temperature_deviation = 0
            
        if isnumber(self.kinetics_data[enzyme][source]["pH"]):
            ph_deviation = abs(self.parameters['pH'] - float(self.kinetics_data[enzyme][source]["pH"]))/self.parameters['pH']
        else:
            ph_deviation = 0

        old_minimum = self.minimum
        deviation = average(temperature_deviation, ph_deviation)
        self.minimum = min(deviation, self.minimum)
#        print('minimum', minimum)
#        print('deviation', deviation)

        if old_minimum == self.minimum:
            return 'a'
        elif deviation == self.minimum:
            return 'w'
        
    def _determine_parameter_value(self, unit, value, param):
        if re.search('m', unit):
            if re.search('(\/m|\/\(m)', unit):
                value /= milli
            else:
                value *= milli
        elif re.search('n', unit):
            if re.search('(\/n|\/\(n)', unit):
                value /= milli
            else:
                value *= nano
        elif re.search('u|U\+00B5', unit):
            if re.search('(\/u|\/U\+00B5|\/\(m|\/\(U\+00B5)', unit):
                value /= micro
            else:
                value *= micro
        else:
            print(f'-->The parameter {param} is possess differently in the BiGG model')
        return value
    
    def _calculate_kinetics(self):        
        previous_column = f'{(self.timestep-1)*self.timestep_value} min'
#        parameter_values = {}
        for enzyme in self.kinetics_data:
            fluxes = []
            for source in self.kinetics_data[enzyme]: 
                
                source_instance = self.kinetics_data[enzyme][source]
                if "SubstitutedRateLaw" in source_instance:     # Statistics of the aggregated data with each condition should be provided in a separate file for provenance of the scraped content.
                    remainder = re.sub('([0-9A-Z/()+.*millimicro])', '', source_instance["SubstitutedRateLaw"])
                    if remainder == '':
                        for param in source_instance["Parameters"]:
                            chemical = source_instance["Parameters"][param]["chemical"]
                            unit = source_instance["Parameters"][param]["unit"]
                            value = source_instance["Parameters"][param]["value"]
                            
                            if param == "A":  
                                A = self._determine_parameter_value(unit, value, chemical)
                            elif param == "B":  
                                B = self._determine_parameter_value(unit, value, chemical)
                            elif param == "C":  
                                C = self._determine_parameter_value(unit, value, chemical)
                            elif param == "D":  
                                D = self._determine_parameter_value(unit, value, chemical)
                            elif param == "S":  
                                S = self._determine_parameter_value(unit, value, chemical)
                            elif param == "I":  
                                I = self._determine_parameter_value(unit, value, chemical)
                            else:
                                print('-->ERROR: The {param} parameter, for the {chemical} chemical, is not captured by < _calculate_kinetics() > function.')
                                continue
                                
                            try:
                                flux = eval(source_instance["SubstitutedRateLaw"])
        #                                print(A, B, source_instance["SubstitutedRateLaw"])
                            except:
                                print('-->ERROR: The kinetic expression {} is not valid'.format(source_instance["SubstitutedRateLaw"]))
                                flux = 0
                                pass
                            
                            add_or_write = self._find_data_match(enzyme, source)
                            if add_or_write == 'a':                                    
                                fluxes.append(flux) 
                            elif add_or_write == 'w':
                                fluxes = [flux]
                    else:
                        print('-->ERROR: The rate law {} is not executable, as the consequence of these excessive characters: {}'.format(source_instance["SubstitutedRateLaw"], remainder))

            flux = average(fluxes)
            self._set_constraints(enzyme, flux)
            self.fluxes.at[enzyme, self.col] = flux 
            if self.printing:
                print(f'{enzyme} flux:', flux)
                if average(fluxes) == 0:
                    print('fluxes:', fluxes)
                print('\n')
                

    def _set_constraints(self, enzyme, flux):           
        # pass the constraint
        enzyme = self.defined_reactions[enzyme]
        enzyme_name = re.sub(' ', '_', enzyme.name)   
        if enzyme_name not in self.constrained:
            self.constrained.append(enzyme_name)
            
            constraint = self.model.problem.Constraint(enzyme.flux_expression, lb=flux, ub=flux, name=f'{enzyme_name}_kinetics')
            self.model.solver.update()
            self.model.add_cons_vars(constraint)
            self.model.solver.update()
        else:
            if flux > self.model.constraints[f'{enzyme_name}_kinetics'].ub:
                self.model.constraints[f'{enzyme_name}_kinetics'].ub = flux
                self.model.constraints[f'{enzyme_name}_kinetics'].lb = flux
            else:                
                self.model.constraints[f'{enzyme_name}_kinetics'].lb = flux
                self.model.constraints[f'{enzyme_name}_kinetics'].ub = flux
            
        if self.printing:
            print(self.model.constraints[f'{enzyme_name}_kinetics'])
                
    def _update_concentrations(self):
        for met in self.model.metabolites:      
            self.concentrations.at[str(met.name), self.col] = 0
            for rxn in self.model.reactions:
                if met in rxn.metabolites:    # flux units: mmol/(g_(dry weight)*hour)
                    stoich = rxn.metabolites[met]
                    flux = self.fluxes.at[str(rxn.name), self.col] 
                    delta_conc = stoich * ((flux * self.timestep_value*(minute/hour) * self.cell_dry_mass) / self.cell_liters) 
                    self.concentrations.at[str(met.name), self.col] += delta_conc
#                    print(met, rxn.metabolites[met], self.fluxes.at[str(rxn.name), col], delta_conc)
                    
    def _execute_cobra(self):
        # execute the COBRA model 
        solution = self.model.optimize()
        self.solutions.append(solution)
        for rxn in self.model.reactions:
            if not isnumber(self.fluxes.at[rxn.name, self.col]):
                self.fluxes.at[rxn.name, self.col] = solution.fluxes[rxn.id]
                            
    def simulate(self, 
                 total_time: float,  # mintues
                 timestep: float,  # minutes
                 initial_concentrations: dict = None,
                 temperature = 25, 
                 p_h = 7, 
                 visualize = True
                 ):
        
        # define the dataframe for the time series content
        self.parameters['timesteps'] = round(total_time/timestep)
        self.timestep_value = timestep
        self.total_time = total_time
        self.cell_dry_mass = 0.2*femto  # Teemu P. Miettinen, Kevin S. Ly, Alice Lam, Scott R. Manalis, https://doi.org/10.1101/2021.12.30.474524 
        self.cell_liters= 1*femto       # Lewis, C. L., Craig, C. C., & Senecal, A. G. (2014). https://doi.org/10.1128/AEM.00117-14
        self.changed = set()
        self.unchanged = set()
        self.minimum = inf
        self.constrained = []
        self.solutions = []
        
        # define the DataFrames
        self.col =  '0 min'
        conc_indices = set(met.name for met in self.model.metabolites)
        self.concentrations = pandas.DataFrame(index = conc_indices, columns = [self.col])
        self.concentrations.index.name = 'metabolite (\u0394mM)'
        
        flux_indices = set(rxn.name for rxn in self.model.reactions)
        self.fluxes = pandas.DataFrame(index = flux_indices, columns = [self.col])
        self.fluxes.index.name = 'enzymes (mmol/g_(dw)/hr)'
        
        # define experimental conditions
        self.parameters['temperature'] = temperature
        self.parameters['pH'] = p_h
        self.variables['elapsed_time'] = 0
        
        # assign the initial concentrations
        for met in self.model.metabolites:
            self.concentrations.at[str(met.name), self.col] = float(0)
            if type(initial_concentrations) is dict:
                if met.name in initial_concentrations:
                    self.concentrations.at[str(met.name), self.col] = initial_concentrations[str(met.name)]
        
        #Simulate for each timestep within total time frame
        if type(self.kinetics_data) is dict:
            # determine the BiGG reactions for which kinetics are predefined
            self.defined_reactions = {}
            for rxn in self.model.reactions:
                if rxn.name in self.kinetics_data:
                    self.defined_reactions[rxn.name] = rxn
                    
            # execute FBA for each timestep
            for self.timestep in range(1,self.parameters['timesteps']+1):
#                if self.timestep == 0:
#                    self._calculate_kinetics() 
#                    continue
                
                # expand the DataFrames
                self.col = f'{self.timestep*self.timestep_value} min'
                self.concentrations[self.col] = [float(0) for ind in conc_indices]
                self.fluxes[self.col] = [nan for ind in flux_indices]
                print('timestep', self.timestep)
                 
                # calculate custom fluxes, constrain the model, and update concentrations
                self._calculate_kinetics()                    
                self._execute_cobra()
                self._update_concentrations()    
        
                self.variables['elapsed_time'] += self.timestep
                
                if self.printing:
                    print(f'\nobjective value for timestep {self.timestep}: ', self.solutions[-1].objective_value)                
        else:
            self.col = f'{total_time} min'
            self.fluxes[self.col] = [nan for ind in flux_indices]
            self._execute_cobra()
            self._update_concentrations()
            if self.printing:
                print(f'\nobjective value: ', self.solutions[-1].objective_value)
            
        # identify the chemicals that dynamically changed in concentrations
        for met in self.model.metabolites:
            first = self.concentrations.at[str(met.name), '0 min']
            final = self.concentrations.at[str(met.name), self.col]
            if first != final:
                self.changed.add(met.name)
            if first == final:
                self.unchanged.add(met.name)
                
        if self.verbose:
            print('\n\nUnchanged metabolite concentrations', '\n', '='*2*len('unchanged metabolites'), '\n', self.unchanged)
            print('\n\nChanged metabolite  concentrations', '\n', '='*2*len('changed metabolites'), '\n', self.changed)
            print(f'\nThe {self.constrained} reactions were constrained in the COBRA model.')     
        elif self.printing:
            if self.jupyter:
                pandas.set_option('max_rows', None)
                display(self.concentrations)
                display(self.fluxes)
            if self.unchanged == set():
                print('\n-->All of the metabolites changed concentration over the simulation')
            else:
                print('\n\nUnchanged metabolite concentrations', '\n', '='*2*len('unchanged metabolites'), '\n', self.unchanged)
        
    def export(self, export_name = None, export_directory = None):
        # define a unique simulation name 
        if export_name is None:
            export_name = '-'.join([re.sub(' ', '_', str(x)) for x in [date.today(), 'dFBA', self.parameters['bigg_model_name'], f'{self.total_time} min']])
        if export_directory is None:
            directory = os.getcwd()
        else:
            directory = os.path.dirname(export_directory)
            
        simulation_number = -1
        while os.path.exists(os.path.join(directory, export_name)):
            simulation_number += 1
            export_name = re.sub('(\-\d+$)', '', export_name)
            export_name = '-'.join([export_name, str(simulation_number)])
            
        self.simulation_path = os.path.join(directory, export_name)
        self.parameters['simulation_path'] = self.simulation_path
        os.mkdir(self.simulation_path)
        
        # export the content to the simulation folder
        self.fluxes.to_csv(os.path.join(self.simulation_path, 'fluxes.csv'))
        self.concentrations.to_csv(os.path.join(self.simulation_path, 'concentrations.csv'))
        
        times = self.fluxes.columns
        with open(os.path.join(self.simulation_path, 'objective_values.txt'), 'w') as obj_val:            
            for sol in self.solutions:
                index = self.solutions.index(sol)
                time = times[index]
                obj_val.write(f'{time}: {sol.objective_value}')              
        
        # export the parameters
        parameters = {'parameter':[], 'value':[]}
        for parameter in self.parameters:
            parameters['parameter'].append(parameter)
            parameters['value'].append(self.parameters[parameter])
            
        parameters_table = pandas.DataFrame(parameters)
        parameters_table.to_csv(os.path.join(self.simulation_path, 'parameters.csv'))