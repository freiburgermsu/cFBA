# -*- coding: utf-8 -*-

# import statements
from scipy.constants import micro, milli, nano, hour, minute, femto
from pprint import pprint
from datetime import date
from numpy import nan
from math import inf 
import cobra
import pandas
import warnings, json, re, os

   
# add the units of logarithm to the Magnesium concentration
def isnumber(string):
    try:
        string = str(string)
        string = string.strip()
        if re.sub('([0-9\.e-])', '',string) == '':
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
                 kinetics_path: str = None, # the path of the kinetics data JSON file
                 kinetics_data: dict = {}, # A dictionary of custom kinetics data
                 verbose = False,
                 printing = False,
                 jupyter = False):
        
        # define object content
        self.bigg_metabolites_ids = json.load(open(os.path.join(os.path.dirname(__file__), 'BiGG_metabolites, parsed.json')))
        self.bigg_metabolites_names = json.load(open(os.path.join(os.path.dirname(__file__), 'BiGG_metabolite_names, parsed.json')))
        self.model = cobra.io.read_sbml_model(bigg_model_path)
        self.verbose = verbose
        self.printing = printing
        self.jupyter = jupyter
        
        # define kinetics of the system
        if type(kinetics_path) is str:
            if not os.path.exists(kinetics_path):
                raise ValueError('The path {kinetics_data} is not a valid path')
            self.kinetics_data  = json.load(open(kinetics_path))
        if kinetics_data != {}:
            for data in kinetics_data:
                self.kinetics_data[data] = kinetics_data[data]
        
        # define the parameter and variable dictionaries
        self.parameters = {}
        self.parameters['bigg_model_name'] = os.path.basename(bigg_model_path)
        
        self.variables = {}
        self.variables['concentrations'] = {}
        self.variables['time_series'] = {}
        
        # define a time-series value for each metabolite in the model
        for metabolite in self.model.metabolites:
            self.variables['time_series'][metabolite.name] = []
            
#    def bigg_metabolite_name(self, bigg_id):
#        if 'bigg_name' in self.bigg_metabolites_ids[bigg_id]:
#            return self.bigg_metabolites_ids[bigg_id]['bigg_name']
#        return self.bigg_metabolites_ids[bigg_id]['name']
        
    def _find_data_match(self,enzyme, source):        
        # define the closest match of the data to the parameterized conditions
        if isnumber(self.kinetics_data[enzyme][source]['metadata']["Temperature"]):
            temperature_deviation = abs(self.parameters['temperature'] - float(self.kinetics_data[enzyme][source]['metadata']["Temperature"]))/self.parameters['temperature']
        else:
            temperature_deviation = 0
            
        if isnumber(self.kinetics_data[enzyme][source]['metadata']["pH"]):
            ph_deviation = abs(self.parameters['pH'] - float(self.kinetics_data[enzyme][source]['metadata']["pH"]))/self.parameters['pH']
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
    
    def _calculate_kinetics(self):        
        previous_column = f'{(self.timestep-1)*self.timestep_value} min'
#        parameter_values = {}
        for enzyme in self.kinetics_data:
            fluxes = []
            for source in self.kinetics_data[enzyme]: 
                uncalculable = False
                source_instance = self.kinetics_data[enzyme][source]
                if "substituted_rate_law" in source_instance:     # Statistics of the aggregated data with each condition should be provided in a separate file for provenance of the scraped content.
                    remainder = re.sub('([0-9A-Za-z/()e\-\+\.\*])', '', source_instance["substituted_rate_law"])
                    if remainder == '':
                        # define the concentrations of each variable value
                        conc_dict = {}
                        for var in self.kinetics_data[enzyme][source]['variables_name']:
                            name = source_instance['variables_name'][var]
                            if len(var) == 1:
                                if not name in self.concentrations.index:
                                    uncalculable = True
                                    break
                                conc_dict[var] = self.concentrations.at[name, self.col]   
                                
                        if uncalculable:
                            warnings.warn(f'MetaboliteError: The {name} variable is not described by the BiGG system')
                            break

                        print(conc_dict)
                        if conc_dict != {}:
                            print(enzyme)
                            print(source_instance["substituted_rate_law"])
                            locals().update(conc_dict)
                            flux = eval(source_instance["substituted_rate_law"])
                            print(flux)
                            
                            # possible combine with previous fluxes, depending upon its match with the desired temperature and pH conditions
                            add_or_write = self._find_data_match(enzyme, source)
                            if add_or_write == 'a':                                    
                                fluxes.append(flux) 
                            elif add_or_write == 'w':
                                fluxes = [flux]
                        else:
                            warnings.warn(f'MetaboliteError: The {enzyme} enzyme possesses chemical names that are not described in the BiGG model')
                    else:
                        warnings.warn('RateLawError: The rate law {} is not executable, as the consequence of these excessive characters: {}'.format(source_instance["substituted_rate_law"], remainder))
                else:    
                    print(f'RateLawError: The {source_instance} does not possess a rate law')
                        
            flux = average(fluxes)
            if isnumber(flux):
                self._set_constraints(enzyme, flux)
                self.fluxes.at[enzyme, self.col] = flux 
                if self.printing:
                    print(f'{enzyme} flux:', flux)
                    if average(fluxes) == 0:
                        print('fluxes:', fluxes)
                    print('\n')
            else:
                warnings.warn(f'FluxError: The {enzyme} reaction flux {source_instance["substituted_rate_law"]} value {flux} is not a number.')

    def _set_constraints(self, enzyme, flux):           
        # pass the constraint
        enzyme = self.defined_reactions[enzyme]
        enzyme_name = re.sub(' ', '_', enzyme.name)   
        if enzyme_name not in self.constrained:
            self.constrained.append(enzyme_name)
            
            constraint = self.model.problem.Constraint(enzyme.flux_expression, lb=flux, ub=flux, name=f'{enzyme_name}_kinetics')
            self.model.solver.update()
            self.model.add_cons_vars(constraint)
            print(flux, constraint)
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
                
    def _define_timestep(self, conc_indices, flux_indices):
        self.col = f'{self.timestep*self.timestep_value} min'
        self.concentrations[self.col] = [float(0) for ind in conc_indices]
        self.fluxes[self.col] = [nan for ind in flux_indices]
        print('timestep', self.timestep)
                            
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
        
        # define a list of metabolite ids
        model_ids = []
        model_names = []
        for met in self.model.metabolites:
            model_names.append(met.name)
            met = re.sub('(_.$)','',met.id)
            model_ids.append(met)
        
        # define the DataFrames
        self.col =  '0 min'
        conc_indices = set(met for met in model_names)
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
        for met in model_names:
            self.concentrations.at[str(met), self.col] = float(0)
        if type(initial_concentrations) is dict:
            for met in model_names:
                if met in initial_concentrations:
                    self.concentrations.at[met, self.col] = initial_concentrations[met]
        else:
            # incorporate initial content from the SABIO scraping
            for enzyme in self.kinetics_data:
                for condition in self.kinetics_data[enzyme]:
                    for var in self.kinetics_data[enzyme][condition]['initial_concentrations_M']:
                        bigg_id = None
                        if self.kinetics_data[enzyme][condition]['variables_name'][var] in self.bigg_metabolites_names:
                            bigg_id = self.bigg_metabolites_names[
                                    self.kinetics_data[enzyme][condition]['variables_name'][var]
                                    ]['id']
                        if bigg_id in model_ids:
                            self.concentrations.at[str(met), self.col] = self.kinetics_data[enzyme][condition]['initial_concentrations_M'][var]/milli
                        else:
                            print(f"The {self.kinetics_data[enzyme][condition]['variables_name'][var]} metabolite is not in the BiGG model")
        
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
            # calculate custom fluxes, constrain the model, and update concentrations
            self._define_timestep(conc_indices, flux_indices)
            self._calculate_kinetics()                    
            self._execute_cobra()
            self._update_concentrations()
        
            self.variables['elapsed_time'] += self.timestep
            if self.printing:
                print(f'\nobjective value for timestep {self.timestep}: ', self.solutions[-1].objective_value)                
        
        # identify the chemicals that dynamically changed in concentrations
        for met in model_names:
            first = self.concentrations.at[met,'0 min']
            final = self.concentrations.at[met, self.col]
            if first != final:
                self.changed.add(met)
            if first == final:
                self.unchanged.add(met)
                
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