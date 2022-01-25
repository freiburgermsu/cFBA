# -*- coding: utf-8 -*-

# import statements
from scipy.constants import micro, milli, nano, hour, minute, femto
from math import inf 
import cobra
import pandas
import json
import re

   
# add the units of logarithm to the Magnesium concentration
def isnumber(string):
    string = str(string)
    try:
        string = string.strip()
        if re.sub('([0-9]|\.)', '',string) == '':
            return True
    except:
        return False
    
def average(num_1, num_2 = None):
    if isnumber(num_1): 
        if num_2 is not None:
            numbers = [num_1, num_2]
            average = sum(numbers) / len(numbers)
            return average
        else:
            return num_1
    elif type(num_1) is list:
        average = None
        summation = total = 0
        for num in num_1:
            if num is not None:
                summation += num
                total += 1
        if total > 0:
            average = summation/total
        return average
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
                 printing = False):
        
        # define object content
        self.bigg_metabolites = json.load(open('BiGG_metabolites, parsed.json'))
        self.model = cobra.io.read_sbml_model(bigg_model_path)
        self.kinetics_data = kinetics_data
        self.verbose = verbose
        self.printing = printing
        
        # define the parameter and variable dictionaries
        self.parameters = {}
        
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
    
    def _calculate_kinetics(self):
        def determine_parameter_value(param):
            chemical = source_instance["Parameters"][param]["chemical"]
            unit = source_instance["Parameters"][param]["unit"]
            if chemical in self.concentrations.index: 
                param_val = self.concentrations.at[chemical, previous_column]
                if re.search('m', unit):
                    param_val *= milli
                elif re.search('n', unit):
                    param_val *= nano
                elif re.search('u|U\+00B5', unit):
                    param_val *= micro        
            else:
                print(f'-->ERROR: The A parameter {chemical} is named differently in the BiGG model')
            return param_val
        
        previous_column = f'{(self.timestep-1)*self.timestep_value} min'
        for enzyme in self.kinetics_data:
            fluxes = []
            for source in self.kinetics_data[enzyme]: 
                source_instance = self.kinetics_data[enzyme][source]
                if "SubstitutedRateLaw" in source_instance:     # Statistics of the aggregated data with each condition should be provided in a separate file for provenance of the scraped content.
                    remainder = re.sub('([0-9ABC/()+.*millimicro])', '', source_instance["SubstitutedRateLaw"])
                    if remainder == '':
                        A = B = C = 0
                        if "A" in source_instance["Parameters"]:  
                            A = determine_parameter_value('A')
                        if "B" in source_instance["Parameters"]:   
                            B = determine_parameter_value('B')
                        if "C" in source_instance["Parameters"]:   
                            C = determine_parameter_value('C')

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


            if fluxes != [] and float(average(fluxes)) != 0:
                if self.printing:
                    print(f'{enzyme} fluxes:', fluxes, '\n')
                self.fluxes.at[enzyme, self.col] = average(fluxes)

    def _set_constraints(self, reaction, constant):           
        # pass the constraint
        reaction_name = re.sub(' ', '_', reaction.name)   
        if self.timestep == 1 and reaction_name not in self.constrained:
            self.constrained.append(reaction_name)
            
            constraint = self.model.problem.Constraint(reaction.flux_expression, lb=constant, ub=constant, name=f'{reaction_name}_kinetics')
            self.model.solver.update()
            self.model.add_cons_vars(constraint)
            self.model.solver.update()
        else:
            reaction.upper_bound = reaction.lower_bound = constant
                
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
        self.cell_dry_mass = 0.2*femto  # Teemu P. Miettinen, Kevin S. Ly, Alice Lam, Scott R. Manalis, https://doi.org/10.1101/2021.12.30.474524 
        self.cell_liters= 1*femto       # Lewis, C. L., Craig, C. C., & Senecal, A. G. (2014). https://doi.org/10.1128/AEM.00117-14
        self.changed = set()
        self.unchanged = set()
        self.minimum = inf
        self.constrained = []
        
        # define a concentrations DataFrame
        indices= set(met.name for met in self.model.metabolites)
        columns = [f'{(step)*self.timestep_value} min' for step in range(self.parameters['timesteps']+1)]
        self.concentrations = pandas.DataFrame(index = indices, columns = columns)
        self.concentrations.index.name = 'metabolite (\u0394mM)'
        
        # define a fluxes DataFrame
        indices= set(rxn.name for rxn in self.model.reactions)
        columns = [f'{(step)*self.timestep_value} min' for step in range(self.parameters['timesteps']+1)]
        self.fluxes = pandas.DataFrame(index = indices, columns = columns)
        self.fluxes.index.name = 'enzymes (mmol/g_(dw)/hr)'
        
        # define experimental conditions
        self.parameters['temperature'] = temperature
        self.parameters['pH'] = p_h
        self.variables['elapsed_time'] = 0
        
        # assign the initial concentrations
        for met in self.model.metabolites:
            self.concentrations.at[str(met.name), '0 min'] = 0
            if type(initial_concentrations) is dict:
                if met.name in initial_concentrations:
                    self.concentrations.at[str(met.name), '0 min'] = initial_concentrations[str(met.name)]
                
        # determine the BiGG reactions for which kinetics are predefined
        self.defined_reactions = []
        for rxn in self.model.reactions:
            if rxn.name in self.kinetics_data:
                self.defined_reactions.append(rxn)
        
        #Simulate for each timestep within total time frame
        solutions = []
        for self.timestep in range(1,self.parameters['timesteps']+1):
            print('timestep', self.timestep)
            self.col = f'{self.timestep*self.timestep_value} min'
            
            # execute the COBRA model 
            solution = self.model.optimize()
            solutions.append(solution)
            print(f'\nobjective value for timestep {self.timestep}: ', solution.objective_value)
            for rxn in self.model.reactions:
                self.fluxes.at[rxn.name, self.col] = solution.fluxes[rxn.id]
            
            # calculate custom fluxes, constrain the model, and update concentrations
            self._calculate_kinetics()
            self._update_concentrations()
            for reaction in self.defined_reactions:
                self._set_constraints(reaction, self.fluxes.at[reaction.name, self.col])       

            self.variables['elapsed_time'] += self.timestep
            
        # identify the chemicals that dynamically changed in concentrations
        for met in self.model.metabolites:
            first = self.concentrations.at[str(met.name), '0 min']
            final = self.concentrations.at[str(met.name), self.col]
            if first != final:
                self.changed.add(met.name)
            if first == final:
                self.unchanged.add(met.name)
                
        if self.verbose:
            print('\n\nUnchanged metabolites', '\n', '='*2*len('unchanged metabolites'), '\n', self.unchanged)
            print('\n\nChanged metabolites', '\n', '='*2*len('changed metabolites'), '\n', self.changed)
            print(f'\nThe {self.constrained} reactions were constrained in the COBRA model.')     
        elif self.printing:
            if self.unchanged == set():
                print('\n-->All of the metabolites changed concentration over the simulation')
            else:
                print('\n\nUnchanged metabolites', '\n', '='*2*len('unchanged metabolites'), '\n', self.unchanged)
                
        return solutions
        
    def export(self,):
        pass