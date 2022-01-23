# -*- coding: utf-8 -*-

# import statements
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
        average = sum(num_1) / len(num_1)
        return average
    else:
        return None
    
#    export_directory: Optional[str] = None, export_name: Optional[str] = None
#    ) -> None
            
# define chemical concentrations
class dFBA():
    def __init__(self, 
                 model_path,                            # BiGG model            
                 reaction_kinetics: dict, # unsure about this structure 
                 verbose = False,
                 printing = False):
        
        # define object content
        self.model = cobra.io.read_sbml_model(model_path)
        self.reaction_kinetics = reaction_kinetics
        self.verbose = verbose
        self.printing = printing
        
        # define the parameter and variable dictionaries
        self.parameters = {}
        self.parameters['calculated_rate_laws'] = {}
        
        self.variables = {}
        self.variables['concentrations'] = {}
        self.variables['time_series'] = {}
        
        # define a time-series value for each metabolite in the model
        for metabolite in self.model.metabolites:
            self.variables['time_series'][metabolite.name] = []
        
    def _find_data_match(self,enzyme, reaction, entry_id):        
        # define the closest match of the data to the parameterized conditions
        minimum = inf
        if isnumber(self.reaction_kinetics[enzyme][reaction][entry_id]["Temperature"]):
            temperature_deviation = abs(self.parameters['temperature'] - float(self.reaction_kinetics[enzyme][reaction][entry_id]["Temperature"]))/self.parameters['temperature']
        else:
            temperature_deviation = 0
            
        if isnumber(self.reaction_kinetics[enzyme][reaction][entry_id]["pH"]):
            ph_deviation = abs(self.parameters['pH'] - float(self.reaction_kinetics[enzyme][reaction][entry_id]["pH"]))/self.parameters['pH']
        else:
            ph_deviation = 0

        old_minimum = minimum
        deviation = average(temperature_deviation, ph_deviation)
        minimum = min(deviation, minimum)
#        print('minimum', minimum)
#        print('deviation', deviation)

        if old_minimum == minimum:
            return 'a'
        elif deviation == minimum:
            return 'w'
    
    def _calculate_kinetics(self, concentrations: dict):
        for enzyme in self.reaction_kinetics:
            fluxes = []
            for reaction in self.reaction_kinetics[enzyme]:   
                for condition in self.reaction_kinetics[enzyme][reaction]:    
                    entry = False
                    condition_instance = self.reaction_kinetics[enzyme][reaction][condition]
                    if "SubstitutedRateLaw" in condition_instance:     # Statistics of the aggregated data with each condition should be provided in a separate file for provenance of the scraped content.
                        remainder = re.sub('([0-9ABC/()+.*])', '', condition_instance["SubstitutedRateLaw"])
                        if remainder == '':
                            A = 0
                            B = 0
                            C = 0

                            if "A" in condition_instance["Parameters"]:   
                                if condition_instance["Parameters"]["A"]["species"] in concentrations: 
                                    A = concentrations[condition_instance["Parameters"]["A"]["species"]]

                            if "B" in condition_instance["Parameters"]:   
                                if condition_instance["Parameters"]["B"]["species"] in concentrations:    
                                    B = concentrations[condition_instance["Parameters"]["B"]["species"]]

                            if "C" in condition_instance["Parameters"]:   
                                if condition_instance["Parameters"]["C"]["species"] in concentrations:    
                                    C = concentrations[condition_instance["Parameters"]["C"]["species"]]

                            try:
                                flux = eval(condition_instance["SubstitutedRateLaw"])
                            except:
                                pass
                            
                            add_or_write = self._find_data_match(enzyme, reaction, entry_id)
                            if add_or_write == 'a':                                    
                                fluxes.append(flux) 
                            elif add_or_write == 'w':
                                fluxes = [flux]
                                
                            entry = True

            if entry:
                if self.verbose:
                    print(enzyme, self.timestep, fluxes)
                self.df.at[enzyme, f'{self.timestep} (min)'] = average(fluxes)

    def _set_constraints(self, reaction, reaction_name, constant, first):
        reaction_name = re.sub('\s', '_', reaction_name)
        if first:
            expression = reaction.flux_expression
            constraint = self.model.problem.Constraint(expression, lb=constant, ub=constant, name=f'{reaction_name}_kinetics')
            self.model.add_cons_vars(constraint)
            self.model.solver.update()
        else:
            if not isnumber(constant):
                print(f'--> ERROR: The constant for {reaction_name} is erronenous.')
            else:
                reaction.upper_bound = reaction.lower_bound = constant
                
    def _update_concentrations(self, solution, first):
        changed = set()
        unchanged = set()
        for metabolite in self.model.metabolites:
            # define the initial concentrations
            if first and metabolite.name in self.parameters['initial_concentrations']:
                before = self.parameters['initial_concentrations'][metabolite.name]
            elif first:
                before = 0
            elif not first: 
                before = self.variables['concentrations'][metabolite.name]
            
            # calculate the change in concentrations from the reaction fluxes
            for reaction in self.model.reactions:
                if metabolite in reaction.metabolites:
                    delta_conc = reaction.metabolites[metabolite] * solution.fluxes[reaction.id]
                    if first:
                        self.variables['concentrations'][metabolite.name] = delta_conc
                    else:
                        self.variables['concentrations'][metabolite.name] += delta_conc
                        if before != self.variables['concentrations'][metabolite.name] and metabolite.name not in changed:
                            changed.add(metabolite.name)
                            if self.verbose:
                                print('The {} concentration changed in timestep {}'.format(metabolite.name, self.timestep))
                        if before == self.variables['concentrations'][metabolite.name] and metabolite.name not in unchanged.union(changed):
                            unchanged.add(metabolite.name)
                            if self.printing:
                                print('--> The {} concentration did not change in timestep {}'.format(metabolite.name, self.timestep))
                            
    def simulate(self, 
                 initial_concentrations: dict,
                 total_time: float, 
                 timestep: float, 
                 temperature = 25, 
                 p_h = 7, 
                 visualize = True
                 ):
        
        # define the dataframe for the time series content
        self.parameters['timesteps'] = round(total_time/timestep)
        indices = [enzyme for enzyme in self.reaction_kinetics]
        columns = [f'{(step+1)*timestep} (min)' for step in range(self.parameters['timesteps'])]
        self.df = pandas.DataFrame(index = indices, columns = columns)
        self.df.index.name = 'enzymes'
        print(self.df)
        
        # define experimental conditions
        self.parameters['initial_concentrations'] = initial_concentrations
        self.parameters['temperature'] = temperature
        self.parameters['pH'] = p_h
        self.variables['elapsed_time'] = 0
        
        #Simulate for each timestep within total time frame
        first = True
        solutions = []
        for self.timestep in range(1,self.parameters['timesteps']+1):
            print('timestep', self.timestep)
            
            if first:
                self._calculate_kinetics(self.parameters['initial_concentrations'])
            else:
                self._calculate_kinetics(self.variables['concentrations'])
            
            #Calcuate and parameterize fluxes from michaelis-menten kinetics
            for reaction in self.model.reactions:
                if any(rxn.lower() == reaction.name for rxn in self.parameters['calculated_rate_laws']):
                    kinetic_flux = self.parameters['calculated_rate_laws'][reaction.name]
#                     print(kinetic_flux)
                    self._set_constraints(reaction, reaction.name,kinetic_flux, first)
                else:
                    if self.verbose and self.timestep == 1:
                        print(f'--> ERROR: The {reaction.name} reaction is not defined in the kinetics data.')                    

            # execute the model and update concentrations
            self.solution = self.model.optimize()
            solutions.append(self.solution)
            print(f'\nobjective value for timestep {self.timestep}: ', self.solution.objective_value, '\n\n')
            self._update_concentrations(self.solution, first)
            
            self.variables['elapsed_time'] += self.timestep
        
            if self.timestep == 1:
                first = False
        
        if self.printing:
            self._visualize()
        
        return self.df, solutions

    def _visualize(self):
        print('='*50)
        print(self.df.to_string())
        
    def export(self,):
        pass