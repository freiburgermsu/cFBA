# import the libraries and the JSON data files        
import pandas as pd
import numpy
import json
import re
import cobra
import parser

sabio_reactions_df = pd.read_csv("sabio_out.csv")

# Opening JSON file
with open('sabio_entryids.json') as json_file:
    entry_id_data = json.load(json_file)

enzymenames = sabio_reactions_df["Enzymename"].unique().tolist()

enzyme_dict = {}

missing_entry_ids = []

parameters = {}

for enzyme in enzymenames:
    sabio_grouped_enzyme_df = sabio_reactions_df.loc[sabio_reactions_df["Enzymename"] == enzyme]

    
    dict_to_append = {}
    

    
    reactions = sabio_grouped_enzyme_df["Reaction"].unique().tolist()
        

            
    for reaction in reactions:
        
        dict_reactions_to_append = {}
         
        sabio_grouped_reactions_df = sabio_grouped_enzyme_df.loc[sabio_grouped_enzyme_df["Reaction"] == reaction]
        
        
    
        entryids = sabio_grouped_reactions_df["EntryID"].unique().tolist()
                
        for entryid in entryids:
 
            entry_ids_df = sabio_grouped_reactions_df.loc[sabio_grouped_reactions_df["EntryID"] == entryid]
            
            dict_entryid_to_append = {}
            
            head_of_df = entry_ids_df.head(1).squeeze()
            
            entry_id_flag = True
            
            parameter_info = {}
            
            try:
                parameter_info = entry_id_data[str(entryid)]
                dict_entryid_to_append["Parameters"] = parameter_info
            except:
                missing_entry_ids.append(str(entryid))
                entry_id_flag = False
                dict_entryid_to_append["Parameters"] = "NaN"
                
            rate_law = head_of_df["Rate Equation"]
            
            bad_rate_laws = ["unknown", "", "-"]
                
            if not rate_law in bad_rate_laws:                    
                dict_entryid_to_append["RateLaw"] = rate_law
                dict_entryid_to_append["SubstitutedRateLaw"] = rate_law
            else:
                dict_entryid_to_append["RateLaw"] = "NaN"
                dict_entryid_to_append["SubstitutedRateLaw"] = "NaN"
                                
            if entry_id_flag:
        
                fields_to_copy = ["Buffer", "Product", "PubMedID", "Publication", "pH", "Temperature", "Enzyme Variant", "UniProtKB_AC", "Organism", "KineticMechanismType", "SabioReactionID"]
                for field in fields_to_copy:  
                    dict_entryid_to_append[field] = head_of_df[field]
                dict_reactions_to_append[entryid] = dict_entryid_to_append
                
                dict_entryid_to_append["Substrates"] = head_of_df["Substrate"].split(";")
                
                
                out_rate_law = rate_law
            
                
                if not rate_law in bad_rate_laws:                    
                    substrates = head_of_df["Substrate"].split(";")
                    
                    to_strip = "1234567890"
                    stripped_string = rate_law
                    for character in to_strip:
                        stripped_string = stripped_string.replace(character, "")
                                            
                    
                    variables = re.split("\^|\*|\+|\-|\/|\(|\)| ", stripped_string)
                    variables = ' '.join(variables).split()
                    
                    start_value_permutations = ["start value", "start val."]
                    
                    substrates_key = {}
                    

                    
                    for var in variables:
                        if var in parameter_info:
                            
                            for permutation in start_value_permutations:
                                try:
                                    if var == "A" or var == "B":
                                        substrates_key[var] = parameter_info[var]["species"]
                                    else:
                                        value = parameter_info[var][permutation]
                                        if value != "-" and value != "" and value != " ":
                                            out_rate_law = out_rate_law.replace(var, parameter_info[var][permutation])
                                except:
                                    pass
                                
                    dict_entryid_to_append["RateLawSubstrates"] = substrates_key
                    

   
                    dict_entryid_to_append["SubstitutedRateLaw"] = out_rate_law
                    
                    '''
                    for var in variables:
                        if var not in parameters:
                            parameters[var] = 1
                        else:
                            parameters[var] += 1
                    '''
    
            
                    
        dict_to_append[reaction] = dict_reactions_to_append
                    
    enzyme_dict[enzyme] = dict_to_append
    
with open('sabio_proccessed.json', 'w') as f:
    json.dump(enzyme_dict, f, indent=4)
            
