# -*- coding: utf-8 -*-
"""
@authors: Ethan Sean Chan, Andrew Philip Freiburger
"""
from scipy.constants import minute, hour
from selenium.webdriver.support.ui import Select
from selenium import webdriver
from itertools import islice
from glob import glob
import datetime
import pandas
import json, time, re, os


class SABIO_scraping():
#     __slots__ = (str(x) for x in [progress_file_prefix, xls_download_prefix, scraped_xls_prefix, scraped_entryids_prefix, sel_xls_download_path, processed_xls, entry_json, scraped_model, bigg_model_name_suffix, sub_directory_path, progress_path, xls_download_path, scraped_xls_path, scraped_entryids_path, xls_csv_file_path, entryids_json_file_path, scraped_model_json_file_path, bigg_model, step_number, cwd])
    
    def __init__(self,
                 export_content: bool = False):
        self.export_content = export_content 
        self.count = 0
        self.parameters = {}
        self.parameters['general_delay'] = 2
        self.variables = {}
        self.paths = {}
        
        # load BiGG model content 
        self.bigg_metabolites = json.load(open('BiGG_metabolites, parsed.json'))
        self.bigg_reactions = json.load(open('BiGG_reactions, parsed.json'))

    #Clicks a HTML element with selenium by id
    def _click_element_id(self,n_id):
        element = self.driver.find_element_by_id(n_id)
        element.click()
        time.sleep(self.parameters['general_delay'])
        
    def _wait_for_id(self,n_id):
        while True:
            try:
                element = self.driver.find_element_by_id(n_id)
                break
            except:
                time.sleep(self.parameters['general_delay'])
        

    #Selects a choice from a HTML dropdown element with selenium by id
    def _select_dropdown_id(self,n_id, n_choice):
        element = Select(self.driver.find_element_by_id(n_id))
        element.select_by_visible_text(n_choice)
        time.sleep(self.parameters['general_delay'])

    def _fatal_error_handler(self,message):
        print("Error: " + message)
        print("Exiting now...")
        exit(0)

    """
    --------------------------------------------------------------------
        STEP 0: GET BIGG MODEL TO SCRAPE AND SETUP DIRECTORIES AND PROGRESS FILE
    --------------------------------------------------------------------    
    """

    def start(self,bigg_model_path,bigg_model_name): # find the BiGG model that will be scraped
        if os.path.exists(bigg_model_path):
            self.model = json.load(open(bigg_model_path))
            
        else:
            print('-> ERROR: The BiGG model file does not exist')
            
        self.bigg_model_name = bigg_model_name 
        if bigg_model_name is None:
            self.bigg_model_name = re.search("([\w+\.?\s?]+)(?=\.json)", bigg_model_path).group()

        # define the paths
        self.paths['cwd'] = os.path.dirname(os.path.realpath(bigg_model_path))
        self.paths['sub_directory_path'] = os.path.join(self.paths['cwd'],f"scraping-{self.bigg_model_name}")
        if not os.path.isdir(self.paths['sub_directory_path']):        
            os.mkdir(self.paths['sub_directory_path'])
                    
        self.variables['scraped_xls'] = {}
        self.variables['scraped_entryids'] = {}
        self.paths['scraped_model_json_file_path'] = os.path.join(self.paths['sub_directory_path'], "scraped_model") + ".json"
        self.paths['sel_xls_download_path'] = os.path.join(self.paths['sub_directory_path'],"downloaded_xls")
        
        self.paths['progress_path'] = os.path.join(self.paths['sub_directory_path'], "current_progress") + '.txt'
        if os.path.exists(self.paths['progress_path']):
            f = open(self.paths['progress_path'], "r")
            self.step_number = int(f.read(1))
            if not re.search('[1-5]',str(self.step_number)):
                self._fatal_error_handler("Progress file malformed. Please delete and restart")
        else:
            self.step_number = 1
            self._progress_update(self.step_number)

        self.paths['xls_download_path'] = os.path.join(self.paths['sub_directory_path'], 'downloaded_xls') 
        if not os.path.isdir(self.paths['xls_download_path']):
            os.mkdir(self.paths['xls_download_path'])

        self.paths['scraped_xls_path'] = os.path.join(self.paths['sub_directory_path'], "scraped_xls") + ".json"
        if os.path.exists(self.paths['scraped_xls_path']):
            f = open(self.paths['scraped_xls_path'], "r")
            self.variables['scraped_xls'] = json.load(f)
            f.close()

        self.paths['scraped_entryids_path'] = os.path.join(self.paths['sub_directory_path'], "scraped_entryids") + ".json"
        if os.path.exists(self.paths['scraped_entryids_path']):
            f = open(self.paths['scraped_entryids_path'], "r")
            self.variables['scraped_entryids'] = json.load(f)
            f.close()

        self.paths['entryids_progress_path'] = os.path.join(self.paths['sub_directory_path'], "entryids_progress") + ".json"
        if os.path.exists(self.paths['entryids_progress_path']):
            f = open(self.paths['entryids_progress_path'], "r")
            try:
                self.variables['entryids_progress'] = json.load(f)
            except:
                print('-> ERROR: The < entryids_progress.json > file is corrupted or empty.')
            f.close()
        else:
            self.variables['entryids_progress'] = {}
        
        # update the step counter
        self.step_number = 1
        self._progress_update(self.step_number)

    """
    --------------------------------------------------------------------
        STEP 1: SCRAPE SABIO WEBSITE BY DOWNLOAD XLS FOR GIVEN REACTIONS IN BIGG MODEL
    --------------------------------------------------------------------    
    """

    def scrape_xls(self,reaction_identifier, search_option):
        self.driver.get("http://sabiork.h-its.org/newSearch/index")
        self._wait_for_id("resetbtn")
        self._click_element_id("resetbtn")
        
        time.sleep(self.parameters['general_delay'])
        
        self._click_element_id("option")
        self._select_dropdown_id("searchterms", search_option)
        text_area = self.driver.find_element_by_id("searchtermField")
        text_area.send_keys(reaction_identifier)  
        
        time.sleep(self.parameters['general_delay']) 
        
        self._click_element_id("addsearch")
        
        time.sleep(self.parameters['general_delay'])

        result_num = ""
        try: 
            result_num_ele = self.driver.find_element_by_id("numberofKinLaw")
            for char in result_num_ele.text:
                if re.search('[0-9]', char):
                    result_num = result_num + char

            result_num = int(result_num)
        except:
            #self.driver.close()
            self.driver.get("http://sabiork.h-its.org/newSearch/index")
            return False

        time.sleep(self.parameters['general_delay'])

        self._select_dropdown_id("max", "100")
        element = Select(self.driver.find_element_by_id("max"))
        element.select_by_visible_text("100")

        time.sleep(self.parameters['general_delay'])

        if result_num > 0 and result_num <= 100:
            self._click_element_id("allCheckbox")
            time.sleep(self.parameters['general_delay'])
        elif result_num > 100:
            self._click_element_id("allCheckbox")
            for i in range(int(result_num/100)):
                element = self.driver.find_element_by_xpath("//*[@class = 'nextLink']")
                element.click()
                time.sleep(self.parameters['general_delay'])
                self._click_element_id("allCheckbox")
                time.sleep(self.parameters['general_delay'])
        else:
            #self.driver.close()
            self.driver.get("http://sabiork.h-its.org/newSearch/index")
            return False

        self.driver.get("http://sabiork.h-its.org/newSearch/spreadsheetExport")
        time.sleep(self.parameters['general_delay']*7.5)
        self._click_element_id("excelExport")
        time.sleep(self.parameters['general_delay']*2.5)
        #self.driver.close()

        return True
    
#    def _expand_shadow_element(self, element):
#        shadow_root = self.driver.execute_script('return arguments[0].shadowRoot', element)
#        return shadow_root
    
    def _split_reaction(self, reaction_string, sabio = False):
        def _parse_stoich(met):
            stoich = ''
            ch_number = 0
            denom = False
            while re.search('[0-9\./]', met[ch_number]): 
                stoich += met[ch_number]
                if met[ch_number] == '/':
                    numerator = stoich
                    denom = True
                if denom:
                    denominator += met[ch_number]
                ch_number += 1
                
            if denom:
                stoich = f'{numerator}/{denominator}'
            return stoich
        
        def met_parsing(met):
    #         print(met)
            met = met.strip()
            met = re.sub('_\w$', '', met)
            if re.search('(\d\s\w|\d\.\d\s|\d/\d\s)', met):
                coefficient = _parse_stoich(met)
                coefficient = '{} '.format(coefficient)
            else:
                coefficient = ''
            met = re.sub(coefficient, '', met)
    #         print(met, coefficient)
            return met, coefficient   
    
        def reformat_met_name(met_name, sabio = False):
            met_name = re.sub(' - ', '-', met_name)
            if not sabio:
                met_name = re.sub(' ', '_', met_name)
            return met_name
        
            
        # parse the reactants and products for the specified reaction string
        if not sabio:
            reaction_split = reaction_string.split('<->')
        else:
            reaction_split = reaction_string.split('=')
        reactants_list = reaction_split[0].split(' + ')
        products_list = reaction_split[1].split(' + ')
        
        # parse the reactants
        reactants = []
        sabio_reactants = []
        for met in reactants_list:
    #         print(met)
            met, coefficient = met_parsing(met)
            reactants.append(coefficient + reformat_met_name(self.bigg_metabolites[met]['name']))
            sabio_reactants.append(coefficient + reformat_met_name(self.bigg_metabolites[met]['name'], True))
    
        # parse the products
        products = []
        sabio_products = []
        for met in products_list:
            if not re.search('[a-z]', met, flags = re.IGNORECASE):
                continue
            met, coefficient = met_parsing(met)
            products.append(coefficient + reformat_met_name(self.bigg_metabolites[met]['name']))
            sabio_products.append(coefficient + reformat_met_name(self.bigg_metabolites[met]['name'], True))
    
    #     compounds = reactants + products
        reactant_string = ' + '.join(reactants)
        product_string = ' + '.join(products)
        if not sabio:
            reaction_string = ' <-> '.join([reactant_string, product_string])
        else:
            reaction_string = ' = '.join([reactant_string, product_string])
        
        # construct the set of compounds in the SABIO format
        sabio_compounds = sabio_reactants + sabio_products
        
        return reaction_string, sabio_compounds
    
    
        return minutes_per_enzyme*reactions_quantity

    def scrape_bigg_xls(self,):
        """
        chrome_options = webdriver.ChromeOptions()
        prefs = {'download.default_directory' : self.paths['cwd'] + self.paths['sel_xls_download_path']}
        chrome_options.add_experimental_option('prefs', prefs)
        self.driver = webdriver.Chrome(chrome_options=chrome_options)
        
        self.driver.get("chrome://settings/security")
        

        
        root = self.driver.find_element_by_tag_name("settings-ui")
        shadow_root = self._expand_shadow_element(root)
        
        root1 = shadow_root.find_element_by_tag_name("settings-main")
        shadow_root1 = self._expand_shadow_element(root1)
        
        root2 = shadow_root1.find_element_by_tag_name("settings-basic-page")
        shadow_root2 = self._expand_shadow_element(root2)
        
        root3 = shadow_root2.find_element_by_tag_name("settings-privacy-page")
        shadow_root3 = self._expand_shadow_element(root3)
        
        root4 = shadow_root3.find_element_by_tag_name("settings-security-page")
        shadow_root4 = self._expand_shadow_element(root4)
        
        root5 = shadow_root4.find_element_by_css_selector("#safeBrowsingStandard")
        shadow_root5 = self._expand_shadow_element(root5)
        
        security_button = shadow_root5.find_element_by_css_selector("div.disc-border")
        
        security_button.click()
        """
        
        """
        fp = webdriver.FirefoxProfile("l2pnahxq.scraper")
        fp.set_preference("browser.download.folderList",2)
        fp.set_preference("browser.download.dir", self.paths['cwd'] + self.paths['sel_xls_download_path'])
        """
        fp = webdriver.FirefoxProfile("l2pnahxq.scraper")
        fp.set_preference("browser.download.folderList", 2)
        fp.set_preference("browser.download.manager.showWhenStarting", False)
        fp.set_preference("browser.download.dir", self.paths["sel_xls_download_path"])
        fp.set_preference("browser.helperApps.neverAsk.saveToDisk", "application/octet-stream")
        self.driver = webdriver.Firefox(firefox_profile=fp, executable_path="geckodriver.exe")         
        self.driver.get("http://sabiork.h-its.org/newSearch/index")
        
        # estimate the completion time
        minutes_per_enzyme = 9.5
        reactions_quantity = len(self.model['reactions'])
        estimated_time = reactions_quantity*minutes_per_enzyme*minute
        estimated_completion = datetime.datetime.now() + datetime.timedelta(seconds = estimated_time)
        
        print(f'Estimated completion of scraping data for {self.bigg_model_name}): {estimated_completion}, in {estimated_time/hour} hours')
        
        self.model_contents = {}
        for reaction in self.model["reactions"]:
            # parse the reaction
            annotations = reaction['annotation']
            reaction_id = reaction['id']
            reaction_name = reaction['name']
            og_reaction_string = self.bigg_reactions[reaction_id]['reaction_string']
            
            reaction_string, compounds = self._split_reaction(og_reaction_string)
            self.model_contents[reaction_name] = {
                'reaction': {
                    'original': og_reaction_string,
                    'substituted': reaction_string,
                },
                'chemicals': compounds,
                'annotations': annotations
            }
    
            # search SABIO for reaction kinetics
            if not reaction["name"] in self.variables['scraped_xls']:
                success_flag = False
                annotation_search_pairs = {"sabiork":"SabioReactionID", "metanetx.reaction":"MetaNetXReactionID", "ec-code":"ECNumber", "kegg.reaction":"KeggReactionID", "rhea":"RheaReactionID"}
                for database in annotation_search_pairs:
                    if not success_flag:
                        if database in annotations:
                            for ID in annotations[database]:
#                                 success_flag = self.scrape_xls(ID, annotation_search_pairs[database])
                                try:
                                    success_flag = self.scrape_xls(ID, annotation_search_pairs[database])
                                except:
                                    success_flag = False
                    else:
                        break

                if not success_flag:
                    try:
                        success_flag = self.scrape_xls(reaction_name, "Enzymename")
                    except:
                        success_flag = False

                json_dict_key = reaction_name.replace("\"", "")
                if success_flag:
                    self.variables['scraped_xls'][json_dict_key] = "yes"
                else:
                    self.variables['scraped_xls'][json_dict_key] = "no"
                    
                self.count += 1
                print("\nReaction: " + str(self.count) + "/" + str(len(self.model["reactions"])), end='\r')

            with open(self.paths['scraped_xls_path'], 'w') as outfile:
                json.dump(self.variables['scraped_xls'], outfile, indent = 4)   
                outfile.close()
                
        if self.export_content:
            with open(f'processed_{self.bigg_model_name}_model.json', 'w') as out:
                json.dump(self.model_contents, out, indent = 3)
                
        # update the step counter
        self.step_number = 2
        self._progress_update(self.step_number)

    """
    --------------------------------------------------------------------
        STEP 2: GLOB EXPORTED XLS FILES TOGETHER
    --------------------------------------------------------------------
    """

    def glob_xls_files(self,):
#         scraped_sans_parentheses_enzymes = glob('./{}/*.xls'.format(self.paths['xls_download_path']))
        total_dataframes = []
        for file in glob(os.path.join(self.paths['xls_download_path'], '*.xls')):
            #file_name = os.path.splitext(os.path.basename(file))[0]
            dfn = pandas.read_excel(file)
            total_dataframes.append(dfn)

        # All scraped dataframes are combined and duplicate rows are removed
        combined_df = pandas.DataFrame()
        combined_df = pandas.concat(total_dataframes)
        combined_df = combined_df.fillna(' ')
        combined_df = combined_df.drop_duplicates()

        # export the dataframe
        self.paths['concatenated_data'] = os.path.join(self.paths['sub_directory_path'], "proccessed-xls") + ".csv"
        combined_df.to_csv(self.paths['concatenated_data'])
        
        # update the step counter
        self.step_number = 3
        self._progress_update(self.step_number)

    """
    --------------------------------------------------------------------
        STEP 3: SCRAPE ADDITIONAL DATA BY ENTRYID
    --------------------------------------------------------------------    
    """

    def _scrape_entry_id(self,entry_id):
        entry_id = str(entry_id)

        fp = webdriver.FirefoxProfile("l2pnahxq.scraper")
        fp.set_preference("browser.download.folderList", 2)
        fp.set_preference("browser.download.manager.showWhenStarting", False)
        fp.set_preference("browser.download.dir", self.paths["sel_xls_download_path"])
        fp.set_preference("browser.helperApps.neverAsk.saveToDisk", "application/octet-stream")
        self.driver = webdriver.Firefox(firefox_profile=fp, executable_path="geckodriver.exe")         

        self.driver.get("http://sabiork.h-its.org/newSearch/index")
        time.sleep(self.parameters['general_delay']*2)

        self._click_element_id("option")
        self._select_dropdown_id("searchterms", "EntryID")
        text_area = self.driver.find_element_by_id("searchtermField")
        text_area.send_keys(entry_id)
        
        time.sleep(self.parameters['general_delay'])
        self._click_element_id("addsearch")
        time.sleep(self.parameters['general_delay'])
        self._click_element_id(entry_id + "img")
        time.sleep(self.parameters['general_delay'])
        
        self.driver.switch_to.frame(self.driver.find_element_by_xpath("//iframe[@name='iframe_" + entry_id + "']"))
        while True:
            try:
                element = self.driver.find_element_by_xpath("//table")                
                break
            except:
                time.sleep(self.parameters['general_delay'])

        element = self.driver.find_element_by_xpath("//table")
        html_source = element.get_attribute('innerHTML')
        table_df = pandas.read_html(html_source)
        reaction_parameters_df = pandas.DataFrame()
        counter = 0
        parameters_json = {}
        for df in table_df:
            try:
                if df[0][0] == "Parameter":
                    reaction_parameters_df = table_df[counter]
            except:
                self.driver.close()
                return parameters_json
            counter += 1
            
#        if self.printing:
#            print(reaction_parameters_df)
        for i in range(len(reaction_parameters_df[0])-2):
            print('i', i)
            parameter_name = reaction_parameters_df[0][i+2]
            inner_parameters_json = {}
            for j in range(len(reaction_parameters_df)-3):
                print('j', j)
                inner_parameters_json[str(reaction_parameters_df[j+1][1])] = reaction_parameters_df[j+1][i+2]

            parameters_json[parameter_name] = inner_parameters_json
            print(parameters_json)
        self.driver.close()
        return parameters_json

    def scrape_entryids(self,):
        if 'concatenated_data' not in self.paths:
            self.paths['concatenated_data'] = os.path.join(self.paths['sub_directory_path'], "proccessed-xls") + ".csv"
        self.sabio_xls_df = pandas.read_csv(self.paths['concatenated_data'])
        entryids = sabio_xls_df["EntryID"].unique().tolist()

        for entryid in entryids:
            if not entryid in self.variables['entryids_progress']:
#                 try:
                self.variables['entryids_progress'][str(entryid)] = self._scrape_entry_id(entryid)
                self.variables['scraped_entryids'][str(entryid)] = "yes"
#                 except:
#                     self.variables['scraped_entryids'][entryid] = "no"
            with open(self.paths["scraped_entryids_path"], 'w') as outfile:
                json.dump(self.variables['scraped_entryids'], outfile, indent = 4)   
                outfile.close()
            with open(self.paths["entryids_progress_path"], 'w') as f:
                json.dump(self.variables['entryids_progress'], f, indent = 4)        
                f.close()
        
        # update the step counter
        self.step_number = 4
        self._progress_update(self.step_number)

    """
    --------------------------------------------------------------------
        STEP 4: COMBINE ENZYME AND ENTRYID DATA INTO JSON FILE
    --------------------------------------------------------------------
    """   

    def combine_data(self,):
        # reviewing the progress of parsing the entry_ids
        with open(self.paths['entryids_progress_path']) as json_file: 
            entry_id_data = json.load(json_file)

        enzymenames = self.sabio_xls_df["Enzymename"].unique().tolist()
        enzyme_dict = {}
        missing_entry_ids = []
        parameters = {}

        for enzyme in enzymenames:
            sabio_grouped_enzyme_df = self.sabio_xls_df.loc[self.sabio_xls_df["Enzymename"] == enzyme]
            dict_to_append = {}
            reactions = sabio_grouped_enzyme_df["Reaction"].unique().tolist()
            for reaction in reactions:
                print(reaction)
                
                # ensure that the reaction chemicals match before accepting kinetic data
                rxn_string, compounds = self._split_reaction(reaction, sabio = True)
                bigg_compounds = self.model_content[enzyme]['chemicals']
                if set(compounds) != set(bigg_compounds):
                    print(compunds, bigg_compounds)
                    continue
                
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

                            stripped_string = re.sub('[0-9]', '', rate_law)

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
                                                if value != "-" and value != "" and value != " ":           # The quantities must be converted to base units
                                                    out_rate_law = out_rate_law.replace(var, parameter_info[var][permutation])
                                        except:
                                            pass

                            dict_entryid_to_append["RateLawSubstrates"] = substrates_key
                            dict_entryid_to_append["SubstitutedRateLaw"] = out_rate_law

                dict_to_append[reaction] = dict_reactions_to_append

            enzyme_dict[enzyme] = dict_to_append

        with open(self.paths["scraped_model_json_file_path"], 'w') as f:
            json.dump(enzyme_dict, f, indent=4)
        
        
    def _progress_update(self, step):
        if not re.search('[0-5]', str(step)):
            print(f'--> ERROR: The {step} step is not acceptable.')
        f = open(self.paths['progress_path'], "w")
        f.write(str(step))
        f.close()

    def main(self,
             bigg_model_path: dict,
             bigg_model_name: str = None
             ):
        self.start(bigg_model_path, bigg_model_name)

        while True:
            if self.step_number == 1:
                self.scrape_bigg_xls()
            elif self.step_number == 2:
                self.glob_xls_files()
            elif self.step_number == 3:
                self.scrape_entryids()
            elif self.step_number == 4:
                self.combine_data()
            elif self.step_number == 5:
                print("Execution complete. Scraper finished.")
                break
            
        # update the step counter
        self.step_number = 5
        self._progress_update(self.step_number)
            
            
scraping = SABIO_scraping()
scraping.main('Ecoli core, BiGG, indented.json', 'test_Ecoli')