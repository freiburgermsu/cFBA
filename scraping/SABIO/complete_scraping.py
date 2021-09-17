# -*- coding: utf-8 -*-
"""
@authors: Ethan Chan, Matthew Freiburger
"""

# Import libraries
import pandas as pd
import glob
import os
from selenium import webdriver
from selenium.webdriver.support.ui import Select
import time
import json
from itertools import islice
import re
from os import path


class scraping():
    def __init__(self):
            
        general_delay = 2
        xls_download_directory = "sabio_downloaded_xls"
        progress_file_prefix = "current-progress-"
        xls_download_prefix = "xls-download-"
        scraped_xls_prefix = "scraped-xls-"
        scraped_entryids_prefix = "scraped-entryids-"
        sel_xls_download_path = ""

        bigg_model_name_suffix = ""
        sub_directory_path = ""
        progress_file_path = ""
        xls_download_path = ""
        scraped_xls_file_path = ""
        scraped_entryids_file_path = ""
        xls_csv_file_path = ""
        entryids_json_file_path = ""
        scraped_model_json_file_path = ""
        bigg_model = {}
        step_number = -1
        cwd = ""


        scraped_xls = {}
        scraped_entryids = {}

    #Clicks a HTML element with selenium by id
    def click_element_id(n_id):
        global driver
        element = driver.find_element_by_id(n_id)
        element.click()
        time.sleep(general_delay)

    #Selects a choice from a HTML dropdown element with selenium by id
    def select_dropdown_id(n_id, n_choice):
        global driver
        element = Select(driver.find_element_by_id(n_id))
        element.select_by_visible_text(n_choice)
        time.sleep(general_delay)

    def fatal_error_handler(message):
        print("Error: " + message)
        print("Exiting now...")
        exit(0)


    """
    --------------------------------------------------------------------
        STEP 0: GET BIGG MODEL TO SCRAPE AND SETUP DIRECTORIES AND PROGRESS FILE
    --------------------------------------------------------------------    
    """

    def start():
        global bigg_model
        global cwd
        global sub_directory_path
        global progress_file_path
        global xls_download_path
        global sel_xls_download_path 
        global scraped_xls_file_path
        global scraped_xls
        global scraped_entryids_file_path
        global scraped_entryids
        global xls_csv_file_path
        global entryids_json_file_path
        global entry_id_json_out
        global scraped_model_json_file_path
        global step_number

        cwd = os.path.dirname(os.path.realpath(__file__))

        #Get name of model to scrape
        while True:
            bigg_model_path = input("Please specify path of BIGG Model JSON file to start scraping: ")

            if path.exists(bigg_model_path) and bigg_model_path[-5:] == ".json":

                try:
                    bigg_model = json.load(open(bigg_model_path))
                    bigg_model_name_suffix = re.search("[^\\/:*?\"<>|\r\n]+$", bigg_model_path).group()[:-5]
                    break
                except:
                    pass

        sub_directory_path = "scraping-"+bigg_model_name_suffix

        if not os.path.isdir(sub_directory_path):        
            os.mkdir(sub_directory_path)

        progress_file_path = sub_directory_path + "/" + progress_file_prefix+bigg_model_name_suffix+".txt"
        xls_download_path = sub_directory_path + "/" + xls_download_prefix+bigg_model_name_suffix
        scraped_xls_file_path = sub_directory_path + "/" + "scraped-xls-" + bigg_model_name_suffix + ".json"
        scraped_entryids_file_path = sub_directory_path + "/" + "scraped-entryids-" + bigg_model_name_suffix + ".json"
        xls_csv_file_path = sub_directory_path + "/" + "proccessed-xls-" + bigg_model_name_suffix + ".csv"
        entryids_json_file_path = sub_directory_path + "/" + "entryids-json-" + bigg_model_name_suffix + ".json"
        scraped_model_json_file_path = sub_directory_path + "/" + "scraped-model-" + bigg_model_name_suffix + ".json"
        sel_xls_download_path = "\\" + sub_directory_path + "\\" + xls_download_prefix+bigg_model_name_suffix

        if os.path.exists(progress_file_path):
            f = open(progress_file_path, "r")

            step_number = int(f.read(1))

            if not step_number in [1, 2, 3, 4, 5]:
                fatal_error_handler("Progress file malformed. Please delete and restart")

        else:
            f = open(progress_file_path, "w")
            f.write("1")
            f.close()
            step_number = 1

        if not os.path.isdir(xls_download_path):
            os.mkdir(xls_download_path)

        if os.path.exists(scraped_xls_file_path):
            f = open(scraped_xls_file_path, "r")
            scraped_xls = json.load(f)
            f.close()

        if os.path.exists(scraped_entryids_file_path):
            f = open(scraped_entryids_file_path, "r")
            scraped_entryids = json.load(f)
            f.close()

        if os.path.exists(entryids_json_file_path):
            f = open(entryids_json_file_path, "r")
            entry_id_json_out = json.load(f)
            f.close()



    """
    --------------------------------------------------------------------
        STEP 1: SCRAPE SABIO WEBSITE BY DOWNLOAD XLS FOR GIVEN REACTIONS IN BIGG MODEL
    --------------------------------------------------------------------    
    """

    def scrape_xls(reaction_identifier, search_option):
        global driver
        global enzymes
        global xls_download_path


        chrome_options = webdriver.ChromeOptions()
        prefs = {'download.default_directory' : cwd + sel_xls_download_path}
        chrome_options.add_experimental_option('prefs', prefs)
        driver = webdriver.Chrome(chrome_options=chrome_options)    
        driver.get("http://sabiork.h-its.org/newSearch/index")

        time.sleep(general_delay)





        click_element_id("option")
        select_dropdown_id("searchterms", search_option)
        text_area = driver.find_element_by_id("searchtermField")
        text_area.send_keys(reaction_identifier)  
        time.sleep(general_delay)  
        click_element_id("addsearch")

        time.sleep(general_delay)

        result_num = ""
        try: 
            result_num_ele = driver.find_element_by_id("numberofKinLaw")
            for char in result_num_ele.text:
                if re.search('[0-9]', char):
                    result_num = result_num + char

            result_num = int(result_num)
        except:
            driver.close()
            return False

        time.sleep(general_delay)

        select_dropdown_id("max", "100")
        element = Select(driver.find_element_by_id("max"))
        element.select_by_visible_text("100")

        time.sleep(general_delay)

        if result_num > 0 and result_num <= 100:
            click_element_id("allCheckbox")
            time.sleep(general_delay)
        elif result_num > 100:
            click_element_id("allCheckbox")
            for i in range(int(result_num/100)):
                element = driver.find_element_by_xpath("//*[@class = 'nextLink']")
                element.click()
                time.sleep(general_delay)
                click_element_id("allCheckbox")
                time.sleep(general_delay)
        else:
            driver.close()
            return False

        driver.get("http://sabiork.h-its.org/newSearch/spreadsheetExport")

        time.sleep(15)

        click_element_id("excelExport")

        time.sleep(5)

        driver.close()

        return True


    def scrape_bigg_xls():
        global scraped_xls_file_path
        global scraped_xls
        global step_number


        for reaction in bigg_model["reactions"]:
            if not reaction["name"] in scraped_xls:
                ids_to_try = reaction["annotation"]

                success_flag = False


                annotation_search_pairs = {"sabiork":"SabioReactionID", "metanetx.reaction":"MetaNetXReactionID", "ec-code":"ECNumber", "kegg.reaction":"KeggReactionID", "rhea":"RheaReactionID"}


                for annotation in annotation_search_pairs:
                    if not success_flag:
                        if annotation in ids_to_try:
                            for id_to_try in ids_to_try[annotation]:
                                try:
                                    success_flag = scrape_xls(id_to_try, annotation_search_pairs[annotation])
                                except:
                                    success_flag = False
                    else:
                        break

                if not success_flag:
                    try:
                        success_flag = scrape_xls(reaction["name"], "Enzymename")
                    except:
                        success_flag = False


                json_dict_key = reaction["name"].replace("\"", "")
                if success_flag:

                    scraped_xls[json_dict_key] = "yes"
                else:
                    scraped_xls[json_dict_key] = "no"

            with open(scraped_xls_file_path, 'w') as outfile:
                json.dump(scraped_xls, outfile, indent = 4)   
                outfile.close()
        f = open(progress_file_path, "w")
        f.write("2")
        f.close()
        step_number = 2


    """
    --------------------------------------------------------------------
        STEP 2: GLOB EXPORTED XLS FILES TOGETHER
    --------------------------------------------------------------------
    """

    def glob_xls_files():
        global step_number
        global xls_csv_file_path

        path_sans_parentheses = './' + xls_download_path + '/'

        sans_parentheses_files = glob.glob(os.path.join(path_sans_parentheses, '*.xls'))

        # create the total list of dataframes
        total_dataframes = []

        for file in sans_parentheses_files:
            #file_name = os.path.splitext(os.path.basename(file))[0]
            dfn = pd.read_excel(file)
            total_dataframes.append(dfn)

        # combine the total set of dataframes
        combined_df = pd.DataFrame()
        combined_df = pd.concat(total_dataframes)

        # replace the NaN values with blank spaces
        combined_df = combined_df.fillna(' ')

        combined_df = combined_df.drop_duplicates()

        # export the dataframe
        combined_df.to_csv(xls_csv_file_path)
        f = open(progress_file_path, "w")
        f.write("3")
        f.close()
        step_number = 3


    """
    --------------------------------------------------------------------
        STEP 3: SCRAPE ADDITIONAL DATA BY ENTRYID
    --------------------------------------------------------------------    
    """

    def scrape_entry_id(entry_id):
        global driver
        global entry_id_json_out

        entry_id = str(entry_id)

        driver = webdriver.Chrome(executable_path=r".\chromedriver.exe")

        driver.get("http://sabiork.h-its.org/newSearch/index")

        time.sleep(general_delay)

        click_element_id("option")

        select_dropdown_id("searchterms", "EntryID")

        text_area = driver.find_element_by_id("searchtermField")
        text_area.send_keys(entry_id)

        time.sleep(general_delay)

        click_element_id("addsearch")

        time.sleep(general_delay)

        click_element_id(entry_id + "img")

        time.sleep(general_delay)

        driver.switch_to.frame(driver.find_element_by_xpath("//iframe[@name='iframe_" + entry_id + "']"))

        element = driver.find_element_by_xpath("//table")
        html_source = element.get_attribute('innerHTML')

        table_df = pd.read_html(html_source)


        reaction_parameters_df = pd.DataFrame()

        counter = 0

        parameters_json = {}

        for df in table_df:
            try:
                if df[0][0] == "Parameter":
                    reaction_parameters_df = table_df[counter]
            except:
                 driver.close()
                 return parameters_json
            counter += 1


        parameter_name = ""

        for i in range(len(reaction_parameters_df[0])-2):
            parameter_name = reaction_parameters_df[0][i+2]


            inner_parameters_json = {}
            for j in range(len(reaction_parameters_df)-3):
                inner_parameters_json[reaction_parameters_df[j+1][1]] = reaction_parameters_df[j+1][i+2]


            parameters_json[parameter_name] = inner_parameters_json

        driver.close()

        return parameters_json


    def scrape_entryids():
        global step_number
        global entry_id_json_out

        sabio_xls_df = pd.read_csv(xls_csv_file_path)
        entryids = sabio_xls_df["EntryID"].unique().tolist()

        for entryid in entryids:
            if not entryid in scraped_entryids:
                try:
                    entry_id_json_out[str(entryid)] = scrape_entry_id(entryid)
                    scraped_entryids[entryid] = "yes"
                except:
                    scraped_entryids[entryid] = "no"
            with open(scraped_entryids_file_path, 'w') as outfile:
                json.dump(scraped_entryids, outfile, indent = 4)   
                outfile.close()
            with open(entryids_json_file_path, 'w') as f:
                json.dump(entry_id_json_out, f, indent = 4)        
                f.close()
        f = open(progress_file_path, "w")
        f.write("4")
        f.close()
        step_number = 4


    """
    --------------------------------------------------------------------
        STEP 4: COMBINE ENZYME AND ENTRYID DATA INTO JSON FILE
    --------------------------------------------------------------------
    """

    def combine_data():
        global scraped_model_json_file_path
        global step_number

        sabio_xls_df = pd.read_csv(xls_csv_file_path)

        # Opening JSON file
        with open(entryids_json_file_path) as json_file:
            entry_id_data = json.load(json_file)

        enzymenames = sabio_xls_df["Enzymename"].unique().tolist()

        enzyme_dict = {}

        missing_entry_ids = []

        parameters = {}

        for enzyme in enzymenames:
            sabio_grouped_enzyme_df = sabio_xls_df.loc[sabio_xls_df["Enzymename"] == enzyme]


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




                dict_to_append[reaction] = dict_reactions_to_append

            enzyme_dict[enzyme] = dict_to_append

        with open(scraped_model_json_file_path, 'w') as f:
            json.dump(enzyme_dict, f, indent=4)
        f = open(progress_file_path, "w")
        f.write("5")
        f.close()
        step_number = 5

    def main():
        global step_number
        start()

        while True:
            if step_number == 1:
                scrape_bigg_xls()
            elif step_number == 2:
                glob_xls_files()
            elif step_number == 3:
                scrape_entryids()
            elif step_number == 4:
                combine_data()
            elif step_number == 5:
                print("Execution complete. Scraper finished.")
                break

    main()
