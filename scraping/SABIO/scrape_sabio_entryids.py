from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import Select
import time
import json
import threading
import numpy
from itertools import islice
import pandas as pd

sabio_reactions_df = pd.read_csv("sabio_out.csv")
entryids = sabio_reactions_df["EntryID"].unique().tolist()


general_delay = 2

bad_entry_ids = []


def click_element_id(n_id):
    global driver
    element = driver.find_element_by_id(n_id)
    element.click()
    time.sleep(general_delay)
    
def select_dropdown_id(n_id, n_choice):
    global driver
    element = Select(driver.find_element_by_id(n_id))
    element.select_by_visible_text(n_choice)
    time.sleep(general_delay)

def chunks(data, SIZE):
    it = iter(data)
    for i in range(0, len(data), SIZE):
        yield {k:data[k] for k in islice(it, SIZE)}


final_json_out = {}

def scrape_entry_id(entry_id):
    global driver
    global final_json_out
    
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
        
    print(parameters_json)
    
    driver.close()
    
    return parameters_json




for entryid in entryids:
    try:
        final_json_out[str(entryid)] = scrape_entry_id(entryid)
    except:
        bad_entry_ids.append(str(entryid))

with open('data.json', 'w') as f:
    json.dump(final_json_out, f)

with open('bad_entry_ids.txt', 'w') as f:
    for entry in bad_entry_ids:
        f.write(entry)
        f.write(",")
    f.close()
    

    

