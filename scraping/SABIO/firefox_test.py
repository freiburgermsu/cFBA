from selenium import webdriver
import os
fp = webdriver.FirefoxProfile("C:/Users/Ethan/Documents/Python/Polar Bar Chart/sabioawesome/l2pnahxq.scraper")
fp.set_preference("browser.download.folderList",2)
fp.set_preference("browser.download.dir", os.getcwd())
driver = webdriver.Firefox(firefox_profile=fp, executable_path="geckodriver.exe")