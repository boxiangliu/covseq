import os
import time
import click
from selenium.webdriver import Chrome, ChromeOptions
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import sys
sys.path.append("download")
from download_gisaid import login, go_to_EpiCov, go_to_EpiCov_browser

LOGIN_URL = "https://www.epicov.org/epi3/frontend#5f0352"
CREDENTIALS_FN = "download/credentials.txt"

def read_credentials(CREDENTIALS_FN):
	with open(CREDENTIALS_FN) as f:
		username = f.readline().split("=")[1].strip()
		password = f.readline().split("=")[1].strip()
	return username, password

def select_all(driver):
	label = driver.find_elements_by_class_name("yui-dt-label")[0]
	input_ = label.find_element_by_tag_name("input")
	input_.click()
	time.sleep(10)

def download_sequence(driver):
	button = driver.find_elements_by_class_name("sys-form-button")[-1]
	assert "Download" in button.text
	button.click()
	time.sleep(5)
	iframe = driver.find_elements_by_tag_name('iframe')[0]
	driver.switch_to.frame(iframe)
	button = driver.find_elements_by_class_name("sys-event-hook")[-1]
	assert "Download" == button.text
	button.click()
	time.sleep(10)

username, password = read_credentials(CREDENTIALS_FN)
driver = Chrome(ChromeDriverManager().install())
driver.get(LOGIN_URL)
login(driver, username, password)
go_to_EpiCov(driver)
go_to_EpiCov_browser(driver)
select_all(driver)