# coding=utf-8
################################################################################
#
# Copyright (c) 2019 Baidu.com, Inc. All Rights Reserved
#
################################################################################
"""Crawl CoronaVirus DNA."""
import os
import time
import click
from selenium.webdriver import Chrome, ChromeOptions
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

LOGIN_URL = 'https://platform.gisaid.org/epi3/frontend#2245d'

def get_driver(download_dir):
	print("Initializing chrome driver...")
	options = ChromeOptions()
	prefs = {'download.default_directory' : download_dir}
	options.add_experimental_option('prefs', prefs)
	# options.add_argument('--headless')
	# options.add_argument('--disable-gpu') 
	main_driver = Chrome(options=options)

	print("Getting to GISAID website...")
	main_driver.get(LOGIN_URL)
	return main_driver


# Login with username and password
def login(main_driver, username, password):
	print("Logging in with username and password...")
	try:
		WebDriverWait(main_driver, 10).until(
			EC.presence_of_element_located((By.CLASS_NAME, 'form_button_submit'))
		)
		WebDriverWait(main_driver, 10).until(
			EC.presence_of_element_located((By.ID, 'elogin'))
		)
		main_driver.find_element_by_id('elogin').send_keys(username)
		main_driver.find_element_by_id('epassword').send_keys(password)
	except Exception as e:
		print(e)
	time.sleep(5)
	main_driver.find_elements_by_class_name('form_button_submit')[0].click()
	time.sleep(5)
	try:
		WebDriverWait(main_driver, 10).until(
			EC.presence_of_element_located((By.ID, 'main_nav'))
		)
	except Exception as e:
		print(e)


def go_to_EpiCov(main_driver):
	print("Going to EpiCov...")
	main = main_driver.find_element_by_id('main_nav')
	ul = main.find_elements_by_tag_name('ul')
	ul0 = ul[0]
	tags = ul0.find_elements_by_tag_name('li')
	beta_cov = tags[2]
	beta_cov.click()
	time.sleep(5)


def go_to_EpiCov_browser(main_driver):
	print("Going to EpiCov browser...")
	try:
		WebDriverWait(main_driver, 10).until(
			EC.presence_of_element_located((By.CLASS_NAME, 'sys-actionbar-action'))
		)
	except Exception as e:
		print(e)

	browse = main_driver.find_elements_by_class_name("sys-actionbar-action")[1]
	browse.click()
	time.sleep(5)


def sort_by_submission_date(main_driver):
	print("Sorting by submission date...")
	submission_date = main_driver.find_element_by_xpath("//a[@class='yui-dt-sortable'][text()='Submission Date']")
	submission_date.click()
	time.sleep(10)
	submission_date = main_driver.find_element_by_xpath("//a[@class='yui-dt-sortable'][text()='Submission Date']")
	submission_date.click()
	time.sleep(10)
	submission_date = main_driver.find_element_by_xpath("//a[@class='yui-dt-sortable'][text()='Submission Date']")
	assert submission_date.get_property("title") == "Click to sort ascending"


def crawl_metadata(main_driver, out_fn):
	table_str = ""
	table = main_driver.find_element_by_tag_name("table")
	for row in table.find_elements_by_tag_name("tr"):
		cells = row.find_elements_by_tag_name("td")
		row_str = "\t".join([cell.text.replace("\n", " ") for cell in cells])
		table_str += f"{row_str}\n"

	with open(out_fn, "w") as f:
		f.write(table_str)


def find_page(main_driver):
	global page_num

	print(f"Finding page {page_num}...")
	time.sleep(5)
	found = False
	current_page_num = int(main_driver.find_element_by_class_name("yui-pg-current-page").text)
	while int(current_page_num) < page_num:
		for page in main_driver.find_elements_by_class_name("yui-pg-page"):
			current_page_num = int(page.text)
			if int(page.text) == page_num:
				page.click()
				time.sleep(5)
				found = True
				break
				
		if not found:
			page.click()
			time.sleep(5)
	print(f"At page {current_page_num}")


def crawl(main_driver, out_dir):
	global page_num
	global finished
	print("Crawling webpage")
	while True:
		table = main_driver.find_elements_by_class_name('sys-datatable')[0]
		trs = table.find_elements_by_tag_name('tr')

		for tr in trs:
			if not tr.text.startswith('hCoV'):
				continue
			accession_id = tr.text.split("\n")[2]

			if os.path.exists(f"{out_dir}/{accession_id}.info"):
				print(f"{accession_id} already downloaded!")
				continue
				# finished = True
			else:
				print(f"Downloading {accession_id}...")

			tr.click()
			time.sleep(5)
			iframe = None
			try:
				WebDriverWait(main_driver, 10).until(
					EC.presence_of_element_located((By.TAG_NAME, 'iframe'))
				)
				iframe = main_driver.find_elements_by_tag_name('iframe')[0]
			except Exception as e:
				print(e)

			main_driver.switch_to.frame(iframe)
			back_btn = None
			crawl_metadata(main_driver, f"{out_dir}/{accession_id}.info")
			# for btn in main_driver.find_elements_by_class_name('sys-form-button'):
				# if 'Download Metadata' in btn.text:
				# 	print("Downloading Metadata...")
				# 	btn.click()
				# 	time.sleep(5)
				# if 'Download FASTA' in btn.text:
				# 	print("Downloading FASTA...")
				# 	btn.click()
				# 	time.sleep(5)
				# if 'Back' in btn.text:
					# back_btn = btn
			# if back_btn is None:
				# raise Exception('Back button missing!')
			back_btn = main_driver.find_element_by_xpath("//button[text()='Back']")
			back_btn.click()
			print("Going back...")
			time.sleep(5)

		next_btn = main_driver.find_elements_by_class_name('yui-pg-next')
		if len(next_btn) != 1:
			break
		next_btn = next_btn[0]
		next_btn.click()
		time.sleep(5)
		print('\n\n')
		print('page: ' + str(page_num))
		print('\n')
		page_num += 1


@click.command()
@click.option("-u", "--username", type=str, required=True, \
	help="Username for GISAID website.")
@click.option("-p", "--password", type=str, required=True, \
	help="Password for GISAID website.")
@click.option("-o", "--out_dir", type=str, required=True, \
	help="Output directory to place the downloaded sequences.")
def main(out_dir, username, password):
	os.makedirs(out_dir, exist_ok=True)
	
	try:
		main_driver = get_driver(out_dir)
		login(main_driver, username, password)
		go_to_EpiCov(main_driver)
		go_to_EpiCov_browser(main_driver)
		sort_by_submission_date(main_driver)
		find_page(main_driver)
		crawl(main_driver, out_dir)
	except Exception as e:
		main_driver.close()
		print(e)
		print(f"Error. Resume on page {page_num}")


if __name__ == "__main__":
	page_num = 1
	finished = False
	username = "lbxjollier" 
	password = "71RwYNz4nljy" 
	out_dir = "/Users/boxiang/Documents/work/Baidu/projects/viraviz/data/gisaid/metadata/"
	while not finished:
		# main(out_dir, username, password)
		main()
