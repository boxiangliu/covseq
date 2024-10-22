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
from download_gisaid import login as gisaid_login, go_to_EpiCov, go_to_EpiCov_browser
import glob
import subprocess
import json

GISAID_URL = "https://www.epicov.org/epi3/frontend#5f0352"
NCBI_URL = "https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=SARS-CoV-2,%20taxid:2697049"
# EMBL_URL = "https://www.ebi.ac.uk/ena/pathogens/covid-19"
EMBL_URL = "ftp://ftp.ebi.ac.uk/pub/databases/covid19dataportal/viral_sequences/sequences/"
CREDENTIALS_FN = "download/credentials.txt"
CONFIG_FN = "download/download_config.json"


def read_config(config_fn):
    with open(config_fn) as f:
        config = json.load(f)
    return config


def read_credentials(CREDENTIALS_FN):
    with open(CREDENTIALS_FN) as f:
        username = f.readline().split("=")[1].strip()
        password = f.readline().split("=")[1].strip()
    return username, password


# def select_all(driver):
# 	label = driver.find_elements_by_class_name("yui-dt-label")[0]
# 	input_ = label.find_element_by_tag_name("input")
# 	input_.click()
# 	time.sleep(10)


# def download_sequence(driver):
# 	button = driver.find_elements_by_class_name("sys-form-button")[-1]
# 	assert "Download" in button.text
# 	button.click()
# 	time.sleep(5)
# 	iframe = driver.find_elements_by_tag_name('iframe')[0]
# 	driver.switch_to.frame(iframe)
# 	button = driver.find_elements_by_class_name("sys-event-hook")[-1]
# 	assert "Download" == button.text
# 	button.click()
# 	time.sleep(10)


def go_to_EpiCov_downloads(main_driver):
    print("Going to EpiCov browser...")
    try:
        WebDriverWait(main_driver, 10).until(
            EC.presence_of_element_located(
                (By.CLASS_NAME, 'sys-actionbar-action'))
        )
    except Exception as e:
        print(e)

    browse = main_driver.find_elements_by_class_name("sys-actionbar-action")[2]
    browse.click()
    time.sleep(5)


def gisaid_download(driver):
    iframe = None
    try:
        WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.TAG_NAME, 'iframe'))
        )
        iframe = driver.find_elements_by_tag_name('iframe')[0]
    except Exception as e:
        print(e)

    driver.switch_to.frame(iframe)

    downicons = driver.find_elements_by_class_name("downicon")
    for downicon in downicons:
        if "nextmeta" in downicon.text:
            downicon.click()
            time.sleep(10)
        elif "nextfasta" in downicon.text:
            downicon.click()
            time.sleep(120)


def gisaid_decompress(config):
    while glob.glob(f"{config['download_dir']}/{config['gisaid_fasta_prefix']}*") == []:
        time.sleep(5)

    compressed_fasta_fn = glob.glob(f"{config['download_dir']}/{config['gisaid_fasta_prefix']}*")[0]
    cmd = f"gunzip {compressed_fasta_fn}"
    output = subprocess.run(cmd, shell=True)

    compressed_metadata_fn = glob.glob(f"{config['download_dir']}/{config['gisaid_metadata_prefix']}*")[0]
    cmd = f"gunzip {compressed_metadata_fn}"
    output = subprocess.run(cmd, shell=True)


def gisaid_move(config):
    src = glob.glob(f"{config['download_dir']}/{config['gisaid_fasta_prefix']}*")[0]
    os.rename(src, config["gisaid_target_fasta"])

    src = glob.glob(f"{config['download_dir']}/{config['gisaid_metadata_prefix']}*")[0]
    os.rename(src, config["gisaid_target_metadata"])


def ncbi_click_download(driver):
    for radio_selection in ["Nucleotide", "CSV format"]:
        download_btn = driver.find_element_by_class_name(
            "ncbi-report-download")
        download_btn.click()
        time.sleep(1)

        radio_btns = driver.find_elements_by_class_name("custom-control")
        for radio_btn in radio_btns:
            if radio_btn.text == radio_selection:
                radio_btn.click()
                time.sleep(1)
                break

        download_btn = driver.find_element_by_class_name("ncbi-download-btn")
        download_btn.click()
        time.sleep(1)

        action_btns = driver.find_elements_by_class_name("ncbi-download-btn")
        for btn in action_btns:
            if btn.text == "Next":
                btn.click()
                time.sleep(1)
                break

        if radio_selection == "CSV format":
            checkboxes = driver.find_elements_by_class_name(
                "ncbi-checkbox-label")
            for checkbox in checkboxes:
                if checkbox.text == "Select All":
                    checkbox.click()
                    break

        action_btns = driver.find_elements_by_class_name("ncbi-download-btn")
        for btn in action_btns:
            if btn.text == "Download":
                btn.click()
                time.sleep(1)
                break


def ncbi_move(config):
    src = f"{config['download_dir']}/{config['ncbi_fasta_fn']}"
    tgt = config["ncbi_target_fasta"]
    while not os.path.exists(src):
        time.sleep(5)
    os.rename(src, tgt)

    src = f"{config['download_dir']}/{config['ncbi_metadata_fn']}"
    tgt = config["ncbi_target_metadata"]
    os.rename(src, tgt)


def embl_download(EMBL_URL, config):
    cmd = f"wget --no-host-directories --cut-dirs=5 --no-parent --directory-prefix={config['download_dir']} -r {EMBL_URL}"
    subprocess.run(cmd, shell=True)


def twice_unzip(src):
    cmd = f"file {src}".replace(".gz","")
    output = subprocess.run(cmd, shell=True, capture_output=True)
    if b"gzip" in output.stdout:
        os.rename(src.replace(".gz",""), src)
        cmd = f"gunzip {src}"
        subprocess.run(cmd, shell=True)


def embl_unzip(config):
    src = glob.glob(f"{config['download_dir']}/{config['embl_fasta_prefix']}*.txt.gz")[-1]
    cmd = f"gunzip {src}"
    subprocess.run(cmd, shell=True)
    twice_unzip(src)

    src = glob.glob(f"{config['download_dir']}/{config['embl_metadata_prefix']}*.txt.gz")[-1]
    cmd = f"gunzip {src}"
    subprocess.run(cmd, shell=True)
    twice_unzip(src)


def embl_move(config):
    src = glob.glob(f"{config['download_dir']}/{config['embl_fasta_prefix']}*.txt")[0]
    os.rename(src, config["embl_target_fasta"])

    src = glob.glob(f"{config['download_dir']}/{config['embl_metadata_prefix']}*.txt")[0]
    os.rename(src, config["embl_target_metadata"])


config = read_config(CONFIG_FN)

username, password = read_credentials(CREDENTIALS_FN)
driver = Chrome(ChromeDriverManager().install())
# driver.get(GISAID_URL)
# gisaid_login(driver, username, password)
# go_to_EpiCov(driver)
# go_to_EpiCov_downloads(driver)
# gisaid_download(driver)
# gisaid_decompress(config)
# gisaid_move(config)


# driver.get(NCBI_URL)
# ncbi_click_download(driver)
# ncbi_move(config)


embl_download(EMBL_URL, config)
embl_unzip(config)
embl_move(config)


driver.close()
