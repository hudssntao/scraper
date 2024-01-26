from selenium import webdriver
from selenium.webdriver.common.by import By

browser = webdriver.Chrome()
browser.get("https://ccr.cancer.gov/staff-directory/sherimay-ablan")

tag = browser.find_element(By.TAG_NAME, "body")
print(tag.text)


print(browser.current_url)
