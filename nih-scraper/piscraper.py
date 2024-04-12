import requests
import bs4 as bs
import re
import tqdm
import webbrowser
import os
import tkinter as tk


def return_all_links(visited_set, blacklist, html):
    links = []
    soup = bs.BeautifulSoup(html, "html.parser")
    last_marked = soup.find_all("ul", {"class": "jumplinks__list"})[-1]

    for tag in last_marked.find_all_next("a", {"class": "pilist__link"}):
        links.append(tag.get("href"))

    return links


def run_search(searched_word="sequenc"):
    visited_set = set()
    root = "https://irp.nih.gov/"
    blacklist = ["jumplinks__link"]  # by class
    total_html = requests.get(
        "https://irp.nih.gov/our-research/principal-investigators/name"
    ).text
    result_list = []
    count = 0

    for link in tqdm.tqdm(
        return_all_links(visited_set, blacklist, total_html), position=0, leave=True
    ):
        if count > 100:
            break
        html = requests.get(root + link).text
        soup = bs.BeautifulSoup(html, "html.parser")

        results = soup.body.find_all(
            string=re.compile(".*{0}.*".format(searched_word)), recursive=True
        )
        if len(results) > 0:
            result_list.append([root + link, len(results)])

        count += 1

    return result_list

def handle_enter(event):
    print("pass")


if __name__ == "__main__":
    window = tk.Tk()
    frame = tk.Frame()
    title = tk.Label(
        text="NIH Private Investigator Scraper",
    )
    
    button = tk.Button(
        text="Start",
        width=25,
        height=5,
        bg="black",
        fg="white",
    )
    
    instructions = tk.Label(
        text="Enter a keyword to search...",
    )
    
    entry = tk.Entry()

    search_term = input("What would you like to search for?: ")
    print("Running search...")
    results = run_search(search_term)

    with open("template.html", "r", encoding="utf-8") as file:
        html = file.read()

    soup = bs.BeautifulSoup(html, "html.parser")
    print("Finalizing results...")

    title = soup.new_tag("h1", id="title")
    title.append(f"Search results for '{search_term}'")
    soup.body.insert(0, title)

    for result in tqdm.tqdm(results, position=0, leave=True):
        link_column = soup.find("tr", {"class": "link"})
        link_tag = soup.new_tag("td")
        link_nested_tag = soup.new_tag("a", href=result[0])
        link_nested_tag.append(result[0])
        link_tag.append(link_nested_tag)
        link_column.append(link_tag)

        num_column = soup.find("tr", {"class": "number"})
        num_tag = soup.new_tag("td")
        num_tag.append(str(result[1]))
        num_column.append(num_tag)

    with open("output.html", "w", encoding="utf-8") as file:
        file.write(str(soup))

    webbrowser.open_new_tab("output.html")


"""Getting all institutes"""
# browser.get('https://www.nih.gov/institutes-nih/list-institutes-centers')

# institutes = browser.find_elements(By.TAG_NAME,"p")

# institute_links = []
# institute_names = []

# for web_element in institutes:
#     try:
#         #because nih is stupid
#         if ("NHGRI" in web_element.get_attribute("innerText")):
#             name = web_element.find_element(By.TAG_NAME, "strong")
#             link = web_element.find_element(By.TAG_NAME, "a")
#         else:
#             link = web_element.find_element(By.TAG_NAME, "a")
#             name = link.find_element(By.TAG_NAME, "strong")
#         institute_links.append(link.get_attribute("href"))
#         institute_names.append(name.get_attribute("innerText"))
#     except Exception:
#         pass
