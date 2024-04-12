import requests
import bs4 as bs
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
import time
from openai import OpenAI
import pandas as pd
# from webdriver_manager.chrome import ChromeDriverManager

# browser = webdriver.Chrome(ChromeDriverManager.install())
browser = webdriver.Chrome()
client = OpenAI()


def check_robots(root):
    page = requests.get(root + "/robots.txt").text
    if "Crawl-delay" in page:
        split_page = page.split("\n")
        for rule in split_page:
            if "Crawl-delay" in rule:
                return rule[-1]
    return 0

def get_links(links, blacklist, processedLinks):
    possibleTags = set()
    for link in links:
        if (
            link.text not in blacklist
            and len(link.text.split(" ")) <= 8
            and link.text.replace("\n", " ") not in processedLinks
            and 'mailto' not in link.get_attribute('outerHTML')
        ):
            possibleTags.add(Clickable(link, link.text.replace("\n", " "), link.get_attribute("href")))
    return possibleTags

def process_tags(
    system_message, prompt_example, assistant_example, blacklist, processedLinks, processedButtons
):
    possibleTags = set()
    buttons = browser.find_elements(By.TAG_NAME, "button")
    links = browser.find_elements(By.TAG_NAME, "a")

    for button in buttons:
        if (
            button.text not in blacklist
            and len(button.text.split(" ")) <= 8
            and button.text not in processedButtons
            and 'submit' not in button.get_attribute('outerHTML')
        ):
            processedButtons.add(button.text)
            button.send_keys(Keys.ENTER)
            time.sleep(.5)
            possibleTags = possibleTags.union(get_links(browser.find_elements(By.TAG_NAME, "a"), blacklist, processedLinks))
    
    #catch any stragglers
    possibleTags = possibleTags.union(get_links(links, blacklist, processedLinks))
    
    #create prompt
    prompt = "---".join(tag.getText() for tag in possibleTags)

    response = client.chat.completions.create(
        model="gpt-4-1106-preview",
        messages=[
            {"role": "system", "content": system_message},
            {"role": "user", "content": prompt_example},
            {"role": "assistant", "content": assistant_example},
            {"role": "user", "content": prompt},
        ],
    )

    ai_tags = response.choices[0].message.content
    # print(prompt)
    # print(ai_tags + "\n")
    if "Done" in ai_tags:
        return processedButtons, "Done"
    if ai_tags != "":
        ai_tags_list = ai_tags.split("---")
        final_list = []
        for tag in possibleTags:
            if tag.getText() in ai_tags_list:
                final_list.append(tag)

        return processedButtons, final_list
    else:
        return processedButtons, None


class Clickable:

    """
    tag: selenium tag
    text: text within selenium tag
    url: page where the button/link is displayed (url)
    """

    def __init__(self, tag, text, url):
        self.tag = tag
        self.text = text
        self.url = url

    def getTag(self):
        return self.tag

    def getText(self):
        return self.text

    def getUrl(self):
        return self.url


def check_if_name(tag):
    if tag.getText() != None:
        response = client.chat.completions.create(
            model="gpt-3.5-turbo-1106",
            messages=[
                {
                    "role": "system",
                    "content": "Check if this is possibly a person's name. If it is, output 'True'. Otherwise, output 'False'. Do not provide any extra explanation",
                },
                {"role": "user", "content": "Jackson Trials, Ph.D."},
                {"role": "assistant", "content": "True"},
                {"role": "user", "content": tag.getText()},
            ],
        )
        return response.choices[0].message.content
    return "Empty Tag"

def scrape_profile(text):
    response = client.chat.completions.create(
        model="gpt-3.5-turbo-1106",
        messages=[
            {
                "role": "system",
                "content": "You are a helpful assistant. The following is likely a researcher's profile. Extract the following pieces of info in this order: researcher's name, email, address, phone number, and a very brief summary of their research. Use the following format: Separate each piece of info with '---'. If you cannot find a piece of info, output 'N' where the info would be.",
            },
            {"role": "user", "content": 
                """Jennifer A. Kanakry, M.D.
                National Cancer Institute
                Building 10-CRC, Room 4-3132
                Bethesda, MD 20892
                240-760-6172
                jennifer.kanakry@nih.gov(link sends email)
                DIRECTOR
                TEAM MEMBER OF:
                Alex Compton, Ph.D.
                INVESTIGATOR
                alex.compton@nih.gov
                RESEARCH SUMMARY
                Dr. Jennifer Kanakry engages in clinical research related to allogeneic hematopoietic cell transplantation for patients with hematologic malignancies and inherited and acquired disorders of the immune system.
                """
             },
            {"role": "assistant", "content": "Jennifer A. Kanakry---jennifer.kanakry@nih.gov---National Cancer Institute, Building 10-CRC, Room 4-3132, Bethesda, MD 20892---240-760-6172---Researches cell transplantation for patients with blood cancers and immune system disorders."},
            {"role": "user", "content": text},
        ],
    )
    return response.choices[0].message.content.split("---")

def search_url(root):
    timeDelay = check_robots(root)

    candidates = [Clickable(None, None, root)] #stack of links wrapped in Clickable class
    with open("blacklist.txt", "r", encoding="utf-8") as f:
        blacklist = f.read().split("\n")
    blacklist = set(blacklist)
    processedLinks = set() #set of processed urls
    processedButtons = set()
    data = {
        'Name': [],
        'Email': [],
        'Address': [],
        'Phone': [],
        'Summary':[]
    }

    while candidates:
        if len(data.get('Name')) > 100:
            break
        #remove from top of stack (depth first search)
        next = candidates.pop()
        url = next.getUrl()
        
        if url not in processedLinks:
            browser.get(url)
            processedLinks.add(url)
            time.sleep(timeDelay)

            #If we've landed a profile!!
            if check_if_name(next) == "True":
                body = browser.find_element(By.TAG_NAME, "body")
                profile = scrape_profile(body.text)
                for count, key in enumerate(data.keys()):
                    data[key].append(profile[count])
                
            else:
                # get a list of relevant tag objects
                processedButtons, tags = process_tags(
                    system_message="""Your role as a data collection assistant is to meticulously examine a website for pages that specifically display staff directories with links to faculty profiles. Upon receiving a list of clickable links, each distinctly separated by '---', you must discern and select only those links that directly correspond to your objective. Preserve the '---' format in your output. Please optimize for directories containing all information. Provide a maximum of two most probable links without additional explanation. If the majority of links refer to individual names, provide all links containing individual names. If all links are irrelevant, simply respond with 'None'""",
                    prompt_example="HOME---Check---Dr. W. Kimryn Rathmell---LinkedIn---About Us",
                    assistant_example="Dr. W. Kimryn Rathmell---Organization---About Us",
                    blacklist=blacklist,
                    processedLinks=processedLinks,
                    processedButtons=processedButtons
                )
                if tags == 'Done':
                    break
                #If this page isn't a dead end, keep the loop going
                if tags != None:
                    candidates += tags
    return data

if __name__ == "__main__":
    data = search_url("https://www.cancer.gov/")    
    df = pd.DataFrame(data)
    df.to_csv("data1.csv")
    # data = search_url("https://www.nei.nih.gov/")    
    # df = pd.DataFrame(data)
    # df.to_csv("data2.csv")
    # data = search_url("https://www.nhlbi.nih.gov/")    
    # df = pd.DataFrame(data)
    # df.to_csv("data3.csv")

    browser.quit()
    
    
    #If none of the presented links explicitly align with these criteria, it is highly preferable not to list any rather than include possibly irrelevant options.