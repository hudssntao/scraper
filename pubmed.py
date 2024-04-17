from Bio import Entrez
import os
import json
import pandas as pd
from xml.dom import minidom    
import tqdm
from dotenv import load_dotenv
import argparse

load_dotenv()

user = {
    'Email': '',
}


data = {
        'Name': [],
        'Email': [],
        'Address': [],
        'Keyword': [],
        'Title': [],
        'Date': []
    }

#for papers outside of the U.S.
non_US = {
        'Name': [],
        'Email': [],
        'Address': [],
        'Keyword': [],
        'Title': [],
        'Date': []
    }

#extract text from an xml node (SHALLOW, does not go beyond one layer)
def getText(nodes):
    text = []
    for node in nodes:
        if node.nodeType == node.TEXT_NODE:
            text.append(node.data)
    return ''.join(text)

#recursively extracts all text starting from the root
def getAllText(root):
    text = []
    for childNode in root.childNodes:
        if childNode.nodeType == childNode.TEXT_NODE:
            text.append(childNode.data)
        else:
            text.append(getAllText(childNode))
    return ' '.join(text)

def parseAffiliation(affiliation_node):
    affiliation_text = getText(affiliation_node[0].childNodes)
    count = 1
    for word in reversed(affiliation_text.split(" ")):
        if "@" in word:
            address = " ".join(affiliation_text.split(" ")[:(count * -1)]).replace('Electronic address:', '')
            if word[-1] == ".":
                return address, word[:-1]
            else:
                return address, word
        count += 1

def searchKeyword(keyword, num_papers):
    print(f"Processing {keyword}...\n")
    
    # Perform the PubMed search
    handle = Entrez.esearch(db="pubmed", term=keyword, sort="relevance", retmax=num_papers)
    record = Entrez.read(handle)
    handle.close()

    # Get the list of PubMed IDs (PMID) for the search results
    pmid_list = record["IdList"]

    data = {
            'Name': [],
            'Email': [],
            'Address': [],
            'Keyword': [],
            'Title': [],
            'Date': []
        }

    #for papers outside of the U.S.
    non_US = {
            'Name': [],
            'Email': [],
            'Address': [],
            'Keyword': [],
            'Title': [],
            'Date': []
        }

    # Fetch information for each article
    for pmid in tqdm.tqdm(pmid_list, leave=True, position=0):
        try:
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="xml")
            article_info = handle.read()
            handle.close()
            
            xmlDoc = minidom.parseString(article_info)
            title = getText(xmlDoc.getElementsByTagName("ArticleTitle")[0].childNodes)
            abstract = getText(xmlDoc.getElementsByTagName("AbstractText")[0].childNodes) if xmlDoc.getElementsByTagName("AbstractText") else ""
            if xmlDoc.getElementsByTagName("ArticleDate"):
                date = getAllText(xmlDoc.getElementsByTagName("ArticleDate")[0]).split(" ") 
                year = date.pop(0)
                date.append(year)
                date = '/'.join(date)
            else:
                date = ""
            
            for author in xmlDoc.getElementsByTagName("Author"):
                affiliation_node = author.getElementsByTagName("Affiliation")
                if affiliation_node and "@" in getText(affiliation_node[0].childNodes):
                    address, email = parseAffiliation(affiliation_node)
                    author_with_email = getText(author.getElementsByTagName("ForeName")[0].childNodes) + " " + getText(author.getElementsByTagName("LastName")[0].childNodes) 
            
                    if any(x in address for x in ['United States', 'US', 'U.S.', 'United States of America', 'USA', 'U.S.A']):
                        data['Name'].append(author_with_email)
                        data['Email'].append(email)
                        data['Address'].append(address)
                        data['Keyword'].append(keyword)
                        data['Title'].append(title)
                        data['Date'].append(date)
                    else:
                        non_US['Name'].append(author_with_email)
                        non_US['Email'].append(email)
                        non_US['Address'].append(address)
                        non_US['Keyword'].append(keyword)
                        non_US['Title'].append(title)
                        non_US['Date'].append(date)
        except Exception as e:
            print(f'Error processing keyword "{keyword}": {e}')
            print("\n")
    
    print("\n")
    return data, non_US

def getKeywords(filename):
    try: 
        with open(filename, 'r') as file:
            raw_words = file.read();
            keywords = raw_words.split('\n')
        return keywords;
    except:
        return []
    
def searchPubmed(keywords, data, non_US, num_papers):
    for keyword in (keywords):
            tempData, tempNon_US = searchKeyword(keyword, num_papers)
            for key in data:
                data[key] += tempData[key]
                non_US[key] += tempNon_US[key]
        
    df = pd.DataFrame(data)
    non_US_df = pd.DataFrame(non_US).drop_duplicates(subset=['Email'])
    df.to_csv("./data/US_DATA.csv")
    non_US_df.to_csv("./data/INTERNATIONAL_DATA.csv")

if __name__ == '__main__':
    quit = False
    NCBI_API_KEY = "";
    os.makedirs('data', exist_ok=True)

    print("Welcome to the pubmed scraper!")
    print("Before using this tool, you must do the following: ")
    print("1. Register for a NCBI account at https://www.ncbi.nlm.nih.gov/")
    print("2. Click on your email at the top right of the dashboard")
    print("3. Click 'Account Settings'")
    print("3. Scroll down and copy down your NCBI API key")
    input("Press ENTER to continue: ")
    print('\n')

    if not os.path.isfile('./.env'):
        NCBI_API_KEY = input("Please enter your NCBI API key: ")
        with open('.env', 'w') as env:
            env.write(f'NCBI_API_KEY={NCBI_API_KEY}')
    else:
        with open('.env', 'r+') as env:
            NCBI_API_KEY = env.read()
            modify_api_key = input(f"Enter 'Y' if you would like to modify your current API key '{NCBI_API_KEY}': ")
            if modify_api_key.upper() == 'Y':
                NCBI_API_KEY = input("Please enter your NCBI API key: ")
                env.write(f'NCBI_API_KEY={NCBI_API_KEY}')

    if not os.path.isfile('user.json'):
        user['Email'] = input("Please enter your email: ")
        with open('user.json', 'w') as user_file:
            json.dump(user, user_file)
    else:
        with open('user.json', 'r+') as user_file:
            user = json.load(user_file)
            modify_email = input(f"Enter 'Y' if you would like to modify your current email '{user['Email']}': ")
            if modify_email.upper() == 'Y':
                user['Email'] = input("Please enter your email: ")
                json.dump(user, user_file)


    Entrez.api_key = os.getenv('NCBI_API_KEY')
    Entrez.email = user['Email']

    keywords = getKeywords('keywords.txt')
    
    if (len(keywords) == 0):
        print(f"Cannot open file 'keywords.txt'")
        quit = True

    print('\n')
    print("Note: You can modify your keywords in 'keywords.txt'")
    print('\n')


    while (quit == False):
        print("Commands (enter exactly as shown): quit, search, view_keywords, single_word_search")
        user_input = input("What would you like to do: ")

        if (user_input == 'quit'):
            quit = True
        elif (user_input == 'search'):
            num_papers = input("How many papers would you like to search for each keyword?: ")
            print("Searching for papers...\n")
            searchPubmed(keywords, data, non_US, num_papers)
            print("\nDone!\n")
            print("Data has been saved in the 'data' folder.\n")
        elif (user_input == 'view_keywords'):
            print("Keywords: ")
            for keyword in keywords:
                print(keyword)
            print('\n')
        elif (user_input == 'single_word_search'):
            keyword = input("Enter a keyword: ")
            num_papers = input("How many papers would you like to search for?: ")
            print("Searching for papers...\n")
            searchPubmed([keyword], data, non_US, num_papers)
            print("\nDone!\n")
            print("Data has been saved in the 'data' folder.\n")
        else:
            print("Invalid command. Please try again.")
            print('\n')
            print('\n')
            print('\n')
        

        
                   