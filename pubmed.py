from Bio import Entrez
import os

Entrez.api_key = os.environ['NCBI_API_KEY']

# Provide your email address to NCBI
Entrez.email = "hudsontao@gmail.com"

# Define the search term
search_term = "BPA"

# Perform the PubMed search
handle = Entrez.esearch(db="pubmed", term=search_term, retmax=100)
record = Entrez.read(handle)
handle.close()

# Get the list of PubMed IDs (PMID) for the search results
pmid_list = record["IdList"]
pmid_list.append("35231536")

# Fetch information for each article
for pmid in pmid_list:
    handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
    article_info = handle.read( )
    handle.close()

    isAD = False

    for line in article_info.split('\n'):
        if isAD == True:
            if line.startswith("   "):
                print(line)
                continue
            else:
                isAD = False
        if line.startswith("FAU"):
            author = line[6:]
            print("Author:", author)
        elif line.startswith("AD"):
            affiliation = line[6:]
            print("Affiliation:", affiliation)
            isAD = True

    print("-----------------------------------------------------------------------------------")
