from Bio import Entrez
import os
import pandas as pd
from xml.dom import minidom    
import tqdm

Entrez.api_key = os.environ['NCBI_API_KEY']

# Provide your email address to NCBI
Entrez.email = "hudsontao@gmail.com"

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

if __name__ == '__main__':
    os.makedirs('data', exist_ok=True)
    
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
    
    keywords = [
        "PCR+Master+Mix",
        "Real-time+PCR+reagents",
        "qPCR+master+mix",
        "Taq+DNA+polymerase",
        "PCR+amplification+kits",
        "Multiplex+PCR",
        "High-fidelity+PCR",
        "Nucleic+Acid+Gel+Stains",
        "SYBR+Green",
        "Propidium+Iodide",
        "GelRed",
        "Nucleic+acid+visualization",
        "Agarose+gel+electrophoresis",
        "DNA/RNA+staining",
        "Nucleic+Acid+Extractor",
        "Automated+DNA/RNA+isolation",
        "Magnetic+bead-based+extraction",
        "Viral+RNA+extraction",
        "Plasmid+DNA+isolation",
        "Sample+preparation+technologies",
        "Next+Generation+Sequencing+(NGS)",
        "High-throughput+sequencing+methods",
        "Illumina+sequencing+technology",
        "Oxford+Nanopore+sequencing",
        "Single-cell+sequencing",
        "Whole-genome+sequencing+(WGS)",
        "RNA-seq",
        "Metagenomics",
        "Gene+Synthesis",
        "Synthetic+biology",
        "Gene+assembly",
        "Oligonucleotide+synthesis",
        "CRISPR-Cas9+gene+editing",
        "Synthetic+genomics",
        "Protein+engineering",
        "Antibody",
        "Monoclonal+antibodies",
        "Polyclonal+antibodies",
        "Antibody+production+and+purification",
        "Recombinant+antibodies",
        "Therapeutic+antibodies",
        "Antibody-drug+conjugates+(ADCs)",
        "Immunoassays",
        "Flow+cytometry+antibodies",
        "Transfection+Reagent",
        "Lipofection",
        "Electroporation",
        "Non-viral+gene+delivery+systems",
        "siRNA+delivery",
        "CRISPR-Cas9+delivery+systems",
        "Lipid+nanoparticles",
        "Calcium+phosphate+transfection",
        "Viral+vector+production",
        "Lentiviral+Packaging",
        "Lentiviral+vectors",
        "Viral+transduction",
        "Gene+therapy+vectors",
        "Stable+cell+line+generation",
        "Viral+vector+purification"
    ]
    
    for keyword in (keywords):
        tempData, tempNon_US = searchKeyword(keyword, 100)
        for key in data:
            data[key] += tempData[key]
            non_US[key] += tempNon_US[key]
    
    df = pd.DataFrame(data)
    non_US_df = pd.DataFrame(non_US).drop_duplicates(subset=['Email'])
    df.to_csv("./data/data.csv")
    non_US_df.to_csv("./data/other_data.csv")
                   