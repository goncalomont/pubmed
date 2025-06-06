from Bio import Entrez
from Bio import Medline

# Set your email for NCBI Entrez
Entrez.email = "REPLACE_WITH_YOUR_EMAIL@example.com"  # Replace with a valid email

# Search PubMed for articles (e.g., about "cancer")
handle = Entrez.esearch(db="pubmed", term="cancer", retmax="100")
record = Entrez.read(handle)
handle.close()
idlist = record["IdList"]

# Fetch details for these IDs
handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")
records = Medline.parse(handle)

records = list(records)  # Convert iterator to list to save multiple times

# Save titles, abstracts, and PMIDs to a file
with open("pubmed_documents.txt", "w") as f:
    for record in records:
        pmid = record.get("PMID", "No PMID found")
        title = record.get("TI", "No title found")
        abstract = record.get("AB", "No abstract found")
        f.write(f"PMID: {pmid}\n")
        f.write(f"Title: {title}\n")
        f.write(f"Abstract: {abstract}\n")
        f.write("---\n")

print("Fetched 100 documents and saved PMIDs, titles, and abstracts to pubmed_documents.txt")
