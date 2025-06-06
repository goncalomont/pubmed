from Bio import Entrez, Medline
import time
import xml.etree.ElementTree as ET
import csv
import requests
from bs4 import BeautifulSoup

# Set your email for NCBI Entrez
Entrez.email = "REPLACE_WITH_YOUR_EMAIL@example.com"  # Replace with a valid email

def parse_pubmed_documents(filepath="pubmed_documents.txt"):
    """
    Parses the pubmed_documents.txt file and returns a list of articles.
    Each article is a dictionary with 'PMID', 'Title', and 'Abstract'.
    """
    articles = []
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            content = f.read()

        raw_articles = content.strip().split("---\n")

        for raw_article in raw_articles:
            if not raw_article.strip():
                continue

            article_data = {}
            lines = raw_article.strip().split("\n")

            for line in lines:
                if line.startswith("PMID:"):
                    article_data["PMID"] = line.replace("PMID:", "").strip()
                elif line.startswith("Title:"):
                    article_data["Title"] = line.replace("Title:", "").strip()
                elif line.startswith("Abstract:"):
                    # Handle multi-line abstracts if they were to occur
                    current_abstract = line.replace("Abstract:", "").strip()
                    if "Abstract" in article_data:
                         article_data["Abstract"] += " " + current_abstract
                    else:
                        article_data["Abstract"] = current_abstract

            if "PMID" in article_data and "Title" in article_data and "Abstract" in article_data:
                articles.append(article_data)
            else:
                print(f"Warning: Could not parse a record properly. Content: {raw_article}")

    except FileNotFoundError:
        print(f"Error: File not found - {filepath}")
        return []
    except Exception as e:
        print(f"Error parsing {filepath}: {e}")
        return []
    return articles

def fetch_article_metadata(pmid):
    """
    Fetches detailed metadata for a given PMID from PubMed using efetch.
    Returns the XML response as a string.
    """
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="full", retmode="xml")
        xml_data = handle.read()
        handle.close()
        return xml_data
    except Exception as e:
        print(f"Error fetching metadata for PMID {pmid}: {e}")
        return None

def extract_pmcid_and_link(xml_data, pmid):
    """
    Tries to extract PMCID and a potential full-text link from the XML metadata.
    """
    pmcid = None
    full_text_link = None

    if not xml_data:
        return pmcid, full_text_link

    try:
        root = ET.fromstring(xml_data)
        # Find PMCID
        # Common locations for PMCID
        pmcid_element = root.find(".//ArticleId[@IdType='pmc']")
        if pmcid_element is not None and pmcid_element.text:
            pmcid = pmcid_element.text
            if not pmcid.startswith("PMC"):
                pmcid = "PMC" + pmcid # Ensure it has PMC prefix
            full_text_link = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
        else:
            # Fallback: sometimes PMCID is within MedlineCitation/Article/ArticleIdList
            article_id_list = root.find(".//MedlineCitation/Article/ArticleIdList")
            if article_id_list is not None:
                for aid in article_id_list.findall("ArticleId"):
                    if aid.get("IdType") == "pmc" and aid.text:
                        pmcid = aid.text
                        if not pmcid.startswith("PMC"):
                             pmcid = "PMC" + pmcid
                        full_text_link = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
                        break
    except ET.ParseError as e:
        print(f"Error parsing XML for PMID {pmid}: {e}")
    except Exception as e:
        print(f"An unexpected error occurred during XML parsing for PMID {pmid}: {e}")

    return pmcid, full_text_link

def main():
    parsed_articles = parse_pubmed_documents()
    if not parsed_articles:
        print("No articles parsed from pubmed_documents.txt. Exiting.")
        return

    print(f"Found {len(parsed_articles)} articles in pubmed_documents.txt. Fetching metadata...\n")

    # Process all articles
    articles_to_process = parsed_articles
    # print(f"Processing the first {len(articles_to_process)} articles for this test run.\n") # Comment out or remove test line

    processed_articles_data = [] # To store data for CSV writing

    for i, article_info in enumerate(articles_to_process):
        pmid = article_info.get("PMID")
        title = article_info.get("Title", "No title found")
        abstract = article_info.get("Abstract", "No abstract found")

        if not pmid:
            print("Skipping article with no PMID.")
            continue

        print(f"Processing article {i+1}/{len(articles_to_process)}: PMID {pmid}")
        xml_data = fetch_article_metadata(pmid)

        pmcid = None
        full_text_link = None
        full_text_content = "" # Default to empty string

        if xml_data:
            pmcid, full_text_link = extract_pmcid_and_link(xml_data, pmid)
            if pmcid:
                print(f"  PMCID: {pmcid} (Potential Open Access)")
                print(f"  Full-text link: {full_text_link}")

                if full_text_link:
                    try:
                        print(f"    Fetching full text from: {full_text_link}")
                        headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'}
                        response = requests.get(full_text_link, headers=headers, timeout=10)
                        response.raise_for_status() # Raise HTTPError for bad responses (4XX or 5XX)

                        content_type = response.headers.get('Content-Type', '').lower()

                        if 'text/html' in content_type:
                            soup = BeautifulSoup(response.content, 'lxml')

                            article_body = soup.find('div', class_='article-content') # Common on PMC
                            if not article_body:
                                article_body = soup.find('div', class_='rendered_body') # Another common one
                            if not article_body:
                                article_body = soup.find('body')

                            if article_body:
                                full_text_content = article_body.get_text(separator='\n\n', strip=True)
                                full_text_content = full_text_content[:5000] # Limit length
                                print(f"    Successfully extracted ~{len(full_text_content)} chars of text content.")
                            else:
                                print("    Could not find main article content body/div in HTML.")

                        elif 'text/plain' in content_type:
                            full_text_content = response.text[:5000]
                            print(f"    Successfully extracted ~{len(full_text_content)} chars of plain text.")

                        elif 'application/pdf' in content_type:
                            print(f"    Link is a PDF, skipping direct text extraction: {full_text_link}")
                            full_text_content = "[PDF content not extracted]"

                        elif 'application/xml' in content_type or 'text/xml' in content_type:
                             # Basic XML text extraction if it's not the PMC HTML page but some other XML
                            try:
                                xml_root = ET.fromstring(response.content)
                                # Attempt to get all text nodes, crude but better than nothing
                                text_parts = [elem.text for elem in xml_root.iter() if elem.text]
                                full_text_content = " ".join(text_parts).strip()
                                full_text_content = full_text_content[:5000]
                                if full_text_content:
                                    print(f"    Successfully extracted ~{len(full_text_content)} chars from XML response.")
                                else:
                                    print("    XML response, but no text content found.")
                            except ET.ParseError:
                                print("    Could not parse XML from direct link.")
                            full_text_content = full_text_content if full_text_content else "[XML content not easily parsable to plain text]"

                        else:
                            print(f"    Unsupported content type: {content_type}")
                            full_text_content = f"[Unsupported content type: {content_type}]"

                    except requests.exceptions.Timeout:
                        print(f"    Timeout while fetching: {full_text_link}")
                        full_text_content = "[Fetching timed out]"
                    except requests.exceptions.HTTPError as http_err:
                        print(f"    HTTP error occurred: {http_err} for {full_text_link}")
                        full_text_content = f"[HTTP error: {http_err.response.status_code}]"
                    except requests.exceptions.RequestException as req_err:
                        print(f"    Request error occurred: {req_err} for {full_text_link}")
                        full_text_content = "[Request error]"
                    except Exception as e:
                        print(f"    An unexpected error occurred during full text fetching/parsing: {e}")
                        full_text_content = "[Error during processing]"
            else:
                print(f"  PMCID: Not found. Cannot attempt full text download without PMCID link.")
                full_text_content = "[No PMCID link for full text]"
        else:
            print(f"  Could not fetch metadata for PMID {pmid}.")
            full_text_content = "[Metadata fetch failed]"

        is_open_access = True if pmcid and full_text_link else False

        processed_articles_data.append({
            "PMID": pmid,
            "Title": title,
            "Abstract": abstract,
            "IsOpenAccess": is_open_access,
            "FullTextContent": full_text_content if full_text_content else "[No content extracted or not applicable]"
        })

        print("-" * 20)

        # Respect NCBI API usage guidelines
        time.sleep(0.5)

    # Write to CSV
    csv_file_path = "pubmed_articles.csv"
    csv_header = ["PMID", "Title", "Abstract", "IsOpenAccess", "FullTextContent"]

    try:
        with open(csv_file_path, "w", newline="", encoding="utf-8") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_header)
            writer.writeheader()
            for article_data in processed_articles_data:
                writer.writerow(article_data)
        print(f"\nSuccessfully wrote {len(processed_articles_data)} articles to {csv_file_path}")
    except IOError:
        print(f"Error: Could not write to CSV file {csv_file_path}")
    except Exception as e:
        print(f"An unexpected error occurred during CSV writing: {e}")

if __name__ == "__main__":
    main()
