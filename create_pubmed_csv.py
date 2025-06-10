import logging
from Bio import Entrez
import pandas as pd
from bs4 import BeautifulSoup, XMLParsedAsHTMLWarning
import warnings
from dateutil import parser
import re
from calendar import month_name, month_abbr
import requests
import time

# Configure logging at the start
logging.basicConfig(
    level=logging.INFO,  # Set minimum log level to INFO to capture progress messages
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('pubmed_parser.log'),  # Log to a file
        logging.StreamHandler()  # Log to console
    ]
)
logger = logging.getLogger(__name__)

# Suppress XMLParsedAsHTMLWarning
warnings.filterwarnings("ignore", category=XMLParsedAsHTMLWarning)

# Set email (required by NCBI and Unpaywall APIs)
Entrez.email = "casa5@ferring.com"
UNPAYWALL_EMAIL = "casa5@ferring.com"

# Function to fetch license and open access information from Unpaywall
def get_unpaywall_info(doi, email):
    if not doi:
        logger.warning("DOI not provided.")
        return None, None, None, None
    url = f"https://api.unpaywall.org/v2/{doi}?email={email}"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            try:
                result = response.json()
                if not result or not isinstance(result, dict):
                    logger.warning(f"Empty or invalid API response for DOI {doi}")
                    return None, None, None, None
                oa_status = result.get("is_oa", False)
                oa_status_str = "Open Access" if oa_status else "Closed Access"
                oa_type = result.get("oa_status", "")
                # Extract license from best_oa_location or oa_locations
                license = None
                best_oa_location = result.get("best_oa_location", {})
                if best_oa_location:
                    license = best_oa_location.get("license", None)
                if not license:
                    # Fallback to other oa_locations
                    for location in result.get("oa_locations", []):
                        if location.get("license"):
                            license = location.get("license")
                            break
                license = license or ""
                best_oa_url = best_oa_location.get("url", "") if best_oa_location else ""
                if oa_status and not license:
                    logger.warning(f"No license was found for DOI {doi} (OA: {oa_status})")
                return oa_status_str, oa_type, license, best_oa_url
            except ValueError:
                logger.warning(f"Invalid JSON response for DOI {doi}")
                return None, None, None, None
        else:
            logger.warning(f"HTTP Error {response.status_code} when querying Unpaywall API for DOI {doi}")
            return None, None, None, None
    except Exception as e:
        logger.error(f"Error querying Unpaywall for DOI {doi}: {e}")
        return None, None, None, None
    finally:
        time.sleep(0.1)  # Respect Unpaywall rate limits (10 requests/second)

# Set Query filters based on MeSH and Article Type
MeSH_QUERY = ('("Interleukin-11"[MeSH Terms] OR "IL-11"[Title/Abstract] OR "IL11"[Title/Abstract] OR "interleukin 11 receptor"[Title/Abstract] OR "IL11RA"[MeSH Terms] OR "IL11RA"[Title/Abstract]) OR '
              '("fibrosis"[MeSH Major Topic] OR "pulmonary fibrosis"[MeSH Terms] OR "non-alcoholic steatohepatitis"[MeSH Terms] OR "NASH"[Title/Abstract] OR "crohn disease"[MeSH Terms] OR "inflammatory bowel diseases"[MeSH Terms] OR "eosinophilic esophagitis"[MeSH Terms]) OR '
              '("immunology"[MeSH Major Topic] OR "cytokines"[MeSH Terms] OR "inflammation"[MeSH Major Topic]) OR '
              '("transcriptomic*"[Title/Abstract] OR "proteomic*"[Title/Abstract] OR "single cell"[Title/Abstract]) OR '
              '("fibroblast*"[Title/Abstract] OR "myofibroblast*"[Title/Abstract] OR "macrophage*"[Title/Abstract])'
             )

SPECIES_FILTER = '"animals"[MeSH Terms]'

ARTICLE_TYPE_FILTERS = (
    "clinical trial[pt] OR randomized controlled trial[pt] OR meta-analysis[pt] OR "
    "systematic review[pt] OR systematic[sb] OR observational study[pt] OR review[pt] OR "
    "case reports[pt] OR practice guideline[pt] OR \"pubmed books\"[sb]"
)
query = f"({MeSH_QUERY}) AND ({SPECIES_FILTER}) AND ({ARTICLE_TYPE_FILTERS})"

# Step 1: Search PubMed
try:
    logger.info("Running PubMed search...")
    search_handle = Entrez.esearch(db="pubmed", term=query, retmax=200, sort="relevance")
    search_results = Entrez.read(search_handle)
    search_handle.close()
    pmids = search_results["IdList"]
    logger.info(f"Found {len(pmids)} PMIDs: {pmids}")
    if not pmids:
        logger.warning("No results found for the query provided.")
except Exception as e:
    logger.error(f"Error searching PubMed: {str(e)}")
    pmids = []

ARTICLE_FILTER_MAP = {
    "Books & Documents": ["Book", "Book Chapter"],
    "Clinical Trial": [
        "Clinical Trial", "Clinical Trial, Phase I", "Clinical Trial, Phase II",
        "Clinical Trial, Phase III", "Clinical Trial, Phase IV",
        "Controlled Clinical Trial", "Pragmatic Clinical Trial"
    ],
    "Randomized Controlled Trial": ["Randomized Controlled Trial"],
    "Meta-Analysis": ["Meta-Analysis", "Network Meta-Analysis"],
    "Observational Study": ["Observational Study", "Observational Study, Veterinary"],
    "Review": ["Review"],
    "Case Reports": ["Case Reports"],
    "Practice Guideline": ["Practice Guideline"],
}

def detect_filters(pub_types, title, abstract, is_book_xml=False):
    """Return a semicolon-delimited list that mirrors PubMed’s article-type sidebar."""
    hits = set()
    # publication-type matches
    for filt, pts in ARTICLE_FILTER_MAP.items():
        if any(pt.lower() in (p.lower() for p in pub_types) for pt in pts):
            hits.add(filt)

    # hybrid hedge for Systematic Review
    if ("systematic review" in (p.lower() for p in pub_types) or
        re.search(r"\bsystematic (literature )?(review|meta[- ]review|mapping review)\b",
                  f"{title} {abstract}", flags=re.I)):
        hits.add("Systematic Review")

    if is_book_xml:
        hits.add("Books & Documents")

    return "; ".join(sorted(hits)) or "Uncategorized"

# Step 2: Fetch full records in XML format

data = []
if pmids:
    try:
        fetch_handle = Entrez.efetch(db="pubmed", id=",".join(pmids), rettype="xml", retmode="xml")
        articles = Entrez.read(fetch_handle)
        fetch_handle.close()
    except Exception as e:
        logger.error(f"Error fetching or reading records: {e}")
        articles = {'PubmedArticle': []}
else:
    articles = {'PubmedArticle': []}
    

# Step 3: Parse Articles

data = []
for article in articles['PubmedArticle']:
    article_data = {}
    try:
        citation = article["MedlineCitation"]
        article_data["PMID"] = str(citation.get("PMID", ""))

        # Article ID Metadata (DOI, DOI link, PMC)
        article_data["DOI"] = None
        article_data["DOI_Link"] = None
        article_data["PMC_Link"] = None
        for aid in article["PubmedData"]["ArticleIdList"]:
            id_type = aid.attributes["IdType"]
            if id_type == "doi":
                article_data["DOI"] = str(aid)
                article_data["DOI_Link"] = f"https://doi.org/{aid}"
            elif id_type == "pmc":
                pmc_id = str(aid)
                article_data["PMC_Link"] = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}/"

     # Article bibliographic Metadata
        article_data["Title"] = citation["Article"].get("ArticleTitle", "")

        if "Abstract" in citation["Article"]:
            abstract_text = citation["Article"]["Abstract"].get("AbstractText", [])
            abstract_parts = []
            if isinstance(abstract_text, list):
                for part in abstract_text:
                    label = part.attributes.get("Label", "") if hasattr(part, "attributes") else ""
                    text = str(part)
                    if text:
                        soup = BeautifulSoup(text, "html.parser")
                        clean_text = soup.get_text().strip()
                        if label and clean_text:
                            abstract_parts.append(f"{label}: {clean_text}")
                        elif clean_text:
                            abstract_parts.append(clean_text)
            elif isinstance(abstract_text, str) and abstract_text.strip():
                soup = BeautifulSoup(abstract_text, "html.parser")
                clean_text = soup.get_text().strip()
                abstract_parts.append(clean_text)
            article_data["Abstract"] = " ".join(abstract_parts) if abstract_parts else ""
        else:
            article_data["Abstract"] = ""

        if "Journal" not in citation["Article"] or citation["Article"]["Journal"] is None:
            logger.warning(f"Journal missing or None for PMID {article_data['PMID']}. Setting related fields to empty.")
            article_data["Journal"] = ""
            article_data["Volume"] = ""
            article_data["Issue"] = ""
            article_data["PubDate"] = ""
            article_data["PublicationYear"] = ""
            article_data["PublicationMonth"] = ""
            article_data["PublicationDay"] = ""
        else:
            article_data["Journal"] = citation["Article"]["Journal"].get("Title", "")

            journal_issue = citation["Article"]["Journal"].get("JournalIssue", {})
            article_data["Volume"] = journal_issue.get("Volume", "")  
            article_data["Issue"] = journal_issue.get("Issue", "")

        # Publication Date Handling
        pub_date = journal_issue.get("PubDate", {})
        year = pub_date.get("Year", "")
        month = pub_date.get("Month", "")
        day = pub_date.get("Day", "")

        # PubDate
        if year:
            if month:
                try:
                    month_num = int(month) if month.isdigit() else parser.parse(month, fuzzy=True).month
                    if day:
                        day_num = int(day)
                        article_data["PubDate"] = f"{year} {month_abbr[month_num]} {day_num:02d}"  # Ex.: "2025 Jun 05"
                    else:
                        article_data["PubDate"] = f"{year} {month_abbr[month_num]}"  # Ex.: "1998 Oct"
                except (ValueError, TypeError):
                    article_data["PubDate"] = year  
            else:
                article_data["PubDate"] = year  # Only year if month missing
        else:
            article_data["PubDate"] = ""

        # Publication Year
        article_data["PublicationYear"] = year if year else ""

        # Publication Month (long)
        if month:
            try:
                month_num = int(month) if month.isdigit() else parser.parse(month, fuzzy=True).month
                article_data["PublicationMonth"] = month_name[month_num]  # Ex.: "October"
            except (ValueError, TypeError):
                article_data["PublicationMonth"] = ""
        else:
            article_data["PublicationMonth"] = ""

        # Publication Day
        article_data["PublicationDay"] = f"{int(day):02d}" if day and day.isdigit() else ""

        # Language and Country
        article_data["Language"] = citation["Article"].get("Language", [""])[0]
        article_data["Country"] = citation.get("MedlineJournalInfo", {}).get("Country", "")

        # Authors, Keyword & Publication Type
        if "AuthorList" in citation["Article"]:
            authors = citation["Article"]["AuthorList"]
            author_names = [
                f"{a['ForeName']} {a['LastName']}" 
                for a in authors if "ForeName" in a and "LastName" in a
            ]
            article_data["Authors"] = "; ".join(author_names) if author_names else ""
        else:
            article_data["Authors"] = ""

        author_keywords = []
        if "KeywordList" in citation:
            for kwlist in citation["KeywordList"]:
                for kw in kwlist:
                    if isinstance(kw, str):
                        author_keywords.append(kw)
        article_data["AuthorKeywords"] = "; ".join(author_keywords) if author_keywords else ""

        pt_list = citation["Article"].get("PublicationTypeList", [])
        article_data["PublicationTypes"] = "; ".join(str(pt) for pt in pt_list) if pt_list else ""

        # MeSH Terms Indexing with Qualifiers
        mesh_terms = []
        if "MeshHeadingList" in citation:
                for mh in citation["MeshHeadingList"]:
                    descriptor = str(mh["DescriptorName"])
                    qualifiers = mh.get("QualifierName", [])
                    if qualifiers:
                         for qualifier in qualifiers:
                            mesh_terms.append(f"{descriptor}/{str(qualifier)}")
                    else:
                        mesh_terms.append(descriptor)
        article_data["MeSH_Terms"] = "; ".join(mesh_terms) if mesh_terms else ""

        # License data from Unpaywall
        oa_status_str, oa_type, license, best_oa_url = get_unpaywall_info(article_data["DOI"], UNPAYWALL_EMAIL)
        article_data["AccessStatus"] = oa_status_str or "Unknown"
        article_data["OA_Type"] = oa_type or ""
        article_data["License"] = license or ""
        article_data["FreePDF_Link"] = best_oa_url or ""

        # Filter articles based on license
        if license and license.lower() in ["cc-by-nc", "cc-by-nd"]:
            logger.info(f"Skipping article {article_data.get('PMID', 'Unknown')} due to license: {license}")
            continue

        data.append(article_data)

    except Exception as e:
        logger.error(f"Error processing article PMID {article_data.get('PMID', 'Unknown')}: {e}")
        continue

# Step 4: Create DataFrame
df = pd.DataFrame(data)

column_order = [
    "PMID", "Title", "Abstract", "Journal", "PubDate", "PublicationYear",
    "PublicationMonth", "PublicationDay", "Authors", "AuthorKeywords", "DOI",
    "DOI_Link", "PMC_Link", "Volume", "Issue", "Pagination", "ELocationID",
    "Language", "Country", "PublicationTypes", "MeSH_Terms",
    "AccessStatus", "OA_Type", "License", "FreePDF_Link"
]

df = df[[col for col in column_order if col in df.columns]]

# Step 5: Log results and save
logger.info(f"Processed {len(df)} articles successfully.")
# Optionally, log the DataFrame preview (limit to first few rows for brevity)
logger.debug(f"DataFrame preview:\n{df.head().to_string()}")

#display(df)

# Save to CSV
df.to_csv("pubmed_articles_Metadata_with_license.csv", index=False)