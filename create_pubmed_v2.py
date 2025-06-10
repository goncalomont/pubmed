#!/usr/bin/env python3
# ---------------------------------
#  PubMed harvesting & enrichment
#  *includes Books & Documents*
# ---------------------------------

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

# ---------------------------------------------------------------------
# 0.  Logging
# ---------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('pubmed_parser.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)
warnings.filterwarnings("ignore", category=XMLParsedAsHTMLWarning)

# ---------------------------------------------------------------------
# 1.  Entrez + Unpaywall credentials
# ---------------------------------------------------------------------
Entrez.email = "casa5@ferring.com"
UNPAYWALL_EMAIL = "casa5@ferring.com"

# ---------------------------------------------------------------------
# 2.  Helper: Unpaywall OA / licence look-up
# ---------------------------------------------------------------------
def get_unpaywall_info(doi: str, email: str):
    if not doi:
        logger.warning("DOI not provided.")
        return None, None, None, None
    url = f"https://api.unpaywall.org/v2/{doi}?email={email}"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            result = response.json()
            if not result or not isinstance(result, dict):
                logger.warning(f"Empty or invalid API response for DOI {doi}")
                return None, None, None, None
            oa_status_str = "Open Access" if result.get("is_oa", False) else "Closed Access"
            oa_type = result.get("oa_status", "")
            licence = (result.get("best_oa_location") or {}).get("license", "")
            if not licence:
                for loc in result.get("oa_locations", []):
                    if loc.get("license"):
                        licence = loc["license"]; break
            best_url = (result.get("best_oa_location") or {}).get("url", "")
            if result.get("is_oa") and not licence:
                logger.warning(f"No licence found for DOI {doi}")
            return oa_status_str, oa_type, licence or "", best_url
        else:
            logger.warning(f"HTTP {response.status_code} from Unpaywall for DOI {doi}")
            return None, None, None, None
    except Exception as exc:
        logger.error(f"Unpaywall error for DOI {doi}: {exc}")
        return None, None, None, None
    finally:
        time.sleep(0.1)          # rate-limit 10 req/s

# ---------------------------------------------------------------------
# 3.  Query construction
# ---------------------------------------------------------------------
MeSH_QUERY = (
    '("Interleukin-11"[MeSH Terms] OR "IL-11"[Title/Abstract] OR "IL11"[Title/Abstract] '
    'OR "interleukin 11 receptor"[Title/Abstract] OR "IL11RA"[MeSH Terms] OR "IL11RA"[Title/Abstract]) OR '
    '("fibrosis"[MeSH Major Topic] OR "pulmonary fibrosis"[MeSH Terms] OR "non-alcoholic steatohepatitis"[MeSH Terms] '
    'OR "NASH"[Title/Abstract] OR "crohn disease"[MeSH Terms] OR "inflammatory bowel diseases"[MeSH Terms] '
    'OR "eosinophilic esophagitis"[MeSH Terms]) OR '
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

# ---------------------------------------------------------------------
# 4.  PubMed ESearch
# ---------------------------------------------------------------------
try:
    logger.info("Running PubMed search…")
    handle = Entrez.esearch(db="pubmed", term=query, retmax=200, sort="relevance")
    search_results = Entrez.read(handle)
    handle.close()
    pmids = search_results.get("IdList", [])
    logger.info(f"Found {len(pmids)} PMIDs")
    if not pmids:
        logger.warning("No results for the query.")
except Exception as exc:
    logger.error(f"ESearch error: {exc}")
    pmids = []

# ---------------------------------------------------------------------
# 5.  Article-type filter mapper
# ---------------------------------------------------------------------
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
    """Return semicolon-delimited label(s) emulating PubMed’s sidebar."""
    hits = set()
    # 5.1 publication-type matches
    low_pub = [p.lower() for p in pub_types]
    for filt, pts in ARTICLE_FILTER_MAP.items():
        if any(pt.lower() in low_pub for pt in pts):
            hits.add(filt)
    # 5.2 hybrid hedge for Systematic Review
    if ("systematic review" in low_pub or
        re.search(r"\bsystematic (literature )?(review|meta[- ]review|mapping review)\b",
                  f"{title} {abstract}", flags=re.I)):
        hits.add("Systematic Review")
    # 5.3 structural hedge for Books & Docs
    if is_book_xml:
        hits.add("Books & Documents")
    return "; ".join(sorted(hits)) or "Uncategorized"



# ---------------------------------------------------------------------
# 6.  Fetch PubMed XML (single call handles articles + books)
# ---------------------------------------------------------------------
articles = {"PubmedArticle": [], "PubmedBookArticle": []}
if pmids:
    try:
        fh = Entrez.efetch(db="pubmed", id=",".join(pmids),
                           rettype="xml", retmode="xml")
        articles = Entrez.read(fh)                     # Biopython handle → dict  :contentReference[oaicite:1]{index=1}
        fh.close()
    except Exception as exc:
        logger.error(f"EFetch error: {exc}")            # EFetch params per NCBI guide  :contentReference[oaicite:2]{index=2}

# ---------------------------------------------------------------------
# 7.  Parse records (journal + book)
# ---------------------------------------------------------------------
data = []
for art_type, article_list in articles.items():        # iterate over both keys
    is_book_xml = (art_type == "PubmedBookArticle")
    for article in article_list:
        article_data = {}
        try:
            # ---------- 7.1  Core citation block ----------
            container = (article["BookDocument"] if is_book_xml
                         else article["MedlineCitation"])
            article_data["PMID"] = str(container.get("PMID", ""))
            # ---------- 7.2  IDs ----------
            article_data["DOI"] = article_data["DOI_Link"] = None
            article_data["PMC_Link"] = None
            for aid in article["PubmedData"]["ArticleIdList"]:
                if aid.attributes["IdType"] == "doi":
                    article_data["DOI"] = str(aid)
                    article_data["DOI_Link"] = f"https://doi.org/{aid}"
                elif aid.attributes["IdType"] == "pmc":
                    article_data["PMC_Link"] = (
                        f"https://www.ncbi.nlm.nih.gov/pmc/articles/{aid}/"
                    )

            # ---------- 7.3  Titles, abstracts ----------
            if is_book_xml:
                # book title lives directly under <BookDocument><BookTitle>
                article_data["Title"] = container["Book"]["BookTitle"]
                article_data["Abstract"] = ""          # Books rarely have abstracts
            else:
                article_data["Title"] = container["Article"].get("ArticleTitle", "")
                abs_field = container["Article"].get("Abstract", {})
                text_fragments = abs_field.get("AbstractText", [])
                if isinstance(text_fragments, list):
                    clean_parts = []
                    for frag in text_fragments:
                        label = getattr(frag, "attributes", {}).get("Label", "")
                        cleaned = BeautifulSoup(str(frag), "html.parser").get_text()
                        clean_parts.append(
                            f"{label}: {cleaned}" if label else cleaned)
                    article_data["Abstract"] = " ".join(clean_parts)
                else:
                    article_data["Abstract"] = (
                        BeautifulSoup(str(text_fragments), "html.parser")
                        .get_text().strip()
                        if text_fragments else ""
                    )

            # ---------- 7.4  Journal / Book metadata ----------
            if is_book_xml:
                pub_date = container["Book"].get("PubDate", {})
                article_data["Journal"] = container["Book"].get("Publisher", "")
                article_data["Volume"] = article_data["Issue"] = ""
            else:
                journal = container["Article"]["Journal"]
                article_data["Journal"] = journal.get("Title", "")
                issue = journal.get("JournalIssue", {})
                article_data["Volume"] = issue.get("Volume", "")
                article_data["Issue"] = issue.get("Issue", "")
                pub_date = issue.get("PubDate", {})

            # ---------- 7.5  Date helpers ----------
            year = pub_date.get("Year", "")
            month = pub_date.get("Month", "")
            day = pub_date.get("Day", "")
            if year:
                if month:
                    try:
                        m_num = int(month) if month.isdigit() else parser.parse(month, fuzzy=True).month
                        article_data["PubDate"] = f"{year} {month_abbr[m_num]}{f' {int(day):02d}' if day else ''}"
                    except Exception:
                        article_data["PubDate"] = year
                else:
                    article_data["PubDate"] = year
            else:
                article_data["PubDate"] = ""
            article_data["PublicationYear"] = year
            try:
                article_data["PublicationMonth"] = month_name[int(month)] if month.isdigit() else (
                    month_name[parser.parse(month, fuzzy=True).month] if month else "")
            except Exception:
                article_data["PublicationMonth"] = ""
            article_data["PublicationDay"] = f"{int(day):02d}" if day.isdigit() else ""

            # ---------- 7.6  Language / Country ----------
            if is_book_xml:
                article_data["Language"] = container.get("Language", [""])[0]
                article_data["Country"] = ""
            else:
                article_data["Language"] = container["Article"].get("Language", [""])[0]
                article_data["Country"] = container.get("MedlineJournalInfo", {}).get("Country", "")

            # ---------- 7.7  Authors ----------
            if is_book_xml:
                authors = container["Book"].get("AuthorList", [])
            else:
                authors = container["Article"].get("AuthorList", [])
            names = [f"{a.get('ForeName','')} {a.get('LastName','')}".strip()
                     for a in authors if "LastName" in a]
            article_data["Authors"] = "; ".join(names)

            # ---------- 7.8  Keywords ----------
            kws = []
            if not is_book_xml and "KeywordList" in container:
                for kwlist in container["KeywordList"]:
                    kws.extend([kw for kw in kwlist if isinstance(kw, str)])
            article_data["AuthorKeywords"] = "; ".join(kws)

            # --- 7.9  Publication Types ---
            if is_book_xml:
                pt_list = container.get("PublicationTypeList", [])
            else:
                pt_list = container["Article"].get("PublicationTypeList", [])
            article_data["PublicationTypes"] = "; ".join(map(str, pt_list))

            # ---------- 7.10  MeSH ----------
            mesh_terms = []
            if not is_book_xml and "MeshHeadingList" in container:
                for mh in container["MeshHeadingList"]:
                    desc = str(mh["DescriptorName"])
                    quals = [str(q) for q in mh.get("QualifierName", [])]
                    if quals:
                        mesh_terms.extend([f"{desc}/{q}" for q in quals])
                    else:
                        mesh_terms.append(desc)
            article_data["MeSH_Terms"] = "; ".join(mesh_terms)

            # ---------- 7.11  OA / licence ----------
            oa, oa_type, licence, pdf_url = get_unpaywall_info(
                article_data["DOI"], UNPAYWALL_EMAIL)      # Unpaywall API  :contentReference[oaicite:3]{index=3}
            article_data.update({
                "AccessStatus": oa or "Unknown",
                "OA_Type": oa_type or "",
                "License": licence or "",
                "FreePDF_Link": pdf_url or ""
            })
            if licence.lower() in {"cc-by-nc", "cc-by-nd"}:
                logger.info(f"Skip PMID {article_data['PMID']} due to licence {licence}")
                continue

            # ---------- 7.12  Article-type filters ----------
            article_data["ArticleTypeFilters"] = detect_filters(
                [str(pt) for pt in pt_list],
                article_data["Title"],
                article_data["Abstract"],
                is_book_xml=is_book_xml
            )

            data.append(article_data)

        except Exception as exc:
            logger.error(f"Parsing error on PMID {article_data.get('PMID','?')}: {exc}")
            continue

# ---------------------------------------------------------------------
# 8.  DataFrame + CSV
# ---------------------------------------------------------------------
df = pd.DataFrame(data)

column_order = [
    "PMID", "Title", "Abstract", "Journal", "PubDate",
    "PublicationYear", "PublicationMonth", "PublicationDay",
    "Authors", "AuthorKeywords", "DOI", "DOI_Link", "PMC_Link",
    "Volume", "Issue", "Pagination", "ELocationID",
    "Language", "Country", "PublicationTypes", "ArticleTypeFilters",
    "MeSH_Terms", "AccessStatus", "OA_Type", "License", "FreePDF_Link"
]
df = df[[c for c in column_order if c in df.columns]]
logger.info(f"Processed {len(df)} records; writing CSV")
df.to_csv("pubmed_articles_Metadata_with_license.csv", index=False)
