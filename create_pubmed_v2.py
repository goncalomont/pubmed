#!/usr/bin/env python3
# -----------------------------------------------------
#  PubMed + PMC harvester
#  Part 1: setup, search, article-type mapper
# -----------------------------------------------------

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

# ------------------------------------------------------------------
# 0.  Logging
# ------------------------------------------------------------------
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

# ------------------------------------------------------------------
# 1.  Credentials
# ------------------------------------------------------------------
Entrez.email = "casa5@ferring.com"          # required by NCBI
UNPAYWALL_EMAIL = "casa5@ferring.com"       # required by Unpaywall

# ------------------------------------------------------------------
# 2.  Helper: Unpaywall OA / licence
# ------------------------------------------------------------------
def get_unpaywall_info(doi: str, email: str):
    if not doi:
        logger.warning("DOI not provided.")
        return None, None, None, None
    url = f"https://api.unpaywall.org/v2/{doi}?email={email}"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code != 200:
            logger.warning(f"Unpaywall HTTP {r.status_code} for DOI {doi}")
            return None, None, None, None
        obj = r.json()
        if not obj or not isinstance(obj, dict):
            logger.warning(f"Invalid Unpaywall JSON for DOI {doi}")
            return None, None, None, None
        is_oa = obj.get("is_oa", False)
        status = "Open Access" if is_oa else "Closed Access"
        oa_type = obj.get("oa_status", "")
        licence = (obj.get("best_oa_location") or {}).get("license", "")
        if not licence:
            for loc in obj.get("oa_locations", []):
                if loc.get("license"):
                    licence = loc["license"]; break
        pdf_url = (obj.get("best_oa_location") or {}).get("url", "")
        if is_oa and not licence:
            logger.warning(f"No licence found for DOI {doi}")
        return status, oa_type, licence or "", pdf_url
    except Exception as exc:
        logger.error(f"Unpaywall error for DOI {doi}: {exc}")
        return None, None, None, None
    finally:
        time.sleep(0.1)        # 10 req/s limit

# ------------------------------------------------------------------
# 3.  Search query
# ------------------------------------------------------------------
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

# ------------------------------------------------------------------
# 4.  ESearch
# ------------------------------------------------------------------
try:
    logger.info("Running PubMed search …")
    h = Entrez.esearch(db="pubmed", term=query, retmax=200, sort="relevance")
    res = Entrez.read(h); h.close()
    pmids = res.get("IdList", [])
    logger.info(f"{len(pmids)} PMIDs retrieved")
    if not pmids:
        logger.warning("Query returned zero records.")
except Exception as exc:
    logger.error(f"ESearch failed: {exc}")
    pmids = []

# ------------------------------------------------------------------
# 5.  Article-type → Sidebar label mapper
# ------------------------------------------------------------------
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
    """Return semicolon-delimited sidebar labels (Clinical Trial; Systematic Review; …)."""
    hits = set()
    low_pub = [p.lower() for p in pub_types]
    for label, pts in ARTICLE_FILTER_MAP.items():
        if any(pt.lower() in low_pub for pt in pts):
            hits.add(label)
    # Hybrid hedge for Systematic Review
    if ("systematic review" in low_pub or
        re.search(r"\bsystematic (literature )?(review|meta[- ]review|mapping review)\b",
                  f"{title} {abstract}", flags=re.I)):
        hits.add("Systematic Review")
    if is_book_xml:
        hits.add("Books & Documents")
    return "; ".join(sorted(hits)) or "Uncategorized"
# -----------------------------------------------------
# 6.  Fetch PubMed XML (journal + book in one call)
# -----------------------------------------------------
articles = {"PubmedArticle": [], "PubmedBookArticle": []}
if pmids:
    try:
        fh = Entrez.efetch(db="pubmed", id=",".join(pmids),
                           rettype="xml", retmode="xml")
        articles = Entrez.read(fh)      # Biopython parses to dict
        fh.close()
    except Exception as exc:
        logger.error(f"EFetch error: {exc}")

# -----------------------------------------------------
# 7.  Parse records
# -----------------------------------------------------
data = []
for art_type, art_list in articles.items():             # iterate over both keys
    is_book_xml = (art_type == "PubmedBookArticle")
    for article in art_list:
        article_data = {}
        try:
            # 7.1  Core citation container
            container = article["BookDocument"] if is_book_xml else article["MedlineCitation"]
            article_data["PMID"] = str(container.get("PMID", ""))

            # 7.2  IDs (DOI, PMC)
            article_data.update({"DOI": None, "DOI_Link": None, "PMC_Link": None})
            for aid in article["PubmedData"]["ArticleIdList"]:
                idt = aid.attributes["IdType"]
                if idt == "doi":
                    article_data["DOI"] = str(aid)
                    article_data["DOI_Link"] = f"https://doi.org/{aid}"
                elif idt == "pmc":
                    pmc_id = str(aid)
                    article_data["PMC_Link"] = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}/"

            # 7.3  Titles & abstracts
            if is_book_xml:
                article_data["Title"] = container["Book"].get("BookTitle", "")
                article_data["Abstract"] = ""
            else:
                article_data["Title"] = container["Article"].get("ArticleTitle", "")
                abs_field = container["Article"].get("Abstract", {})
                atext = abs_field.get("AbstractText", [])
                parts = []
                if isinstance(atext, list):
                    for seg in atext:
                        label = getattr(seg, "attributes", {}).get("Label", "")
                        txt = BeautifulSoup(str(seg), "html.parser").get_text().strip()
                        if txt:
                            parts.append(f"{label}: {txt}" if label else txt)
                elif isinstance(atext, str):
                    parts.append(BeautifulSoup(atext, "html.parser").get_text().strip())
                article_data["Abstract"] = " ".join(parts)

            # 7.4  Journal / Book metadata
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

            # 7.5  Date helpers
            year = pub_date.get("Year", "")
            month = pub_date.get("Month", "")
            day = pub_date.get("Day", "")
            if year:
                if month:
                    try:
                        mnum = int(month) if month.isdigit() else parser.parse(month, fuzzy=True).month
                        article_data["PubDate"] = f"{year} {month_abbr[mnum]}{f' {int(day):02d}' if day else ''}"
                    except Exception:
                        article_data["PubDate"] = year
                else:
                    article_data["PubDate"] = year
            else:
                article_data["PubDate"] = ""
            article_data["PublicationYear"] = year
            try:
                article_data["PublicationMonth"] = (
                    month_name[int(month)] if month.isdigit()
                    else month_name[parser.parse(month, fuzzy=True).month] if month else ""
                )
            except Exception:
                article_data["PublicationMonth"] = ""
            article_data["PublicationDay"] = f"{int(day):02d}" if day and day.isdigit() else ""

            # 7.6  Language / Country
            if is_book_xml:
                article_data["Language"] = container.get("Language", [""])[0]
                article_data["Country"] = ""
            else:
                article_data["Language"] = container["Article"].get("Language", [""])[0]
                article_data["Country"] = container.get("MedlineJournalInfo", {}).get("Country", "")

            # 7.7  Authors
            alist = (container["Book"].get("AuthorList", []) if is_book_xml
                     else container["Article"].get("AuthorList", []))
            names = [f"{a.get('ForeName','')} {a.get('LastName','')}".strip()
                     for a in alist if "LastName" in a]
            article_data["Authors"] = "; ".join(names)

            # 7.8  Keywords
            kws = []
            if not is_book_xml and "KeywordList" in container:
                for kwl in container["KeywordList"]:
                    kws.extend([kw for kw in kwl if isinstance(kw, str)])
            article_data["AuthorKeywords"] = "; ".join(kws)

            # 7.9  Publication-type list (correct path)
            pt_list = (container.get("PublicationTypeList", []) if is_book_xml
                       else container["Article"].get("PublicationTypeList", []))
            article_data["PublicationTypes"] = "; ".join(map(str, pt_list))

            # 7.10  MeSH terms
            mesh_terms = []
            if not is_book_xml and "MeshHeadingList" in container:
                for mh in container["MeshHeadingList"]:
                    desc = str(mh["DescriptorName"])
                    quals = [str(q) for q in mh.get("QualifierName", [])]
                    mesh_terms.extend([f"{desc}/{q}" for q in quals]) if quals else mesh_terms.append(desc)
            article_data["MeSH_Terms"] = "; ".join(mesh_terms)

            # 7.11  OA / licence (Unpaywall)
            oa, oa_type, licence, pdf_url = get_unpaywall_info(article_data["DOI"], UNPAYWALL_EMAIL)
            article_data.update({
                "AccessStatus": oa or "Unknown",
                "OA_Type": oa_type or "",
                "License": licence or "",
                "FreePDF_Link": pdf_url or ""
            })
            if licence and licence.lower() not in {"cc0", "cc-by", "cc-by-sa", "cc-by-nc-nd"}: #avoid "cc-by-nc", "cc-by-nd"
                logger.info(f"Skipping PMID {article_data['PMID']} due to licence {licence}")
                continue

            # ---------- 7.12  Full-text from PMC ----------
            full_text = ""
            
            pmc_link = article_data.get("PMC_Link")            # e.g. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9183209/
            if pmc_link:
                m = re.search(r"PMC(\d+)", pmc_link)           # capture digits only → 9183209
                if m:
                    pmcid = m.group(1)
                    try:
                        pmc_handle = Entrez.efetch(
                            db="pmc",
                            id=pmcid,                          # numeric UID (no “PMC” prefix)
                            rettype="full",
                            retmode="xml"
                        )
                        pmc_xml = pmc_handle.read()
                        pmc_handle.close()
            
                        # ----- parser: try fast lxml first, fallback to built-in -----
                        try:
                            soup = BeautifulSoup(pmc_xml, "lxml-xml")
                        except Exception:
                            soup = BeautifulSoup(pmc_xml, "xml")   # built-in xml parser
            
                        body = soup.find("body")
                        if body:
                            full_text = "\n\n".join(
                                p.get_text(" ", strip=True) for p in body.find_all("p")
                            )
                        else:
                            logger.warning(f"No <body> tag in PMC XML for {pmcid}")
            
                    except Exception as exc:
                        logger.warning(f"Full-text fetch failed for PMC{pmcid}: {exc}")
            
            # store result (empty string if nothing retrieved)
            article_data["FullText"] = full_text

            # 7.13  Article-type filter labels
            article_data["ArticleTypeFilters"] = detect_filters(
                [str(pt) for pt in pt_list],
                article_data["Title"],
                article_data["Abstract"],
                is_book_xml=is_book_xml
            )

            data.append(article_data)

        except Exception as exc:
            logger.error(f"Parsing error PMID {article_data.get('PMID','?')}: {exc}")
            continue

# -----------------------------------------------------
# 8.  DataFrame + CSV
# -----------------------------------------------------
df = pd.DataFrame(data)
column_order = [
    "PMID", "Title", "Abstract", "Journal", "PubDate",
    "PublicationYear", "PublicationMonth", "PublicationDay",
    "Authors", "AuthorKeywords", "DOI", "DOI_Link", "PMC_Link",
    "Volume", "Issue", "Pagination", "ELocationID",
    "Language", "Country", "PublicationTypes", "ArticleTypeFilters",
    "MeSH_Terms", "AccessStatus", "OA_Type", "License",
    "FreePDF_Link", "FullText"           # new column
]
df = df[[c for c in column_order if c in df.columns]]
logger.info(f"Processed {len(df)} records – writing CSV")
df.to_csv("pubmed_articles_Metadata_with_license.csv", index=False)
