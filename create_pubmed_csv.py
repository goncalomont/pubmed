from Bio import Entrez, Medline
import time
import xml.etree.ElementTree as ET
import csv
import requests
from bs4 import BeautifulSoup

# Set your email for NCBI Entrez
Entrez.email = "apitest@example.com"  # Replace with a valid email

def search_pubmed_and_fetch_details(search_term, max_articles=20):
    """
    Searches PubMed for a given term, fetches article details (PMID, Title, Abstract)
    for the specified maximum number of articles.
    Returns a list of dictionaries: [{'PMID': pmid, 'Title': title, 'Abstract': abstract}, ...]
    """
    articles = []
    print(f"Searching PubMed for '{search_term}' (max {max_articles} articles)...")
    try:
        # Use Entrez.esearch to get PMIDs
        handle = Entrez.esearch(db="pubmed", term=search_term, retmax=str(max_articles))
        search_results = Entrez.read(handle)
        handle.close()
        pmids = search_results["IdList"]

        if not pmids:
            print("No articles found for the search term.")
            return articles

        print(f"Found {len(pmids)} PMIDs. Fetching details...")

        for i, pmid in enumerate(pmids):
            print(f"  Fetching details for PMID {pmid} ({i+1}/{len(pmids)})...")
            try:
                # Use Entrez.efetch to get article details in Medline format
                fetch_handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
                medline_text = fetch_handle.read()
                fetch_handle.close()

                # Parse Medline text
                # Medline.parse returns an iterator, so we expect one record
                medline_records = Medline.parse(medline_text)
                record = next(medline_records, None)

                if record:
                    title = record.get("TI", "No Title Available")
                    abstract = record.get("AB", "No Abstract Available")
                    articles.append({"PMID": pmid, "Title": title, "Abstract": abstract})
                else:
                    print(f"    Warning: Could not parse Medline record for PMID {pmid}")
                    articles.append({"PMID": pmid, "Title": "Error parsing Medline", "Abstract": "Error parsing Medline"})

                time.sleep(0.34) # NCBI API rate limit (3 requests per second without API key)

            except Exception as e_fetch:
                print(f"    Error fetching or parsing details for PMID {pmid}: {e_fetch}")
                # Optionally add a placeholder or skip
                articles.append({"PMID": pmid, "Title": "Error fetching details", "Abstract": str(e_fetch)})
                time.sleep(0.34) # Still sleep to avoid overwhelming the server after an error

        print(f"Successfully fetched details for {len(articles)} articles.")

    except Exception as e_search:
        print(f"Error during PubMed search or initial fetch: {e_search}")

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

def parse_license_string(license_text):
    """
    Parses a license string or URL to a short license identifier.
    e.g., "https://creativecommons.org/licenses/by/4.0/" -> "cc-by"
    """
    if not license_text:
        return None
    license_text = license_text.lower()
    if "creativecommons.org/licenses/by-nc-nd/" in license_text:
        return "cc-by-nc-nd"
    elif "creativecommons.org/licenses/by-nc-sa/" in license_text:
        return "cc-by-nc-sa"
    elif "creativecommons.org/licenses/by-nd/" in license_text:
        return "cc-by-nd"
    elif "creativecommons.org/licenses/by-sa/" in license_text:
        return "cc-by-sa"
    elif "creativecommons.org/licenses/by-nc/" in license_text:
        return "cc-by-nc"
    elif "creativecommons.org/licenses/by/" in license_text: # Must be after more specific "by-" versions
        return "cc-by"
    elif "creativecommons.org/publicdomain/zero/" in license_text or "creativecommons.org/publicdomain/mark/" in license_text:
        return "cc0"
    # Add more specific parsing rules if needed
    # For now, return a generic part if it's a URL or a simplified string
    if "http" in license_text:
        try:
            # Try to get path part
            path = license_text.split("://")[1].split("/")
            if len(path) > 2: # e.g. creativecommons.org/licenses/by/4.0 -> by
                return path[2]
        except IndexError:
            pass # Fall through

    # Fallback for non-URL or unparsed URLs
    # Remove common terms and simplify
    simplified = license_text.replace("license", "").replace("public", "").replace("domain","").replace("creative commons","").strip()
    simplified = "".join(c for c in simplified if c.isalnum() or c == '-').strip('-')
    return simplified if simplified else None

# Helper function to check full text availability
def is_full_text_available(full_text_content):
    """
    Checks if the full_text_content is non-empty and not a known placeholder string.
    Returns True if text is considered available, False otherwise.
    """
    if not full_text_content or full_text_content.strip() == "":
        return False
    placeholder_strings = [
        "[No PMCID link for full text]",
        "[Fetching timed out]",
        "[Metadata fetch failed]",
        "[No content extracted or not applicable]",
        "[Unsupported content type", # Check for prefix
        "[PDF content not extracted]",
        "[HTTP error", # Check for prefix
        "[Request error]",
        "[Error during processing]",
        "[XML content not easily parsable to plain text]"
    ]
    for placeholder in placeholder_strings:
        if placeholder.startswith("[") and placeholder.endswith("]") and full_text_content == placeholder:
            return False
        if full_text_content.startswith(placeholder): # For placeholders like "[Unsupported content type: ...]"
            return False
    return True

# Helper function to check abstract availability
def is_abstract_available(abstract_content):
    """
    Checks if the abstract_content is non-empty and not the 'No Abstract Available' placeholder.
    Returns True if abstract is considered available, False otherwise.
    """
    if not abstract_content or abstract_content.strip() == "" or abstract_content == "No Abstract Available":
        return False
    return True

def extract_pmcid_and_link(xml_data, pmid):
    """
    Tries to extract PMCID, a potential full-text link, and license information from the XML metadata.
    Returns: (pmcid, full_text_link, license_from_xml)
    """
    pmcid = None
    full_text_link = None
    license_from_xml = None

    if not xml_data:
        return pmcid, full_text_link, license_from_xml

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

    # Try to find license information
    try:
        if root is not None: # Ensure root was parsed
            # Look for <license> tag often under <permissions>
            permissions_node = root.find(".//permissions")
            if permissions_node is not None:
                license_node = permissions_node.find(".//license")
                if license_node is not None:
                    # Check for <license_ref> first (often a URL)
                    license_ref_node = license_node.find(".//license_ref")
                    if license_ref_node is not None and license_ref_node.text:
                        license_from_xml = parse_license_string(license_ref_node.text)
                    # If no <license_ref>, check text content of <license> itself
                    elif license_node.text and license_node.text.strip():
                         license_from_xml = parse_license_string(license_node.text.strip())

                # Fallback: Check for license_ref directly under permissions if no license tag
                if not license_from_xml:
                    license_ref_node = permissions_node.find(".//license_ref")
                    if license_ref_node is not None and license_ref_node.text:
                        license_from_xml = parse_license_string(license_ref_node.text)

            # Broader search if not found under <permissions>
            if not license_from_xml:
                license_node = root.find(".//license") # More general search
                if license_node is not None:
                    license_ref_node = license_node.find(".//license_ref")
                    if license_ref_node is not None and license_ref_node.text:
                        license_from_xml = parse_license_string(license_ref_node.text)
                    elif license_node.text and license_node.text.strip():
                        license_from_xml = parse_license_string(license_node.text.strip())

            if not license_from_xml:
                # Some articles might have license info directly in <article-meta>
                article_meta = root.find(".//article-meta")
                if article_meta is not None:
                    permissions = article_meta.find("permissions")
                    if permissions is not None:
                        license_p = permissions.find("license/p") # e.g. JATS format <license><p>...</p></license>
                        if license_p is not None and license_p.text:
                             license_from_xml = parse_license_string(license_p.text)
                        if not license_from_xml: # Check for <license_ref> inside <license> within <permissions>
                            license_ref_in_permissions = permissions.find("license/license_ref")
                            if license_ref_in_permissions is not None and license_ref_in_permissions.text:
                                license_from_xml = parse_license_string(license_ref_in_permissions.text)


            if license_from_xml:
                print(f"  License from XML (parsed): {license_from_xml}")

    except Exception as e:
        # Don't let license extraction errors stop PMCID/link extraction
        print(f"Error extracting license from XML for PMID {pmid}: {e}")


    return pmcid, full_text_link, license_from_xml

def extract_license_from_html(soup, pmid):
    """
    Extracts license information from HTML soup by looking for Creative Commons links.
    Returns a short license string (e.g., "cc-by", "cc0") or None.
    """
    if not soup:
        return None

    license_found = None
    try:
        # Common Creative Commons patterns
        # Example: <a href="http://creativecommons.org/licenses/by/4.0/">
        # Example: <a rel="license" href="https://creativecommons.org/publicdomain/zero/1.0/">

        # Prioritize links with rel="license"
        license_links = soup.find_all('a', rel='license')
        if not license_links:
            # Fallback: search all 'a' tags that contain 'creativecommons.org' in their href
            license_links = soup.find_all('a', href=lambda href: href and 'creativecommons.org' in href.lower())

        for link in license_links: # This loop will now only process pre-filtered CC links or rel=license links
            href = link.get('href', '').lower()
            # Ensure we only parse actual Creative Commons links here, or rel="license"
            if 'creativecommons.org' in href or link.get('rel') == ['license']:
                parsed_license = parse_license_string(href)
                if parsed_license:
                    # Take the first valid one found
                    license_found = parsed_license
                    print(f"    License from HTML (PMID {pmid}): {href} -> {license_found}")
                    break

        # If no specific CC link found through href, check for common license text in footer or license sections
        # This is a more heuristic approach and might need refinement
        if not license_found:
            # Look for elements that might contain license text
            # Common class names: 'license', 'footer-license', 'copyright-license'
            # Common tag types: div, p, span, footer
            possible_license_elements = soup.find_all(['div', 'p', 'span', 'footer'],
                                                      class_=['license', 'footer-license', 'copyright-license', 'licenses'])

            # Also check for elements with id containing 'license'
            for id_val in ['license', 'licenses', 'copyright-license']:
                el = soup.find(id=lambda x: x and id_val in x.lower())
                if el:
                    possible_license_elements.append(el)

            for element in possible_license_elements:
                text_content = element.get_text(separator=" ").lower()
                # Simple check for common license phrases if not a URL
                if "cc-by-nc-nd" in text_content or "creative commons attribution-noncommercial-noderivatives" in text_content:
                    license_found = "cc-by-nc-nd"
                    break
                elif "cc-by-nc-sa" in text_content or "creative commons attribution-noncommercial-sharealike" in text_content:
                    license_found = "cc-by-nc-sa"
                    break
                elif "cc-by-nd" in text_content or "creative commons attribution-noderivatives" in text_content:
                    license_found = "cc-by-nd"
                    break
                elif "cc-by-sa" in text_content or "creative commons attribution-sharealike" in text_content:
                    license_found = "cc-by-sa"
                    break
                elif "cc-by-nc" in text_content or "creative commons attribution-noncommercial" in text_content: # must be after -nd and -sa
                    license_found = "cc-by-nc"
                    break
                elif "cc-by" in text_content or "creative commons attribution" in text_content : # must be after other by- variants
                    license_found = "cc-by"
                    break
                elif "cc0" in text_content or "public domain zero" in text_content or "public domain mark" in text_content:
                    license_found = "cc0"
                    break
            if license_found and not any(href_link in str(element) for href_link in ["creativecommons.org", "license"]): # Avoid re-logging if found via text but was in a link
                 print(f"    License from HTML text (PMID {pmid}): {license_found}")


    except Exception as e:
        print(f"    Error extracting license from HTML for PMID {pmid}: {e}")

    return license_found

def main():
    ALLOWED_LICENSES = ["cc0", "cc-by", "cc-by-sa", "cc-by-nc-nd"] # Define allowed licenses

    # Call the new function to search PubMed and fetch initial details
    search_term = "open access genomics AND human" # Example search term
    max_results = 20 # Fetch up to 20 articles

    # Replace parse_pubmed_documents with search_pubmed_and_fetch_details
    parsed_articles = search_pubmed_and_fetch_details(search_term, max_articles=max_results)

    if not parsed_articles:
        print(f"No articles fetched from PubMed for search term '{search_term}'. Exiting.")
        return

    print(f"Fetched {len(parsed_articles)} articles from PubMed for search term '{search_term}'. Now processing for further metadata...\n")

    # Process all articles
    articles_to_process = parsed_articles
    # print(f"Processing the first {len(articles_to_process)} articles for this test run.\n") # Comment out or remove test line

    processed_articles_data = [] # To store data for CSV writing
    filtered_out_count = 0 # Counter for filtered articles

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

        article_license = None # Initialize article_license for each article
        if xml_data:
            pmcid, full_text_link, license_from_xml = extract_pmcid_and_link(xml_data, pmid)
            article_license = license_from_xml # XML takes precedence

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
                            html_soup = BeautifulSoup(response.content, 'lxml') # Renamed to html_soup to avoid conflict

                            # Try to extract license from HTML if not found in XML
                            if not article_license and html_soup:
                                license_from_html = extract_license_from_html(html_soup, pmid)
                                if license_from_html:
                                    article_license = license_from_html

                            article_body = html_soup.find('div', class_='article-content') # Common on PMC
                            if not article_body:
                                article_body = html_soup.find('div', class_='rendered_body') # Another common one
                            if not article_body:
                                article_body = html_soup.find('body')

                            if article_body:
                                full_text_content = article_body.get_text(separator='\n\n', strip=True)
                                full_text_content = full_text_content[:5000] # Limit length
                                print(f"    Successfully extracted ~{len(full_text_content)} chars of text content.")
                            else:
                                print("    Could not find main article content body/div in HTML.")
                                # Still try to get license even if main body not found
                                if not article_license and html_soup: # Check again, in case html_soup was valid but body wasn't
                                    license_from_html = extract_license_from_html(html_soup, pmid)
                                    if license_from_html:
                                        article_license = license_from_html


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

        is_open_access = True if pmcid and full_text_link else False # Keep this as is

        final_license_value = article_license # Initialize with the detected license

        # Start of new filtering logic based on license and content availability
        # License Handling:
        # The script checks if a license is found. If not, it's marked 'unavailable'.
        # If a license is found, it's checked against ALLOWED_LICENSES.
        # Articles with unallowed licenses are filtered out.
        if not article_license or article_license.strip() == "":
            # Case: License is missing or empty. Mark as 'unavailable'.
            # The article is not immediately excluded at this point.
            # It can still be included if it has available content (full text or abstract).
            final_license_value = "unavailable"
            print(f"  Article {pmid}: License not found or empty, marked as 'unavailable'.")
        elif article_license not in ALLOWED_LICENSES:
            # Case: License is present but not in the allowed list. Exclude the article.
            filtered_out_count += 1
            print(f"  Article {pmid} filtered out. Unallowed License: '{article_license}'.")
            print("-" * 20)
            time.sleep(0.5) # Respect API limits even when skipping
            continue # Exclude article, proceed to the next one.

        # Content Availability Check:
        # After license check, the article's content (full text and abstract) is assessed.
        # The article must have either available full text or an available abstract to be included.
        full_text_is_avail = is_full_text_available(full_text_content)
        abstract_is_avail = is_abstract_available(abstract)

        if not full_text_is_avail and not abstract_is_avail:
            # Case: Neither full text nor abstract is available. Exclude the article.
            filtered_out_count += 1
            reason = "No full text and no abstract available"
            if full_text_content and full_text_content not in ["[No PMCID link for full text]", "[Metadata fetch failed]"]:
                reason = f"Full text placeholder '{full_text_content}' and no abstract"
            elif abstract == "No Abstract Available" and not full_text_is_avail :
                 reason = "No abstract available and no usable full text"
            print(f"  Article {pmid} filtered out. Reason: {reason}.")
            print("-" * 20)
            time.sleep(0.5) # Respect API limits even when skipping
            continue # Exclude article, proceed to the next one.

        # Article Inclusion Criteria:
        # An article is included if it meets the following conditions:
        # 1. Its license is either in ALLOWED_LICENSES OR was marked as 'unavailable' (originally missing/empty).
        # AND
        # 2. It has available full text OR an available abstract.
        processed_articles_data.append({
            "PMID": pmid,
            "Title": title,
            "Abstract": abstract, # Store original abstract
            "IsOpenAccess": is_open_access,
            "FullTextContent": full_text_content if full_text_content else "[No content extracted or not applicable]",
            "License": final_license_value # Use the final_license_value
        })
        print(f"  Article {pmid} (License: '{final_license_value}') added to dataset. Full text available: {full_text_is_avail}, Abstract available: {abstract_is_avail}.")
        print("-" * 20)

        # Respect NCBI API usage guidelines
        time.sleep(0.5)

    # Write to CSV
    csv_file_path = "pubmed_articles.csv"
    csv_header = ["PMID", "Title", "Abstract", "IsOpenAccess", "FullTextContent", "License"] # Add "License" to header

    try:
        with open(csv_file_path, "w", newline="", encoding="utf-8") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_header)
            writer.writeheader()
            for article_data in processed_articles_data:
                writer.writerow(article_data)
        print(f"\nSuccessfully wrote {len(processed_articles_data)} articles to {csv_file_path}")
        if filtered_out_count > 0:
            # Update print statement for more general filtering
            print(f"Filtered out {filtered_out_count} articles due to license or content availability criteria.")
    except IOError:
        print(f"Error: Could not write to CSV file {csv_file_path}")
    except Exception as e:
        print(f"An unexpected error occurred during CSV writing: {e}")

if __name__ == "__main__":
    main()
