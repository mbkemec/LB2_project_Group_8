import requests
from requests.adapters import HTTPAdapter, Retry
import json
import re

retries = Retry(total=5, backoff_factor=2, status_forcelist=[500,502,503,504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def get_next_link(headers):
    if "Link" in headers:
        re_next_link = re.compile(r'<(.+)>; rel="next"')
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)
    return None

def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url, headers={"Accept-Encoding": "identity"}, timeout=30)
        response.raise_for_status()
        total = response.headers.get("x-total-results", "unknown")
        yield response, total
        batch_url = get_next_link(response.headers)

def filter_entry(entry):
    return True
"""
    for feature in entry.get("features", []):
        if feature["type"] == "SIGNAL":
            return False
"""


def extract_fields(entry):
    accession = entry["primaryAccession"]
    organism = entry["organism"]["scientificName"]

    lineage = entry["organism"].get("lineage", [])
    kingdom = "Other"
    if "Metazoa" in lineage:
        kingdom = "Metazoa"
    elif "Fungi" in lineage:
        kingdom = "Fungi"
    elif "Viridiplantae" in lineage or "Plantae" in lineage:
        kingdom = "Plants"

    length = entry["sequence"]["length"]


    tm_first90 = False
    for f in entry.get("features", []):
        if f["type"] in ["Transmembrane"]:
            if f["description"] == "Helical":
                start = f["location"]["start"].get("value")
                if start is not None and start <= 90:
                    tm_first90 = True
                    break

    return (accession, organism, kingdom, length, tm_first90)

def get_dataset(search_url, filter_function, extract_function, output_file_name):
    filtered_json = []
    sequences = []
    n_total, n_filtered = 0, 0
    for batch, total in get_batch(search_url):
        batch_json = json.loads(batch.text)
        for entry in batch_json.get("results", []):
            n_total += 1
            if filter_function(entry):
                n_filtered += 1
                filtered_json.append(entry)
                sequences.append((entry["primaryAccession"],entry["sequence"]["value"]))
    print(f"Total entry: {n_total}, filtered: {n_filtered}")

    with open(output_file_name, "w", encoding="utf-8") as ofs:
        print("Accession\tOrganism\tKingdom\tLength\tTransmembrane_Helix", file=ofs)
        for entry in filtered_json:
            fields = extract_function(entry)
            print(*fields, sep="\t", file=ofs)
    with open(fasta_file,"w",encoding="utf-8") as ofs:
        for acc, seq in sequences:
            ofs.write(f">{acc}\n{seq}\n")


url = "https://rest.uniprot.org/uniprotkb/search?format=json&query=%28%28fragment%3Afalse%29+AND+%28taxonomy_id%3A2759%29+AND+%28length%3A%5B40+TO+*%5D%29+NOT+%28ft_signal%3A*%29+AND+%28%28cc_scl_term_exp%3ASL-0091%29+OR+%28cc_scl_term_exp%3ASL-0191%29+OR+%28cc_scl_term_exp%3ASL-0173%29+OR+%28cc_scl_term_exp%3ASL-0209%29+OR+%28cc_scl_term_exp%3ASL-0204%29+OR+%28cc_scl_term_exp%3ASL-0039%29%29AND+%28reviewed%3Atrue%29+AND+%28existence%3A1%29%29&size=200"

output_file = "negative.tsv"
fasta_file = "negative.fasta"
get_dataset(url, filter_entry, extract_fields, output_file)
