#!/usr/bin/env python3

import requests
from requests.adapters import HTTPAdapter, Retry
import json


retries = Retry(total=5, backoff_factor=1, status_forcelist=[500,502,503,504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


import re
def get_next_link(headers):
    if "Link" in headers:
        # next link regex
        re_next_link = re.compile(r'<(.+)>; rel="next"')
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)
    return None

def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url, headers={"Accept-Encoding": "identity"})
        response.raise_for_status()
        total = response.headers.get("x-total-results", "unknown")
        yield response, total
        batch_url = get_next_link(response.headers)


def filter_entry(entry):
    for feature in entry.get("features", []):
        if feature["type"] == "Signal":
            start = feature["location"]["start"].get("value")
            end = feature["location"]["end"].get("value")
            description = feature.get("description","")
            if not start or not end: # For null and None
                continue
            sp_length = int(end) - int(start) + 1
            if sp_length > 13 and not description:
                return True
    return False


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

    # Signal peptide cleavage site
    cleavage = ""
    for f in entry.get("features", []):
        if f["type"] == "Signal":
            start = f["location"]["start"].get("value")
            end = f["location"]["end"].get("value")
            if start is not None and end is not None:
                cleavage = int(end)
            #cleavage = f["location"]["end"]["value"]
            break

    return (accession, organism, kingdom, length, cleavage)


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
        # TSV header
        print("Accession\tOrganism\tKingdom\tLength\tCleavage_site", file=ofs)
        for entry in filtered_json:
            fields = extract_function(entry)
            print(*fields, sep="\t", file=ofs)
    with open(fasta_file,"w",encoding="utf-8") as ofs:
        for acc, seq in sequences:
            ofs.write(f">{acc}\n{seq}\n")


#batch_size = 500
url = "https://rest.uniprot.org/uniprotkb/search?format=json&query=%28%28taxonomy_id%3A2759%29+AND+%28reviewed%3Atrue%29+AND+%28existence%3A1%29+AND+%28fragment%3Afalse%29+AND+%28ft_signal_exp%3A*%29+AND+%28length%3A%5B40+TO+*%5D%29%29&size=500"

output_file = "positive.tsv"
fasta_file = "positive.fasta"
get_dataset(url, filter_entry, extract_fields, output_file)
