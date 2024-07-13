#!/usr/bin/env python3

import requests
import pandas as pd
import re
import argparse

def get_gene_info(ensembl_id):
    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{ensembl_id}?content-type=application/json"
    
    response = requests.get(server + ext, headers={"Content-Type": "application/json"})
    if not response.ok:
        response.raise_for_status()
    
    gene_info = response.json()
    return gene_info

def get_xrefs(ensembl_id):
    server = "https://rest.ensembl.org"
    ext = f"/xrefs/id/{ensembl_id}?content-type=application/json"
    
    response = requests.get(server + ext, headers={"Content-Type": "application/json"})
    if not response.ok:
        response.raise_for_status()
    
    xrefs = response.json()
    return xrefs

def get_entrez_id(xrefs):
    for xref in xrefs:
        if xref['dbname'] == 'EntrezGene':
            return xref['primary_id']
    return None

def get_gene_go_terms(ncbi_id):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id={ncbi_id}&rettype=xml"
    response = requests.get(url)
    
    if response.status_code != 200:
        print(f"Error: Unable to fetch data (status code: {response.status_code})")
        return []
    
    xml_content = response.content.decode()
    
    go_terms = re.findall(r'<Dbtag_db>GO</Dbtag_db>\s*<Dbtag_tag>\s*<Object-id>\s*<Object-id_id>(\d+)</Object-id_id>', xml_content)
    go_terms = [term.zfill(7) for term in go_terms]
    
    return go_terms

def get_kegg_pathways(ncbi_id):
    url = f"http://rest.kegg.jp/link/pathway/oas:{ncbi_id}"
    response = requests.get(url)
    response.raise_for_status()
    

    pathways = [line.split('\t')[1] for line in response.text.strip().split('\n') if line]
    return list(set(pathways))

def fetch_gene_info(ensembl_ids):
    results = []
    
    for ensembl_id in ensembl_ids:
        if isinstance(ensembl_id, float) or pd.isna(ensembl_id):
            continue
        ensembl_id = ensembl_id.strip()
        xrefs = get_xrefs(ensembl_id)
        ncbi_id = get_entrez_id(xrefs)
        
        if ncbi_id:
            go_terms = get_gene_go_terms(ncbi_id)
            kegg_pathways = get_kegg_pathways(ncbi_id)
            
            if not go_terms:
                go_terms = ['N/A']
            if not kegg_pathways:
                kegg_pathways = ['N/A']
            
            results.append({
                'Ensembl ID': ensembl_id,
                'NCBI ID': ncbi_id,
                'GO Terms': ', '.join(go_terms),
                'KEGG Pathways': ', '.join(kegg_pathways)
            })
        else:
            results.append({
                'Ensembl ID': ensembl_id,
                'NCBI ID': 'Not Found',
                'GO Terms': 'N/A',
                'KEGG Pathways': 'N/A'
            })
    
    return results

def main(input_file, output_file):
    data = pd.read_csv(input_file, delimiter='\t')
    ensembl_ids = data['gene_id'].tolist()
    gene_info = fetch_gene_info(ensembl_ids)
    
    df = pd.DataFrame(gene_info)
    df.to_csv(output_file, index=False)
    # print(df)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch Gene Information using Ensembl IDs")
    parser.add_argument("input_file", type=str, help="Path to the input file containing Ensembl IDs")
    parser.add_argument("output_file", type=str, help="Path to the output CSV file to save gene information")

    args = parser.parse_args()
    main(args.input_file, args.output_file)
