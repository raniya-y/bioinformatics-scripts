# LAB TASK: Gene Variant Analysis (BRCA1)

import requests

def get_gene_info(gene_name="BRCA1"):
    server = "https://rest.ensembl.org"
    ext = f"/lookup/symbol/homo_sapiens/{gene_name}?expand=1"
    r = requests.get(server+ext, headers={"Content-Type": "application/json"})
 
    if not r.ok:
       r.raise_for_status()

    return r.json()

def get_variants(gene_id):
    server = "https://rest.ensembl.org"
    ext = f"/overlap/id/{gene_id}?feature=variation"
    r = requests.get(server+ext, headers={"Content-Type": "application/json"})
  
    if not r.ok:
        r.raise_for_status()
    return r.json()

coding_terms = {
	    "missense_variant", "synonymous_variant", "stop_gained", "stop_lost",
	    "start_lost", "frameshift_variant", "inframe_insertion", "inframe_deletion",
	    "splice_acceptor_variant", "splice_donor_variant", "protein_altering_variant",
	    "coding_sequence_variant"
	}
	
noncoding_terms = {
	    "intron_variant", "upstream_gene_variant", "downstream_gene_variant",
	    "5_prime_UTR_variant", "3_prime_UTR_variant", "regulatory_region_variant",
	    "TF_binding_site_variant", "non_coding_transcript_exon_variant",
	    "nc_transcript_variant", "intergenic_variant"
	}

gene_info = get_gene_info("BRCA1")
print("Gene:", gene_info["display_name"])
print("Location:", gene_info["seq_region_name"], gene_info["start"], "-", gene_info["end"])

variants = get_variants(gene_info["id"])
print(f"Found {len(variants)} variants in BRCA1")

coding, noncoding, others = [], [], []

for v in variants[:900]:  # limit to first 900 for readability
    cons = v.get("consequence_type", "N/A")
    if cons in coding_terms:
        coding.append((v["id"], cons))
    elif cons in noncoding_terms:
        noncoding.append((v["id"], cons))
    else:
        others.append((v["id"], cons))
        
print(f"\n=== Coding Consequences Detected ({len(coding)}) ===")
if coding:
    for vid, cons in coding:
        print("Variant:", vid, "| Consequence:", cons)
else:
    print("None detected")

print(f"\n=== Non-Coding Consequences Detected ({len(noncoding)}) ===")
if noncoding:
    for vid, cons in noncoding:
        print("Variant:", vid, "| Consequence:", cons)
else:
    print("None detected")

print(f"\n=== Other/Unclassified Consequences ({len(others)}) ===")
if others:
    for vid, cons in others:
        print("Variant:", vid, "| Consequence:", cons)
else:
    print("None detected")