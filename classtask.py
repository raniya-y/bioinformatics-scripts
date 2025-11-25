# Import required libraries 
from Bio import Entrez
import pandas as pd
from tqdm import tqdm
import time

# Set your email as required by NCBI
Entrez.email = "ryaqub05@gmail.com" 

# Search for all complete M2 gene sequences from Influenza A virus
query = '"Influenza A virus"[Organism] AND "M2"[Gene] AND "complete cds"'
handle = Entrez.esearch(db="nucleotide", term=query, retmax=100000)  # up to 100k records
record = Entrez.read(handle)
handle.close()

# Get list of IDs
ids = record["IdList"]
print(f"Total records found: {len(ids)}")

# Fetch summaries in batches
batch_size = 80
results = []

for start in tqdm(range(0, len(ids), batch_size), desc="Fetching summaries"):
    end = min(start + batch_size, len(ids))
    batch_ids = ids[start:end]
    handle = Entrez.esummary(db="nucleotide", id=",".join(batch_ids))
    summaries = Entrez.read(handle)
    handle.close()
    time.sleep(0.3)  # polite delay to avoid overloading NCBI servers

    for item in summaries:
        accession = item.get("AccessionVersion", "")
        organism = item.get("Organism", "")
        title = item.get("Title", "")
        results.append({
            "Accession": accession,
            "Organism": organism,
            "Title": title
        })

# Convert results to DataFrame
df = pd.DataFrame(results, columns=["Accession", "Organism", "Title"])

# Save to CSV
df.to_csv("influenzaA_M2_accessions.csv", index=False)
print("Saved results to influenzaA_M2_accessions.csv")
