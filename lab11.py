import requests
from Bio.Seq import Seq

server = "https://rest.ensembl.org"
ext = "/sequence/id/ENSG00000012048?content-type=text/plain"
response = requests.get(server + ext, headers={"Content-Type": "text/plain"})

if not response.ok:
    response.raise_for_status()

dna_seq = Seq(response.text)
print("Length:", len(dna_seq))
print("First 100 bases:", dna_seq[:100])
