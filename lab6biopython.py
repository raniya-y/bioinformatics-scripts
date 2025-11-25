
from Bio import Entrez, SeqIO
Entrez.email = "ryaqub05@gmail.com"

# search for BRCA1 gene in humans
search = Entrez.esearch(db="nucleotide", term="BRCA1 Homo sapiens", retmax=1)
result = Entrez.read(search)
search.close()

# get the first accession ID
accession_id =result["IdList"][0]

# fetch the FASTA sequence
fetch = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
record = SeqIO.read(fetch, "fasta")
fetch.close()

# print results
print("Accession Number: ", record.id)
print("Sequence Length: ", len(record.seq)) 
print("FASTA Sequence:\n", record.seq)