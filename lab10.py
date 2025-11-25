from Bio import Entrez, SeqIO

Entrez.email = "ryaqub05@gmail.com" 

# Search for INS gene sequences in humans
handle = Entrez.esearch(db="nucleotide", term="INS Homo sapiens", retmax=20)
record = Entrez.read(handle)
handle.close()

ids = record["IdList"]
valid_count = 0

for acc_id in ids:
    if valid_count >= 5:
        break

    try:
        fetch = Entrez.efetch(db="nucleotide", id=acc_id, rettype="fasta", retmode="text")
        records = list(SeqIO.parse(fetch, "fasta"))
        fetch.close()

        for seq_record in records:
            if seq_record.id and seq_record.seq and len(seq_record.seq) > 0:
                valid_count += 1
                print(f"\n Sequence {valid_count}")
                print("Accession ID:", seq_record.id)
                print("Length:", len(seq_record.seq))
                print("FASTA Format:\n", seq_record.format("fasta").strip())
                break  # Only print one valid record per ID

    except Exception as e:
        print(f"\n Skipped ID {acc_id}: Error fetching sequence. Reason: {e}")
