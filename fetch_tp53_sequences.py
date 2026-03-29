# TP53 sequence retrieval and BLAST analysis
# used this for my MSc project at IBAB
# reference: Biopython tutorial + JHU genomic data science course

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import time

Entrez.email = "aashrithasr@gmail.com"

# fetching wild type TP53 protein sequence from NCBI
# accession P04637 is the human TP53 uniprot entry
print("fetching TP53 wild type sequence...")

handle = Entrez.efetch(db="protein", id="P04637", rettype="fasta", retmode="text")
wt_record = SeqIO.read(handle, "fasta")
handle.close()

print("sequence fetched:", wt_record.id)
print("length:", len(wt_record.seq))

# save to fasta file
with open("tp53_wt.fasta", "w") as f:
    SeqIO.write(wt_record, f, "fasta")

print("saved to tp53_wt.fasta")

# these are the 10 mutant variants i looked at from ClinVar
# pulled accessions manually from NCBI protein database
variants = ["R175H", "R248W", "R248Q", "R273H", "R273C",
            "R282W", "G245S", "V143A", "R249S", "Y220C"]

print("\nvariants being analysed:")
for v in variants:
    print(" -", v)

# running blastp on the wild type sequence
# note: this takes several minutes, don't close terminal
print("\nrunning BLASTp on wild type TP53...")
print("(this may take 5-10 minutes)")

result_handle = NCBIWWW.qblast("blastp", "nr", wt_record.seq)

with open("blast_results.xml", "w") as blast_out:
    blast_out.write(result_handle.read())

print("blast complete. saved to blast_results.xml")

# parsing results
print("\nparsing top 10 hits...")

with open("blast_results.xml") as f:
    blast_records = NCBIXML.parse(f)
    blast_record = next(blast_records)

for i, alignment in enumerate(blast_record.alignments[:10]):
    hsp = alignment.hsps[0]
    identity = round((hsp.identities / hsp.align_length) * 100, 1)
    print(f"{i+1}. {alignment.title[:55]}")
    print(f"   score={hsp.score}  evalue={hsp.expect}  identity={identity}%")
