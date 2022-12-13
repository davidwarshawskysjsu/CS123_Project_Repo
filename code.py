import Bio
import Bio
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "david.warshawsky@sjsu.edu"
info = Entrez.esearch(db = "gene",term = "NF1")
record = Entrez.read(info) 
id_list = record['IdList']
handle_read_id_info = Entrez.efetch(db="gene",id=int(id_list[0]),rettype="fasta"
)
record_seq = SeqIO.read(handle_read_id_info, "fasta")
stringIO = handle_read_id_info.read().decode("utf-8")

