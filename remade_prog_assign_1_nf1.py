# -*- coding: utf-8 -*-
# NAME: cs223_prog_exercise_2.py

"""
AUTHOR: <your name>

  ============== VARIABLE, FUNCTION, etc. NAMING CONVENTIONS ==================
<ALL CAPITOL LETTERS>:  Indicates a symbol defined by a
        #define statement or a Macro.

   <Capitalized Word>:  Indicates a user defined global var, fun, or typedef.

   <all small letters>:  A variable or built in functions.


========================== MODIFICATION HISTORY ==============================
The format for each modification entry is:

MM/DD/YY:
    MOD:     <a description of what was done>
    AUTHOR:  <who made the mod>
    COMMENT: <any special notes to make about the mod (e.g., what other
              modules/code depends on the mod) >

    Each entry is separated by the following line:

====================== END OF MODIFICATION HISTORY ============================
"""

# IMPORTS
import os
import re
from copy import copy
from Bio import Entrez
from Bio import SeqIO
import xmltodict

# ENTREZ LIST OF DBs
""""
{'DbList': ['pubmed', 'protein', 'nuccore', 'ipg', 'nucleotide', 'structure', 'genome',
            'annotinfo', 'assembly', 'bioproject', 'biosample', 'blastdbinfo', 'books',
            'cdd', 'clinvar', 'gap', 'gapplus', 'grasp', 'dbvar', 'gene', 'gds', 'geoprofiles',
            'homologene', 'medgen', 'mesh', 'ncbisearch', 'nlmcatalog', 'omim', 'orgtrack', 'pmc',
            'popset', 'proteinclusters', 'pcassay', 'protfam', 'pccompound', 'pcsubstance', 
            'seqannot', 'snp', 'sra', 'taxonomy', 'biocollections', 'gtr']}
"""
# HELPER FUNS
def get_valid_ids(list_of_ids, db_name, return_type, query_type):
    """Validate each Entrez ID# that is in the list_of_ids  argument. Invalid
    ids are removed from the list. A valid id# is one that does not generate an error
    for the given DB, return type, and query type."""

    # Return if there are no ids to check
    if len(list_of_ids) <= 0:
        return list_of_ids

    # Iterate through list and check each id
    copy_list_of_ids = copy(list_of_ids)    # Make a copy of the list of ids that will be modified.
    for lid in list_of_ids:
        if query_type == 'efetch':
            try:
                Entrez.efetch(db = db_name, id = int(lid), rettype = return_type)
            except:
                print("get_valid_ids: The Entrez.efetch id# " + str(lid) + " for NCBI " + str(db_name) +
                      " DB, and return type = " + str(return_type) + " is not valid.")
                copy_list_of_ids.remove(lid)

    # When we reach here, copy_list_of_ids will only have valid ids
    return copy_list_of_ids


# MAIN FUNCTION
def main():
    start = 31094927
    end = 31377677
    first_SNP_location = 31357033
    second_SNP_location = 31095311
    for record in SeqIO.parse("gene1.fasta", "fasta"):
        print(record.id)
        # print(record.seq)
        # print(len(record.seq))
    # NC_000017.11:31094927-31377677
    # first SNP location C>T 31357033
    # second SNP location 31095311 T>A
    # with open('gene.fna') as f:
    #     # NF1 gene
    #     lines = f.readlines()
    #     lines = [line.replace("\n","").strip() for line in lines]
    #     lines = [line.replace(" ","").strip() for line in lines]
    #     gene = "".join(lines)
    #     print(len(gene))
    # """The main program collects """

    # Init Entrez email addr to use
    # Entrez.email = "David.Warshawsky@sjsu.edu"

    # Extract Entrez info for the specified ACC number
    # print("\n\n#####")
    # acc_num = "NC_000017.11"
    # print("PROCESSING ACC #: ", acc_num)

    # Entrez.esearch for info related to the current accession number
    # handle = Entrez.esearch(db='gene', term='human[organism] ' + str(acc_num))

    # Entrez.read the result returned from NCBI search for the current accession number
    # handle_read = Entrez.read(handle)

    # Get Idlist from previous read of accession number info
    # handle_read_ids = handle_read['IdList']
    # print("handle_read_ids: ", handle_read_ids)

    # Make sure IDs in Idlist are valid. Return a list of only valid IDs.
    # valid_handle_read_ids = get_valid_ids(handle_read_ids, 'nucleotide', 'fasta', 'efetch')
    # if valid_handle_read_ids == []:
    #     print("   epigen_pipeline_get_entrez_info: No valid IDs found for ACC #" + str(acc_num) + ".")
    #     return
    #
    # # If we get here, valid_handle_read_ids has usable ID nums.
    # for acc_handle_id in valid_handle_read_ids:
    #     # Entrez.efetch based on ID#
    #     handle_read_id_info = Entrez.efetch(db="nucleotide", id=int(acc_handle_id), rettype="fasta")
    #
    #     # SeqIO.read  fasta info
    #     record_seq = SeqIO.read(handle_read_id_info, "fasta")
    #
    #     # Print the fasta sequences
    #     print()
    #     # Print sequence ACC number and sequence
    #     print("record_seq.id: ", record_seq.id)
    #     print("record_seq.seq: ", record_seq.seq)
    #     print(len(record_seq.seq))
    #     sum += len(record_seq.seq)
    #     print("==========")


    print()
    print("sum length: " + str(sum))
    print ("############### DONE ###############")
    return










#####################################################################################
if __name__ == "__main__":
    main()
else:
    print("cs223_prog_assign_2.py : Is intended to be executed and not imported.")


