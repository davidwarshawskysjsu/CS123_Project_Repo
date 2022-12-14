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
import copy
from copy import copy

import Bio.SeqRecord
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParam import ProtParamData
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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


def get_nucleotide_in_NF1(location,sequence):
    start = 31094927
    relative_location = location - start
    if relative_location < 0:
        return None
    return sequence[relative_location]

def bio_seq_analytics():
    pass


def complete_nucleotide_codons_N(record):
    seq = record.seq
    seq_list = list(seq)
    result = len(seq_list) / 3
    while not (result).is_integer():
        seq_list.append("N")
        result = len(seq_list)/3
    return SeqRecord(Seq(''.join(seq_list)))

def copy_record(record):
    seq_list = list(record.seq)
    return SeqRecord(Seq(''.join(seq_list)))

def verify_nucleotide(record,location,nucleotide):
    seq = record.seq
    if seq[location] != nucleotide:
        return False
    return True

def change_nucleotide(record,location,nucleotide):
    seq_list = list(record.seq)
    seq_list[location] = nucleotide
    return SeqRecord(Seq(''.join(seq_list)))

def verify_all_nucleotides(record):
    seq_len = len(record.seq)
    count = 0
    count += sum([1 for x in record.seq if (x == "A" or x == "G" or x == "C" or x == "T")])
    if count == seq_len:
        print("all nucleotides defined")
    else:
        print("not all nucleotides defined")

def compare_AA_seqs(*args):
    print("Molecular weight:")
    for aa_seq in args:
        print("| %0.2f" % aa_seq.molecular_weight(),end=" | ")

    print("\nAromaticity:")
    for aa_seq in args:
        print("| %0.2f" % aa_seq.aromaticity(),end=" | ")

    print("\nInstability Index:")
    for aa_seq in args:
        print("| %0.2f" % aa_seq.instability_index(),end=" | ")

    print("\nIsoelectric Point:")
    for aa_seq in args:
        print("| %0.2f" % aa_seq.isoelectric_point(),end=" | ")

    print("\nSecondary Structure Fraction:")
    for aa_seq in args:
        sec_struc_first = aa_seq.secondary_structure_fraction()
        print("| %0.2f" % sec_struc_first[0],end=" | ")  # helix

    print("\nMolar Extinction coefficient:")
    for aa_seq in args:
        first_epsilon_prot = aa_seq.molar_extinction_coefficient()  # [reduced, oxidized]
        print("| " + str(first_epsilon_prot[0]),end=" | ")  # with reduced cysteines

def clean_AA_seq(seq):
    seq = str(seq).replace("*", "")
    seq = seq.replace("X","")
    return seq


# MAIN FUNCTION
def main():
    start = 31094927
    end = 31377677
    length = end-start
    # print(length)
    # print(length/3)
    initial_record = None
    first_SNP_location = 31357033
    first_SNP_record = None
    second_SNP_location = 31095311
    second_SNP_record = None


    # Get the FASTA sequence of NF1
    for record in SeqIO.parse("gene1.fasta", "fasta"):
        initial_record = record
    verify_all_nucleotides(initial_record)

    first_SNP_record = copy_record(initial_record)
    second_SNP_record = copy_record(initial_record)


    if verify_nucleotide(first_SNP_record,first_SNP_location-start,"C"):
        first_SNP_record = change_nucleotide(first_SNP_record,first_SNP_location-start,"T")
        print("FIRST SNP SEQ replace C with T successfully")
    if verify_nucleotide(second_SNP_record,second_SNP_location-start,"T"):
        second_SNP_record = change_nucleotide(second_SNP_record,second_SNP_location-start,"A")
        print("SECOND SNP SEQ replace T with A successfully")

    initial_complete = complete_nucleotide_codons_N(initial_record)
    first_SNP_record = complete_nucleotide_codons_N(first_SNP_record)
    second_SNP_record = complete_nucleotide_codons_N(second_SNP_record)

    initial_complete_seq_2_Protein = initial_complete.seq.translate()
    first_SNP_seq_2_Protein = first_SNP_record.seq.translate()
    second_SNP_seq_2_Protein = second_SNP_record.seq.translate()
    # print(first_SNP_seq_2_Protein)
    # print(second_SNP_seq_2_Protein)

    # print(len(initial_complete_seq_2_Protein))
    # print(len(first_SNP_seq_2_Protein))
    # print(len(second_SNP_seq_2_Protein))
    # print(type(first_SNP_seq_2_Protein))
    initial_clean = clean_AA_seq(initial_complete_seq_2_Protein)
    first_clean = clean_AA_seq(first_SNP_seq_2_Protein)
    second_clean = clean_AA_seq(second_SNP_seq_2_Protein)
    # print(first_clean)

    initial = ProteinAnalysis(initial_clean)
    one = ProteinAnalysis(first_clean)
    two = ProteinAnalysis(second_clean)

    compare_AA_seqs(initial,one,two)
    # NC_000017.11:31094927-31377677
    # first SNP location C>T 31357033
    # second SNP location 31095311 T>A
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
    return










#####################################################################################
if __name__ == "__main__":
    main()
else:
    print("cs223_prog_assign_2.py : Is intended to be executed and not imported.")


