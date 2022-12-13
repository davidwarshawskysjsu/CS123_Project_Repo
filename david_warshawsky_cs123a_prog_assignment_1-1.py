# -*- coding: utf-8 -*-
# NAME: cs223_prog_exercise_2.py

# IMPORTS
import os
import re
import sys
from copy import copy
from Bio import Entrez
from Bio import SeqIO
import xmltodict
import itertools
from collections import OrderedDict


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
                # print("get_valid_ids: The Entrez.efetch id# " + str(lid) + " for NCBI " + str(db_name) +
                #       " DB, and return type = " + str(return_type) + " is not valid.")
                copy_list_of_ids.remove(lid)

    # When we reach here, copy_list_of_ids will only have valid ids
    return copy_list_of_ids

# Gets the accession numbers from the assignment text file and returns them as a list of strings
def get_accession_nums():
    acc_nums = []
    acc_file = open("acc_nums.txt", "r")
    acc_lines = acc_file.readlines()
    for line in acc_lines:
        if line.strip() == "":
            continue
        else:
            acc_nums.append(line.strip())
    return acc_nums



# Returns valid handle read ids which can be used to get the fafsa files.
def get_valid_handle_read_ids(acc_num):
    # # Entrez.esearch for info related to the current accession number
    handle = Entrez.esearch(db='gene', term=str(acc_num))
    #
    #
    # # Entrez.read the result returned from NCBI search for the current accession number
    handle_read = Entrez.read(handle)
    #
    # # Get Idlist from previous read of accession number info
    handle_read_ids = handle_read['IdList']
    # print("handle_read_ids: ", handle_read_ids)
    #
    # # Make sure IDs in Idlist are valid. Return a list of only valid IDs.
    valid_handle_read_ids = get_valid_ids(handle_read_ids, 'nucleotide', 'fasta', 'efetch')
    return valid_handle_read_ids

# Converts the dictionary of species and percentage differences to a sorted list from most likely to be
# the species because it is the smallest percentage difference in absolute value which means it matches
# best to highest percentage difference which means it matches least.
def convert_dict_to_list(dictionary):
    sorted_dict = sorted(dictionary.items(),key=lambda x: x[1],reverse=False)
    the_list = []
    for i in range(len(sorted_dict)):
        the_list.append([k for k,v in sorted_dict if v == sorted_dict[i][1]])
    values = OrderedDict((tuple(x), x) for x in the_list).values()
    return list(values)


def main():
    """The main program collects """
    acc_nums = get_accession_nums()
    # # Init Entrez email addr to use
    Entrez.email = "David.Warshawsky@sjsu.edu"
    #
    # # Extract Entrez info for the specified ACC number
    for acc_num in acc_nums:
        print(acc_num)
        valid_handle_read_ids = get_valid_handle_read_ids(acc_num)
        # print(valid_handle_read_ids)
        if valid_handle_read_ids == []:
            print("   epigen_pipeline_get_entrez_info: No valid IDs found for ACC #" + str(acc_num) + ".")
            return
        #
        # # If we get here, valid_handle_read_ids has usable ID nums.
        # This is for testing just one of them.
        # valid_handle_read_ids = [valid_handle_read_ids[0]]
        sum = 0
        for acc_handle_id in valid_handle_read_ids:
            print("ACC: ", str(acc_num))

            # Entrez.efetch based on ID#
            handle_read_id_info = Entrez.efetch(db="gene", id=int(acc_handle_id), rettype="fasta")
            # SeqIO.read  fasta info
            record_seq = SeqIO.read(handle_read_id_info, "fasta")
            print("record_seq.seq: ", record_seq.seq)
            print("sequence length" + str(len(record_seq)))
            print(len(record_seq.seq))
            sum += len(record_seq.seq)
        print(sum)
        #
    print()
    print ("############### DONE ###############")
    return



def retrieve_annotation(id_list):

    """Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
    submit the data to NCBI) and esummary to retrieve the information.
    Returns a list of dictionaries with the annotations."""
    Entrez.email = "David.Warshawsky@sjsu.edu"
    print(id_list)
    request = Entrez.epost("gene", id=id_list[0])
    try:
        result = Entrez.read(request)
    except RuntimeError as e:
        # FIXME: How generate NAs instead of causing an error with invalid IDs?
        print("An error occurred while retrieving the annotations.")
        print("The error returned was %s" % e)
        sys.exit(-1)

    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    data = Entrez.esummary(db="gene", webenv=webEnv, query_key=queryKey)
    annotations = Entrez.read(data)

    print("Retrieved %d annotations for %d genes" % (len(annotations), len(id_list)))

    return annotations

def print_data(annotation):
    for gene_data in annotation:
        gene_id = gene_data["Id"]
        gene_symbol = gene_data["NomenclatureSymbol"]
        gene_name = gene_data["Description"]
        print("ID: %s - Gene Symbol: %s - Gene Name: %s" % (
            gene_id,
            gene_symbol,
            gene_name,
        ))





#####################################################################################
if __name__ == "__main__":
    id_list = get_accession_nums()
    print_data(retrieve_annotation(id_list))
    # main()
else:
    print("cs223_prog_assign_2.py : Is intended to be executed and not imported.")


