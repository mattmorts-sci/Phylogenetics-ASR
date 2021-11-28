"""
Created by Matt Mortimer
26 Oct 21
matthew.mortimer@anu.edu.au
ORCID id: https://orcid.org/0000-0002-8135-9319
Python 3

Version 1.3.2 (211129)
"""

from Bio import SeqIO
from bioservices import UniProt
from datetime import datetime, timedelta
from modules.log import log
import pandas as pd
from modules.progress import update_progress
from time import time
from modules.utilities import dict_to_fasta

date = datetime.now().strftime("%y%m%d")


def annotations(file, data_source, project, delim="|", el_num=1):
    """
    Downloads sequence annotations from uniprot for all sequences in the
    fasta file.
    Takes 3 arguments:
        a fasta file
        data source name
        project name
    This splits the fasta header by '|' and takes the second element as
    the UniProt Entry. If the fasta header is different, argument defults
    delim (delimeter) and el_num (element number) will need to be changed
    *Modified from code contributed by Matt Spence*
    """

    # Start point for run time calc.
    start = time()

    # This variable forms the column names in the resulting csv, it
    # will need to be modified depending on the data set.
    annotations = "id, length, mass, genes(PREFERRED), reviewed, lineage(SUPERKINGDOM), lineage(KINGDOM), lineage(PHYLUM), lineage(CLASS), lineage(ORDER), lineage(FAMILY), lineage(GENUS), lineage(SPECIES), lineage-id(SPECIES), database(chebi), database(chebi(Catalytic activity)), database(chebi(Cofactor)), database(pdb), database(pfam)"
    columns = "id\t\
length\t\
mass\t\
genes\t\
reviewed\t\
Superkingdom\t\
Kingdom\t\
Phylum\t\
Class\t\
Order\t\
Family\t\
Genus\t\
Species\t\
Species_id\t\
database(chebi)\t\
chebi_Catalytic_activity\t\
chebi_Cofactor\t\
database_pdb\t\
database_pfam"

    # Specifies the BioServices service used, creates a counter
    service = UniProt()
    anno_count = 0

    # Creates the output file, will overwrite if this code is run more than
    # once

    tot = len(tuple(SeqIO.parse(file, "fasta")))

    with open(
        f"{project}/output/{date}_{data_source}_filt_\
annotations.txt",
        "w+",
    ) as anno_test:
        anno_test.write(columns)  # Writes the column names to the csv
        # Parses the fasta file and isolates the UniProt Entry from the header
        # as i. Passes i to the search string which is assigned to variable n
        with open(file, "r") as seqs:
            for seq_record in SeqIO.parse(seqs, "fasta"):
                i = seq_record.id.split(delim)
                try:
                    n = service.search(i[el_num], columns=annotations, frmt="tab")
                except IndexError:
                    print("IndexError: Fasta header doesn't match expected format")
                    print("Caught error at:")
                    print(i)
                # Executes the search and writes to file.
                # First splits the response into 2 elements using '\n
                # dumps the first element which would produce duplicate
                # column names in as a new row for each record. The
                # second element is the record annotation which is
                # written to file.
                else:
                    try:
                        anno_test.write(n.split("\n")[1] + "\n")
                    except IndexError:
                        print("IndexError: Annotation doesn't match expected format")
                        print("Last annotation was:")
                        print(f"Accession: {n}")
                        print(f"Sequence number : {anno_count}")
                    else:
                        anno_count += 1
                        progress = anno_count / tot
                        remaining = tot - anno_count
                        update_progress(progress, anno_count, remaining, start)

    # Run time calc
    end = time()
    elapsed_time = end - start
    f_elapsed_time = str(timedelta(seconds=elapsed_time))
    run_time = f_elapsed_time.split(".")[0]

    # Printing and logging the summary
    print(f"Run time: {run_time} hh:mm:ss")

    print(
        f"Annotations were retrieved or attempted for {anno_count} sequences \
from {file} and outputted to '{project}/output/{date}_{data_source}_filt_\
annotations.txt'"
    )
    log(
        f"Annotations were retrieved or attempted for {anno_count} sequences \
from {file}. Run time took {run_time} hh:mm:ss."
    )


def indexing(anno_file, seq_file, project, delim="|", el_num=1):
    """
    Creates a master index by adding the sequences to the annotation
    data.
    Takes 3 arguments:
        The annotation file
        A fasta file
        The project name
    Data is outputted to a master index file.
    This function is specific to the data set and may need to be
    tweaked depending on the format of the annotation file and fasta
    file.
    """

    # Dict for storing the fasta sequences
    seq_dict = {}

    # Parses the annotation file and stores as a pandas dataframe
    anno_df = pd.read_csv(anno_file, sep="\t")

    # Parses the fasta file saving the header to the dict key and
    # sequence as the value
    for seq_record in SeqIO.parse(seq_file, "fasta"):
        sequence = str(seq_record.seq).upper()
        # Fasta header is formated, spliting on '_' and keeping
        # the first element
        k = seq_record.id.split(delim)
        seq_dict[k[el_num]] = sequence

    # Creates a dataframe from the above dict
    seq_df = pd.DataFrame.from_dict(seq_dict, orient="index")
    seq_df = seq_df.reset_index()
    seq_df.columns = ["Entry", "Sequence"]

    # Merges the annotation and fasta dataframes on 'Entry' columns
    # as the key. Inner merge will drop data without a matching 'Entry'
    # value.
    master_index = pd.merge(anno_df, seq_df, how="inner", on="Entry")

    # Writes the above, merged, dataframe to file
    output_file = f"output/{date}_{project}_PFAM_master_index.csv"
    master_index.to_csv(output_file, sep=",", index=False)

    # Prints and logs summary
    print(
        f"Master index containing annotation and sequence data created \
with {len(master_index.index)} records."
    )
    print(f"Outputtted file called '{output_file}'")
    log(
        f"Master index containing annotation and sequence data created \
with {len(master_index.index)} records, from {anno_file} and \
{seq_file}."
    )


def reviewed(file, PROJECT):
    """
    Seperate dataset into two, 'reviewed' and 'unreviewed' (column 'Status')
    Takes two arguments, the file name and the project as strings
    Outputs two fasta files
    """

    # This is important in order to ensure we do not filter out the
    # reviewed sequences.

    index_df = pd.read_csv(file)

    reviewed_df = index_df[
        index_df["Status"] == "reviewed"
    ]  # Make subset df of only reviewed sequences
    unreviewed_df = index_df[index_df["Status"] != "reviewed"]

    # Assigns length of the dataframes (proxy for number of sequencess to
    # variables)
    rev = len(reviewed_df)
    unrev = len(unreviewed_df)

    print(f"There are {rev} reviewed sequences in the dataset")
    print(f"There are {unrev} unreviewed sequences in the dataset")

    # Prepares two dicts to make fasta files for reviewed and unreviewed sequences
    # with the Entry as the key and the sequence as the variable
    reviewed_dict = dict(zip(reviewed_df["Entry"], reviewed_df["Sequence"]))
    unreviewed_dict = dict(zip(unreviewed_df["Entry"], unreviewed_df["Sequence"]))

    # Write the dicts to fasta file

    # Reviewed sequences
    # Variable with the name of the outputted fasta file
    output_file = f"{PROJECT}/output/{date}_reviewed.fasta"
    # Runs the dict_to_fasta function in the utilities module
    dict_to_fasta(reviewed_dict, output_file)

    # Unreviewed sequences
    # Variable with the name of the outputted fasta file
    output_file = f"{PROJECT}/output/{date}_unreviewed.fasta"
    # Runs the dict_to_fasta function in the utilities module
    dict_to_fasta(unreviewed_dict, output_file)

    log(f"{rev} sequences were written to {PROJECT}/output/{date}_reviewed.fasta")
    log(f"{unrev} sequences were written to {PROJECT}/output/{date}_unreviewed.fasta")

    print(f"{rev} sequences were written to {PROJECT}/output/{date}_reviewed.fasta")
    print(f"{unrev} sequences were written to {PROJECT}/output/{date}_unreviewed.fasta")
