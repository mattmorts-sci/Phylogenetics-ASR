"""
Created by Matt Mortimer based off the sequence
cleaner script from BioPython.
26 Oct 21
matthew.mortimer@anu.edu.au
ORCID id: https://orcid.org/0000-0002-8135-9319
Python 3

Version 1.3.3 (211118)
"""

from Bio import SeqIO
from modules.log import log
from datetime import datetime

# Formats the current date, assigns to variable date
DATE = datetime.now().strftime("%y%m%d")


def cleaner(file, data_source, project, por_n=0):

    """
    Modification of the sequence cleaner script from BioPython.
    Takes four arguments (one default):
        The fasta file name as a string with the file extension,
        The name of the data source as a sting
        The project as a string
        por_n value as the acceptable freq of 'X' residues (default 0)
    Outputs a file with the duplicate sequeces and sequences with 'X' residues
        (exceeding the set freq)
    Summary is outputed to log.
    """

    # Creates 2 dicts and set 2 counters to 0
    sequences = {}  # Dict as a hash table
    duplicate = {}
    dup_count = 0
    d_count = 0

    # Uses the Biopython fasta parse
    for seq_record in SeqIO.parse(file, "fasta"):
        # Take the current sequence
        sequence = str(seq_record.seq).upper()

        # If the sequence passed in the test "is it clean?" and it isn't in
        # the hash table, the sequence and its id are going to be in the hash
        if (float(sequence.count("X")) / float(len(sequence))) * 100 <= por_n:
            if sequence not in sequences:
                sequences[sequence] = seq_record.description
            else:  # Places duplicate sequences into the 'duplicate' dict
                duplicate[sequence] = seq_record.id
                dup_count += 1
        else:  # Places sequences with ambiguous residues into the same dict
            duplicate[sequence] = seq_record.id
            dup_count += 1
            # A count is required to check if some sequences were duplicated
            # multiple times

    # Write the clean sequences
    # Create a file with the cleaned sequences
    with open(
        f"{project}/output/{DATE}_{data_source}_cleaned.fasta", "w+"
    ) as output_file:
        # Just read the hash table and write on the file as a fasta format
        # as key and variable
        for sequence in sequences:
            output_file.write(">" + sequences[sequence] + "\n" + sequence + "\n")

    # Creates a file for dropped sequences using the same methods as above but
    # taking the sequences from the 'duplicate' dict.
    with open(
        f"{project}/output/{DATE}_{data_source}_cleaner_dropped_seq.fasta", "w+"
    ) as dropped_file:
        for sequence in duplicate:
            dropped_file.write(">" + duplicate[sequence] + "\n" + sequence + "\n")
            d_count += 1

    # Printing and logging the summary
    log(
        f"Sequence cleaner used to remove duplicate sequences and those with \
invalid amino acids, output was written to {project}/{DATE}_\
{data_source}_cleaned.fasta. {d_count} sequences were dropped \
and written to {project}/output/{DATE}_{data_source}_cleaner_\
dropped_seq.fasta, there were {dup_count} instances of duplicate \
or ambiguous sequences"
    )
    print(
        f"{len(sequences)} were kept, output was written to {project}/{DATE}\
_{data_source}_cleaned.fasta. {d_count} sequences were dropped and \
written to {project}/output/{DATE}_{data_source}_cleaner_dropped_seq.fasta"
    )
    print(f"There were {dup_count} instances of duplicate or ambiguous sequences.")

    seq_file = f"{project}/output/{DATE}_{data_source}_cleaned.fasta"

    print(seq_file)


def reviewed_unreviewed(file, file2, data_source, project):
    """
    Function takes two fasta files, combining both but dropping duplicates but
    only from the second file.
    Takes four arguements:
    file (all sequences are kept)
    file2 (sequences duplicated with the first file are dropped, rest are kept)
    data_source
    project
    """

    sequences = {}
    duplicate = {}
    d_count = 0

    for seq_record in SeqIO.parse(file, "fasta"):
        # Load the dict 'sequences' with all records/seq from the first file
        sequence = str(seq_record.seq).upper()

        sequences[sequence] = seq_record.description

    # Take the sequences from the second file and add to the dict if not already
    # in it. Duplicates are added to the dict 'duplicate'
    for seq_record in SeqIO.parse(file2, "fasta"):
        # Take the current sequence
        sequence = str(seq_record.seq).upper()

        # If the sequence passed in the test "is it clean?" and it isn't in
        # the hash table, the sequence and its id are going to be in the hash
        if sequence not in sequences:
            sequences[sequence] = seq_record.description
        else:  # Places duplicate sequences into the 'duplicate' dict
            duplicate[sequence] = seq_record.id
            d_count += 1

    # Write the clean sequences
    # Create a file with the cleaned sequences
    with open(
        f"{project}/output/{DATE}_{data_source}_rev_unrev_deduped.fasta",
        "w+",
    ) as output_file:
        # Just read the hash table and write on the file as a fasta format
        # as key and variable
        for sequence in sequences:
            output_file.write(">" + sequences[sequence] + "\n" + sequence + "\n")

    # Creates a file for dropped sequences using the same methods as above but
    # taking the sequences from the 'duplicate' dict.
    with open(
        f"{project}/output/{DATE}_{data_source}_rev_unrev_dropped.fasta", "w+"
    ) as dropped_file:
        for sequence in duplicate:
            dropped_file.write(">" + duplicate[sequence] + "\n" + sequence + "\n")
            d_count += 1

    # Printing and logging the summary
    log(
        f"Sequence cleaner used to remove sequences from {file2} that are\
duplicated in {file} and the datasets combined, output was written to \
{project}/output/{DATE}_{data_source}_rev_unrev_deduped.fasta. {d_count} \
sequences were dropped and written to {project}/output/{DATE}_{data_source}\
_rev_unrev_dropped.fasta"
    )
    print(
        f"{len(sequences)} sequences were kept, output was written to {project}/\
{DATE}_{data_source}_rev_unrev_deduped.fasta. {d_count} sequences were \
dropped and written to output/test_joint_cleaned_cleaner_dropped_seq.fasta"
    )

    seq_file = f"{project}/output/{DATE}_{data_source}_rev_unrev_deduped.fasta"
    print("")
    print(seq_file)
