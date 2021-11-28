"""
Created by Matt Mortimer
26 Oct 21
matthew.mortimer@anu.edu.au
ORCID id: https://orcid.org/0000-0002-8135-9319
Python 3

Version 1.0.3 (211129)
"""

from Bio import SeqIO
from datetime import datetime
from modules.log import log


def size_filter(file, project, data_source, logic, seq_len):
    """
    Executes sequence_length_filter taking 5 arguments:
    File name
    Project name
    Data source
    Filtering logic ('greater'/'less')
    Outputs fasta file.
    This function calls on, and handles the sequence_length_filter function
    """

    # Formats the current date and assigns it to varible date
    date = datetime.now().strftime("%y%m%d")

    # Calls the sequence_length_filter function and uses the returned
    # dict to create seq_dict
    seq_dict = sequence_length_filter(file, logic, seq_len)

    # Writes the filtered sequences in seq_dict to file
    filt_file = f"{project}/output/{date}_{data_source}_size_filt.fasta"

    with open(filt_file, "w") as filt:
        for k, v in seq_dict.items():
            filt.write(">" + k + "\n" + v + "\n")

    # Printing and logging the summary
    log(
        f'{len(seq_dict)} filtered sequences were outputed to "{project}/output/{date}\
_{data_source}_size_filt.fasta"'
    )
    print(
        f'{len(seq_dict)} filtered sequences were outputed to "{project}/output/{date}\
_{data_source}_size_filt.fasta"'
    )


def sequence_length_filter(file, logic, seq_len):
    """
    Called by the size_filter function
    Filter sequences by length, returns sequences in a dict call sfilt_dict
    Takes file name, 'greater' or 'less', and sequence length (int)
    as args
    """

    sfilt_dict = {}

    # Filters sequences in fasta file by length SEQ_LEN
    if logic == "greater":
        # Biopython parses fasta file
        for seq_record in SeqIO.parse(file, "fasta"):
            # Checks sequence length
            if len(seq_record.seq) > seq_len:
                # Converts sequence to uppercase string
                sequence = str(seq_record.seq).upper()
                # Adds to dict with header as key and sequence as value
                sfilt_dict[seq_record.id] = sequence
    elif logic == "less":
        # Biopython parses fasta file
        for seq_record in SeqIO.parse(file, "fasta"):
            # Checks sequence length
            if len(seq_record.seq) < seq_len:
                # Converts sequence to uppercase string
                sequence = str(seq_record.seq).upper()
                # Adds to dict with header as key and sequence as value
                sfilt_dict[seq_record.id] = sequence
    else:
        print("Invalid input")

    # Prints and logs the summary
    log(
        f"Seqeunces from {file} were filtered to those {logic} than {seq_len} \
aa"
    )
    print(
        f"Seqeunces from {file} were filtered to those {logic} than \
{seq_len} aa"
    )

    return sfilt_dict
