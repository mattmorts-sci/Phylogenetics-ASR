"""
Created by Matt Mortimer
26 Oct 21
matthew.mortimer@anu.edu.au
ORCID id: https://orcid.org/0000-0002-8135-9319
Python 3

Version 2.1.0 (211117)
"""

from Bio import SeqIO
import pandas as pd
import seaborn as sns


def len_distro(file, bin=5, x_size=15, y_size=6):
    """
    Plots the distribution of sequence length from the input fasta file
    Takes the file name (str), bin size (int) and x & y size (int) of the
    chart as arguments.
    """

    # Create a pandas dataframe with the sequence id and sequence length,
    # then visulise in a chart

    # Set dict for handling sequences
    len_lst = []

    # Parses fasta file using BioPython, assigns sequences to dict
    # with header as key and sequence as value
    for seq_record in SeqIO.parse(f"{file}", "fasta"):
        sequence = str(seq_record.seq).upper()
        len_lst.append(len(sequence))

    # Creates a pandas dataframe from the dict
    len_df = pd.DataFrame(len_lst, columns=["Sequence_length"])

    # Visulise distribution of sequence lengths using seaborn
    sns.set(rc={"figure.figsize": (x_size, y_size)})

    sns.displot(data=len_df, x="Sequence_length", kde=True, binwidth=bin)
