"""
Created by Matt Mortimer
26 Oct 21
matthew.mortimer@anu.edu.au
ORCID id: https://orcid.org/0000-0002-8135-9319
Python 3

Version 1.0.0 (21126)
"""


def dict_to_fasta(fasta_dict, outfile):
    """
    Takes a dictionary where the key is the string for the sequence header
    and the value is the sequence and creates a fasta file.
    Takes two arguments the dict and then the name of the output file.
    """
    from modules.log import log

    count = 0
    # Creates a fasta file from the dict above (fasta_dict)
    with open(outfile, "w") as form:
        for k, v in fasta_dict.items():
            form.write(">" + k + "\n" + v + "\n")
            count += 1

    log(f"{count} sequences were written to {outfile}")


def fasta_headers(file):
    """
    Returns a list of headers from the fasta file, prints the first 5
    Useful for checking header format.
    Takes one argument, the file name as a string.
    """

    with open(file, "r") as fasta:
        header_lst = [
            line.lstrip(">").rstrip() for line in fasta if line.startswith(">")
        ]

    print("First five fasta headers:")
    n = 0
    for i in header_lst:
        if n < 6:
            print(i)
            n += 1
