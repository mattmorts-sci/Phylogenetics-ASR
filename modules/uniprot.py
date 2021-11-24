"""
Created by Matt Mortimer
26 Oct 21
matthew.mortimer@anu.edu.au
ORCID id: https://orcid.org/0000-0002-8135-9319
Python 3

Version 1.0.2 (211129)
"""


def uniprot(pfam_id, output, sp_output, source="UniProt"):
    """
    Put discription in here
    """

    from bioservices import UniProt, PDB
    from Bio import SeqIO
    from datetime import timedelta
    from log import log
    from time import time

    # Start point for run time calc.
    start = time()

    # Set which service to use (default is UniProt)
    if source.lower() == "uniprot":
        service = UniProt()
    elif source.lower() == "pdb":
        service = PDB()
    else:
        print("Source/service must be entered as str, either UniProt or PDB")

    # Establishes the BioServices command to execute
    command = f"database:(type:pfam {pfam_id})"
    ds = service.search(command, frmt="fasta")
    # Executes the command writing to file
    with open(output, "w+") as pfam:
        pfam.write(ds)

    # Find the number of sequences in the outputted file
    seq_num = len(tuple(SeqIO.parse(output, "fasta")))

    # Run time calc
    end = time()
    elapsed_time = end - start
    f_elapsed_time = str(timedelta(seconds=elapsed_time))
    run_time = f_elapsed_time.split(".")[0]

    # Printing and logging the summary
    print(f"Download time: {run_time} hh:mm:ss")

    print(f"{seq_num} sequences downloaded from PFAM family {pfam_id}")
    print(f"Sequences outputted to {output}")

    log(
        f"{seq_num} sequences downloaded (in {run_time} hh:mm:ss) from PFAM \
family {pfam_id}; outputted to {output}"
    )

    # Write all UniProt sequence ids to file
    sp_count = 0

    with open(sp_output, "w+") as ids:
        for seq_record in SeqIO.parse(output, "fasta"):
            sequence = str(seq_record.seq).upper()
            i = seq_record.id.split("|")
            if i[0] == "sp":
                ids.write(">" + seq_record.description + "\n")
                ids.write(sequence + "\n")
                sp_count += 1

    print(f"{sp_count} UniProt ids written to {sp_output}")
    log(f"{sp_count} UniProt ids written to {sp_output}")
