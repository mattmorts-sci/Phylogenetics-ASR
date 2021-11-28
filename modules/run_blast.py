"""
Created by Matt Mortimer using code from Matt Spence
26 Oct 21
matthew.mortimer@anu.edu.au
ORCID id: https://orcid.org/0000-0002-8135-9319
Python 3

Version 1.1.1 (211129)
"""

import pandas as pd
import os as os
from modules.log import log
from datetime import datetime, timedelta
from time import time


def blast(infile, project, E_value_threshold="10e-10", cpus="2"):
    """
    IMPORTANT: This function assumes blast has been added to the path environment
    Runs BLAST locally using a fasta file to generate the database. A all v.
    all search is then preformed.
    Takes the following arguments:
        blast_path (full dir path for blast)
        infile as a fasta file as a string
        project as a string
        E_value_threshold (default = "10e-10")
        cpus; the number of CPU to use (default 2)
    Outputs a network file which can be viewed in Cytoscape.
    *Modified from code contributed by Matt Spence*
    """
    # Start point for run time calc.
    start = time()

    # Format the date
    date = datetime.now().strftime("%y%m%d")

    # Sets the name for the database files
    db_file = f"{project}/output/BLAST/{date}_dataset_db"

    # Make a BLAST database using the non-redundant/annotated sequence dataset
    os.system(f"makeblastdb -in {infile} -dbtype prot -out {db_file}")

    # Perform an all.vs.all BLAST within the database
    outblast = f"{db_file}_blast_{E_value_threshold}"

    # Runs blastp, calling blast+
    os.system(
        f"blastp -db {db_file} -query {infile} \
-outfmt 6 -out {outblast} -evalue {E_value_threshold} -num_threads {cpus}"
    )

    cols = [
        "Query",
        "Target",
        "% Identity",
        "Length",
        "Mistmatch",
        "Gapopen",
        "Qstart",
        "Qend",
        "Sstart",
        "Send",
        "E-value",
        "Bit-score",
    ]

    # Labels columns in a csv and writes to file
    blast_out = pd.read_csv(f"{outblast}", delimiter="\t", names=cols)

    blast_out.to_csv(
        f"{project}/output/BLAST/{date}_dataset_network_{E_value_threshold}\
.csv"
    )

    # Run time calc
    end = time()
    elapsed_time = end - start
    f_elapsed_time = str(timedelta(seconds=elapsed_time))
    run_time = f_elapsed_time.split(".")[0]

    # Printing and logging the summary
    print(f"Run time: {run_time} hh:mm:ss")

    print(
        f"Network file {project}/output/BLAST/{date}_dataset_network_\
{E_value_threshold}.csv generated"
    )
    log(
        f"All v. all BLAST search run against {infile} with an E-value \
threshold of {E_value_threshold}. Network file output/BLAST/{date}_\
{project}_dataset_network_{E_value_threshold}.csv generated with a run time\
of {run_time} hh:mm:ss"
    )
