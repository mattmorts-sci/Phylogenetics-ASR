"""
Modified by Matt Mortimer from rayryeng's code @
https://stackoverflow.com/questions/3160699/python-progress-bar
29 Oct 21
matthew.mortimer@anu.edu.au
ORCID id: https://orcid.org/0000-0002-8135-9319
Python 3

Version 1.0.1 (211129)
"""

import sys

# from datetime import timedelta
from time import time


def update_progress(progress, complete, remaining, start):
    """
    update_progress() : Displays or updates a console progress bar
    Accepts a float between 0 and 1. Any int will be converted to a float.
    A value under 0 represents a 'halt'.
    A value at 1 or bigger represents 100%
    """

    barLength = 20  # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done!\r\n"

    # Varibles for the progress bar
    block = int(round(barLength * progress))
    per_progress = progress * 100
    time_now = time()
    elapsed_time = time_now - start
    # while True:
    #     if remaining > 0:
    time_est = elapsed_time / complete * remaining / 60

    # Formats the progress bar
    text = f'\rPercent compelete: [{"#" * block + "-" * (barLength - block)}]\
 {per_progress:.2f}% {status} Estimated time left {time_est:.2f} min '

    # Writes to system, flush forces it to write immediately(?)
    sys.stdout.write(text)
    sys.stdout.flush()
