"""
Created by Matt Mortimer
26 Oct 21
matthew.mortimer@anu.edu.au
ORCID id: https://orcid.org/0000-0002-8135-9319
Python 3

VersionL 1.2.0 (211128)
"""

from datetime import datetime


def log(log_data):
    """
    A simple logging function. Whatever string argument is passed to the
    function it will append to the 'log.txt' file with a date-time stamp
    Use: insert a string as the argument to the function, the string will
    be added to the log
    """
    #    if __name__ == "__main__":

    # Sets the format for the current date-time and assigns to variable
    # date_time
    date_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Writes the string to the log with the date-time
    with open("log.txt", "a") as logging:
        logging.write(date_time + "\t" + log_data + "\n")
