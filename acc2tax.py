#!/usr/bin/env python

import sys

from progress.bar import ShadyBar, ChargingBar, PixelBar
from Bio import Entrez
from collections import Counter


def getAccessionNumbers(blastout):
    """Get the top 3 results from all reads in the BLAST output file.

    Arguments:
        blastout: BLAST output file.

    Return:
        Dictionary with read names and the top 3 accession numbers.
    """
    with open(blastout, "r") as bfile:
        total_lines = len(bfile.read().split('\n'))
    with open(blastout, "r") as blastfile:
        current_read_id = ""
        accession_numbers = dict()
        count = 0
        print()
        bar = PixelBar(
            'Getting accession numbers:',
            max=total_lines-1,
            suffix='%(percent)d%%'
        )
        for line in blastfile:
            line = line.split('\t')
            read_id = line[0]
            accession = line[1]
            if read_id != current_read_id:
                count = 0
                current_read_id = read_id
            if count <= 2 and current_read_id != "":
                if read_id in accession_numbers.keys():
                    accession_numbers[read_id].append(accession)
                    count += 1
                else:
                    accession_numbers[read_id] = [accession]
                    count += 1
            bar.next()
        bar.finish()
        print(str(len(accession_numbers.keys())) + " reads found.")
    return accession_numbers


def getName(accession_numbers):
    """Gets the organism name based on NCBI accession numbers.
    Reads will only be classified when the first two BLAST hits 
    are identical. When the the top 2 hits are different the 
    reads will be registered as unclassified.

    Arguments:
        accession_numbers: List of NCBI accession numbers
        from the BLAST output.

    Raises:
        IndexError: Only one result found so no top 3 can be selected.

    Returns:
        Dictionary with all found organism names and count.
    """
    acc_num = []
    identified = []
    count = 0
    for dummyread, numbers in accession_numbers.items():
        count += 1
        try:
            if numbers[0] == numbers[1]:
                acc_num.append(numbers[0])
            else:
                identified.append("unclassified")
        except IndexError:
            identified.append("unclassified")
    Entrez.email = 'someuser@mail.com'
    print()
    bar = PixelBar(
        'Getting names:',
        max=len(acc_num),
        suffix='%(percent)d%%'
    )
    sys.stdout.flush()
    for accession in acc_num:
        handle = Entrez.efetch(
            db="nucleotide",
            id=accession,
            rettype="gb",
            retmode="text"
        )
        result = handle.read().split('\n')
        for line in result:
            if 'ORGANISM' in line:
                identified.append(' '.join(line.split()[1:3]))
        bar.next()
    bar.finish()
    print()
    name_count = Counter(identified), len(identified)
    return name_count


if __name__ == "__main__":
    blastout = sys.argv[1]
    accession_numbers = getAccessionNumbers(blastout)
    name_count, total_reads = getName(accession_numbers)
    for organism, count in name_count.items():
        avg = float(count / total_reads * 100)
        print(organism + ": " + str(avg) + "% (" + str(count) + ")")
