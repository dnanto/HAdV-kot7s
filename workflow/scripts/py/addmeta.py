#!/usr/bin/env python3

import logging
import sys
from collections import OrderedDict
from csv import DictReader, DictWriter
from traceback import format_exception

from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid


def excepthook(exc_type, exc_value, exc_traceback):
    logger.error("".join(format_exception(exc_type, exc_value, exc_traceback)))


def main(snakemake):
    delimiter = "\t"
    with open(snakemake.input.meta) as file:
        reader = DictReader(file, delimiter=delimiter)
        meta = OrderedDict((row["strain"], row) for row in reader)
    with open(snakemake.output.meta, "w") as file:
        fields = ("seguid", "unknown", "clade_membership")
        writer = DictWriter(file, (*reader.fieldnames, *fields), delimiter=delimiter)
        writer.writeheader()
        for record in SeqIO.parse(snakemake.input.seqs, "fasta"):
            pct_unk = (record.seq.upper().count("N")) / len(record)
            clade = f"gt{meta[record.id]['genotype']}"
            writer.writerow(
                {
                    **meta[record.id],
                    **dict(zip(fields, (seguid(record.seq), pct_unk, clade))),
                }
            )
    return 0


if __name__ == "__main__":
    logging.basicConfig(
        filename=snakemake.log.log, level=logging.INFO, format="%(asctime)s %(message)s"
    )
    sys.excepthook = excepthook
    sys.exit(main(snakemake))
