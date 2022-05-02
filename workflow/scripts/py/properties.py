#!/usr/bin/env python3

import json
import logging
import sys
from collections import Counter, defaultdict
from csv import DictReader
from traceback import format_exception

from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid


def excepthook(exc_type, exc_value, exc_traceback):
    logger.error("".join(format_exception(exc_type, exc_value, exc_traceback)))


def main(snakemake):
    reference = SeqIO.read(snakemake.input.ref, "fasta")

    with open(snakemake.input.meta) as file:
        reader = DictReader(file, delimiter="\t")
        meta = {row["strain"]: row["clade_membership"] for row in reader}

    records = SeqIO.parse(snakemake.input.msa, "fasta")
    nclade = Counter(meta[ele.id] for ele in records)

    muts = defaultdict(list)

    records = SeqIO.parse(snakemake.input.msa, "fasta")
    for record in records:
        idx = 1
        for ref, alt in zip(str(reference.seq), str(record.seq)):
            if ref != "-" and alt != "-" and ref != alt:
                muts[(idx, alt)].append(meta[record.id])
            # inc ref position
            idx += ref != "-"

    nucMutLabelMap = defaultdict(list)
    nucMutLabelMapReverse = defaultdict(list)
    for key, val in sorted(muts.items()):
        mut = f"{key[0]}{key[1]}"
        if len(mut) > 3:
            for key, val in Counter(val).items():
                if val / nclade[key] > 0.3:
                    nucMutLabelMap[mut].append(key)
                    nucMutLabelMapReverse[key].append(mut)

    with open(snakemake.output.json, "w") as file:
        json.dump(
            dict(
                schemaVersion="1.10.0",
                nucMutLabelMap=nucMutLabelMap,
                nucMutLabelMapReverse=nucMutLabelMapReverse,
            ),
            file,
            sort_keys=True,
            indent=4,
        )

    return 0


if __name__ == "__main__":
    logging.basicConfig(
        filename=snakemake.log.log, level=logging.INFO, format="%(asctime)s %(message)s"
    )
    sys.excepthook = excepthook
    sys.exit(main(snakemake))
