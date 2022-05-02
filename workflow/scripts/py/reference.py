#!/usr/bin/env python3

import logging
import re
import sys
from collections import defaultdict
from copy import deepcopy
from traceback import format_exception

from BCBio import GFF
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation


def excepthook(exc_type, exc_value, exc_traceback):
    logger.error("".join(format_exception(exc_type, exc_value, exc_traceback)))


def clean_name(name):
    return re.sub(r"[^a-zA-Z0-9*_-]+", "_", name)


def main(snakemake):
    record = SeqIO.read(snakemake.input.gen, "genbank")

    # output FASTA
    SeqIO.write(record, snakemake.output.fasta, "fasta")

    gene_part = defaultdict(int)
    gene_list = []
    for feature in record.features:
        if feature.type == "CDS":
            name = clean_name(
                feature.qualifiers.get("gene", feature.qualifiers["product"])[0]
            )
            gene_part[name] += 1
            name = f"{name}-{gene_part[name]}"
            gene_list.append(name)
            feature.qualifiers["locus_tag"] = [name]

    # output genes
    with open(snakemake.output.text, "w") as file:
        print(*gene_list, sep=",", end="", file=file)

    # output GenBank
    SeqIO.write(record, snakemake.output.genbank, "genbank")

    features = []
    for feature in record.features:
        if feature.type == "CDS":
            feature = deepcopy(feature)
            for part in feature.location.parts:
                feature_part = deepcopy(feature)
                feature_part.type = "gene"
                feature_part.qualifiers["gene_name"] = feature.qualifiers["locus_tag"]
                feature_part.location = FeatureLocation(
                    part.start,
                    part.end - len(part) % 3,
                    strand=part.strand,
                    ref=part.ref,
                    ref_db=part.ref_db,
                )
                features.append(feature_part)

    record.features = features

    # output GFF3
    with open(snakemake.output.gff3, "w") as file:
        GFF.write([record], file)

    return 0


if __name__ == "__main__":
    logging.basicConfig(
        filename=snakemake.log.log, level=logging.INFO, format="%(asctime)s %(message)s"
    )
    sys.excepthook = excepthook
    sys.exit(main(snakemake))
