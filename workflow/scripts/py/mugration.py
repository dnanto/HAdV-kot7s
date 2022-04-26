#!/usr/bin/env python3

import json
import sys
from argparse import ArgumentParser
from collections import defaultdict
from csv import DictReader
from pathlib import Path
from traceback import format_exception

from augur.traits import register_arguments, run
from Bio import Phylo


def excepthook(exc_type, exc_value, exc_traceback):
    print(*format_exception(exc_type, exc_value, exc_traceback), sep="", file=log_err)


def get_argv():
    snakemake = getattr(sys.modules[__name__], "snakemake", None)
    if snakemake:
        argv = (
            "--tree",
            snakemake.input.tree,
            "--metadata",
            snakemake.input.meta,
            "--output-node-data",
            snakemake.output[0],
            "--columns",
            *snakemake.params.columns.split(" "),
            "--confidence",
        )
    else:
        argv = sys.argv[1:]
    return argv


def main(argv):
    parser = ArgumentParser()
    register_arguments(parser)
    parser.epilog += "\nValid columns automatically selected..."
    args = parser.parse_args(argv)

    labs = {ele.name for ele in Phylo.read(args.tree, "newick").get_terminals()}
    cols = args.columns
    feats = defaultdict(set)

    with open(args.metadata) as file:
        for row in DictReader(file, delimiter="\t"):
            if row["strain"] in labs:
                for col in cols:
                    feats[col].add(row[col])

    args.columns = []
    for key, val in feats.items():
        if len(val) > 1:
            args.columns.append(key)

    print("columns:", args.columns)
    if args.columns:
        run(args)
    else:
        with open(args.output_node_data, "w") as file:
            json.dump(dict(node={}), fp=file)

    return 0


if __name__ == "__main__":
    snakemake = getattr(sys.modules[__name__], "snakemake", None)
    if snakemake:
        log_out = open(snakemake.log.out, "w")
        log_err = open(snakemake.log.err, "w")
    else:
        path = Path(__file__)
        log_out = open(f"{path.stem}.out.log", "w")
        log_err = open(f"{path.stem}.err.log", "w")
    sys.stdout = log_out
    sys.stderr = log_err
    sys.excepthook = excepthook
    sys.exit(main(get_argv()))
