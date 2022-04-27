#!/usr/bin/env python3

import sys
from Bio import SeqIO

record = SeqIO.read(sys.argv[1], "genbank")
for feature in record.features:
  if feature.type == "CDS":
      feature.qualifiers["locus_tag"] = feature.qualifiers.get("gene", feature.qualifiers["product"])

SeqIO.write(record, sys.stdout, "genbank")
