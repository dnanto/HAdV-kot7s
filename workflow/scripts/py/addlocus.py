#!/usr/bin/env python3

import sys
import re
from Bio import SeqIO

record = SeqIO.read(sys.argv[1], "genbank")
for feature in record.features:
  if feature.type == "CDS":
    feature.qualifiers["locus_tag"] = [
      re.sub(
        r"[^a-zA-Z0-9*_-]+",
        "_",
        feature.qualifiers.get("gene", feature.qualifiers["product"])[0]
      )
    ]

SeqIO.write(record, sys.stdout, "genbank")
