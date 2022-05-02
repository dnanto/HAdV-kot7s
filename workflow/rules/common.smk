from pathlib import Path

from BCBio import GFF
from snakemake.utils import validate


validate(config, "../schemas/config.yaml")

root_cfg = Path("config")
root_dat = Path("data")
root_out = Path(config["out"])

genes = {}
groups = []
for path in Path().glob(config["seqs"]):
    # parse group
    group = path.name.split(".", maxsplit=1)[0]
    groups.append(group)
    # parse genes
    genes[group] = ",".join(
        feat.qualifiers["gene"][0]
        for rec in GFF.parse(config.get("gff3", root_dat / f"{group}.gff3"))
        for feat in rec.features
        if feat.type == "gene"
    )
