from pathlib import Path

from BCBio import GFF
from snakemake.utils import validate


validate(config, "../schemas/config.yaml")

root_cfg = Path("config")
root_dat = Path("data")
root_out = Path("results")
root_aus = Path("auspice")

genes = {}
speciess = []
for path in Path().glob(config["glob"]):
    # parse species
    species = path.name.split(".", maxsplit=1)[0]
    speciess.append(species)
    # parse genes
    genes[species] = ",".join(
        feat.qualifiers["gene"][0]
        for rec in GFF.parse(root_dat / f"{species}.gff3")
        for feat in rec.features
        if feat.type == "gene"
    )
