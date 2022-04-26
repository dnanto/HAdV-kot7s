rule target_auspice:
    input:
        expand(root_aus / "{species}.json", species=speciess),


rule target_bactdate:
    input:
        expand(
            root_out / "{species}" / "bac.{ext}",
            species=speciess,
            ext=("tree", "json", "tsv"),
        ),


rule target_gubbins:
    input:
        expand(root_out / "{species}" / "gub.final_tree.tre", species=speciess),


rule target_iqtree:
    input:
        expand(root_out / "{species}" / "iqt.treefile", species=speciess),


rule target_align:
    input:
        expand(root_out / "{species}" / "seq.aligned.fasta", species=speciess),


rule target_filter:
    input:
        expand(root_out / "{species}" / "seq.fasta", species=speciess),


rule target_addmeta:
    input:
        expand(root_out / "{species}" / "meta.tsv", species=speciess),
