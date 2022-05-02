rule target_nextclade:
    input:
        expand(root_out / "{group}" / "next.tsv", group=groups),


rule target_properties:
    input:
        expand(root_out / "{group}" / "prop.json", group=groups),


rule target_auspice:
    input:
        expand(root_out / "auspice" / "{group}.json", group=groups),


rule target_bactdate:
    input:
        expand(
            root_out / "{group}" / "bac.{ext}",
            group=groups,
            ext=("tree", "json", "tsv"),
        ),


rule target_gubbins:
    input:
        expand(root_out / "{group}" / "gub.final_tree.tre", group=groups),


rule target_iqtree:
    input:
        expand(root_out / "{group}" / "iqt.treefile", group=groups),


rule target_align:
    input:
        expand(root_out / "{group}" / "seq.aligned.fasta", group=groups),


rule target_filter:
    input:
        expand(root_out / "{group}" / "seq.fasta", group=groups),


rule target_addmeta:
    input:
        expand(root_out / "{group}" / "meta.tsv", group=groups),
