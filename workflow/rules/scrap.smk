rule mutation_r:
    message:
        "Reconstruct ancestral sequences and mutations: {wildcards.species}."
    input:
        msa=rules.align.output,
        tree=rules.refine.output.tree,
    output:
        seqs=root_out / "{species}" / "mut.fasta",
        json=root_out / "{species}" / "mut.nt.json",
        tsv=root_out / "{species}" / "mut.seq.model.tsv",
    log:
        out=root_out / "{species}" / "mut.nt.out.log",
        msg=root_out / "{species}" / "mut.nt.msg.log",
    params:
        config["seed"],
    threads: 16
    conda:
        "envs/hadv-r.yaml"
    script:
        "scripts/R/ancestral.R"


rule mutation_augur:
    message:
        "Reconstruct ancestral sequences and mutations: {wildcards.species}."
    input:
        msa=rules.align.output,
        tree=rules.refine.output.tree,
    output:
        json=root_out / "{species}" / "mut.nt.json",
    log:
        out=root_out / "{species}" / "mut.nt.out.log",
        msg=root_out / "{species}" / "mut.nt.msg.log",
    params:
        config["seed"],
    threads: 16
    conda:
        "envs/hadv-augur.yaml"
    shell:
        "augur ancestral --alignment {input.msa:q} --tree {input.tree:q} --output-node-data {output.json:q}"


rule translate:
    message:
        "Translating amino acid sequences: {wildcards.species}."
    input:
        tree=rules.refine.output.tree,
        json=rules.mutation_augur.output.json,
        gff3=rules.files.params.gff3,
    output:
        json=root_out / "{species}" / "mut.aa.json",
    log:
        out=root_out / "{species}" / "mut.aa.out.log",
        err=root_out / "{species}" / "mut.aa.err.log",
    threads: 1
    conda:
        "envs/hadv-augur.yaml"
    shell:
        """
        augur translate \
          --tree {input.tree:q} \
          --ancestral-sequences {input.json:q} \
          --reference-sequence {input.gff3:q} \
          --output-node-data {output:q} \
          1> {log.out:q} \
          2> {log.err:q}
        """
