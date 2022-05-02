#!/usr/bin/env 

species="${1:-A}"
iter="${2:-3}"

for (( i=1; i<=iter; i++ ))
do    
    root="experiment/$species/$i" && \
    
    mkdir -p "$root" && \

    augur filter \
      --metadata data/meta.tsv \
      --sequences "data/$species.fasta.gz" \
      --exclude-where "date=?" \
      --group-by genotype \
      --sequences-per-group 3 \
      --subsample-seed "$RANDOM" \
      --output "$root/ref.fasta" && \

    grep \> "$root/ref.fasta" | \
      tr -d \> | \
      cut -f 1 -d " " > \
      "$root/exclude.txt" && \

    augur filter \
      --metadata data/meta.tsv \
      --sequences "data/$species.fasta.gz" \
      --exclude-where "date=?" \
      --group-by genotype \
      --sequences-per-group 3 \
      --subsample-seed 7080 \
      --output "$root/qry.fasta" && \

    snakefmt -v workflow/rules/*.smk workflow/Snakefile workflow/scripts/py/* && \
    snakemake \
      --printshellcmds \
      --config \
        seed="$RANDOM" \
        seqs="$root/ref.fasta" \
        fasta="data/$species.fasta" \
        gff3="data/$species.gff3" \
        genbank="data/$species.gb" \
        sample="$root/qry.fasta" \
        qc=config/qc.json \
        out="$root" \
      --configfile config/test.yaml \
      --set-threads iqtree=4 \
      --cores 10 \
      target_nextclade > \
      "$root/snakemake.log"
done
