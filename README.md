# HAdV-kot7s
HAdV ποΈποΈποΈποΈποΈποΈποΈ Keeper of the Seven Species

# Quickstart

Create environment:
```
mamba env create -f workflow/envs/hadv.yaml && \
conda activate hadv
```

Run the workflow for HAdV-A:
```
snakemake --config seqs=data/A.fasta.gz --printshellcmds --configfile config/test.yaml --set-threads iqtree=4 --cores 10 target_auspice
```

Run the server:
```
HOST="localhost" auspice view --datasetDir auspice
```

# Cite
```
NegrΓ³n, D. A. (2021). Molecular Clock Analysis of Human Adenovirus [Ph.D., George Mason University]. In ProQuest Dissertations and Theses. http://www.proquest.com/docview/2572612100/abstract/9972F7D348C34013PQ/1
```
