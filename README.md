# In utero caffeine exposure and the adult mouse heart: RNA-seq reanalysis

Code for the manuscript *"In Utero Caffeine Exposure Alters Cardiac Gene Expression and Metabolic and Extracellular Matrix Pathways in Adult Mice: A Transcriptomic Reanalysis"* (Zirui Xu, Vanderbilt University).

The pipeline reanalyzes the public bulk RNA-seq dataset **GEO GSE79013** (Fang et al., *Sci Rep* 2016): left ventricles of one-year-old male CD-1 mice exposed to caffeine (20 mg/kg, i.p., E6.5–9.5) or vehicle *in utero*, n = 3 per group. One command downloads the data and reproduces every table and figure in the paper.

## What it does

| Step | Tool | Output |
|---|---|---|
| 1. Sample metadata (GEO + SRA run info) | curl, NCBI E-utilities / ENA | `metadata/samples.tsv` |
| 2. Reference index (Ensembl 110, GRCm39 cDNA) | kallisto 0.46.1 | `ref/` |
| 3. Download + quantification (single-end, 100 bootstraps) | SRA Toolkit, kallisto | `quant/` |
| 4. Differential expression (gene level; Wald and likelihood ratio tests) | sleuth 0.30.1 | `results/TableS2_*.csv`, Tables 1–2 |
| 5. Pathway analysis (preranked GSEA; Hallmark, Reactome, BioCarta, WikiPathways) | fgsea 1.28, msigdbr 7.5.1 | `results/TableS3_*.csv` |
| 6. Figures 1–3 (PDF, 600-dpi TIFF, PNG) | ggplot2, patchwork | `figures/` |

`results/SUMMARY.txt` records software versions, mapping rates, gene counts and pathway directions for each run.

## Requirements

- macOS or Linux, ~25 GB free disk space, ~8 GB RAM, internet access
- [conda](https://github.com/conda-forge/miniforge) (Miniforge recommended)

## Run

```bash
# one-time setup
conda env create -f environment.yml
# on Apple Silicon Macs: CONDA_SUBDIR=osx-64 conda env create -f environment.yml
#                        (kallisto 0.46.1 is only built for x86_64; requires Rosetta)
conda activate caffeine-rnaseq

# full analysis (~1–2 h, mostly download time)
bash run_all.sh
```

The script is restartable: completed steps are skipped. Options can be set as environment variables, e.g. `THREADS=8 BOOTSTRAPS=100 ENSEMBL_RELEASE=110 bash run_all.sh`.

To redo only the statistics and figures after quantification: `Rscript analysis.R`.

## Statistical notes

- sleuth β is a natural-log fold change; tables report log2FC = β / ln 2.
- Primary differential-expression calls use the Wald test (FDR ≤ 0.05). Likelihood ratio test results are reported alongside as a more conservative check; with n = 3 per group, no gene reaches LRT FDR ≤ 0.05.
- GSEA ranks genes by −log10(Wald P) × sign(β), so a positive NES means enrichment among genes upregulated in the caffeine group. The sign convention was verified on simulated data with known up- and downregulated pathways.

## Data source

Fang X, et al. In utero caffeine exposure induces transgenerational effects on the adult heart. *Sci Rep* 2016;6:34106. https://doi.org/10.1038/srep34106 — GEO: [GSE79013](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79013)

## License

MIT
