#!/usr/bin/env bash
# =============================================================================
# Reanalysis of GSE79013 (Fang et al. 2016): in utero caffeine vs vehicle,
# adult male F1 left ventricle, bulk RNA-seq.
#   1. sample metadata  (GEO + SRA run info)       -> metadata/
#   2. reference        (Ensembl GRCm39 cDNA)      -> ref/
#   3. quantification   (kallisto, 100 bootstraps) -> quant/
#   4. statistics, GSEA, figures, tables (R)       -> results/ figures/
# Every step is skipped if its output already exists, so the script can be
# re-run safely after an interruption.
#
# Usage:   conda activate caffeine-rnaseq
#          bash run_all.sh
# Options (environment variables):
#   THREADS=4  BOOTSTRAPS=100  ENSEMBL_RELEASE=110  FRAG_LEN=180  FRAG_SD=20
# =============================================================================
set -euo pipefail

GSE="GSE79013"
THREADS="${THREADS:-4}"
BOOTSTRAPS="${BOOTSTRAPS:-100}"
ENSEMBL_RELEASE="${ENSEMBL_RELEASE:-110}"
FRAG_LEN="${FRAG_LEN:-180}"
FRAG_SD="${FRAG_SD:-20}"
SEED=42
# Testing hooks (leave unset for the real run)
REF_FASTA="${REF_FASTA:-}"          # use a local cDNA FASTA instead of downloading
FASTQ_DIR="${FASTQ_DIR:-}"          # use local <SRR>.fastq(.gz) files instead of SRA

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
mkdir -p metadata ref quant results figures logs tmp
log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a logs/pipeline.log; }

# -----------------------------------------------------------------------------
# 1. Metadata
# -----------------------------------------------------------------------------
if [[ ! -s metadata/samples.tsv ]]; then
  log "Fetching GEO sample information for $GSE"
  curl -fsSL "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${GSE}&targ=gsm&form=text&view=brief" \
    -o metadata/geo_samples.soft
  # one row per GSM: gsm, title, SRX experiment accession
  awk -F' = ' '
    /^\^SAMPLE/        { if (g!="") print g "\t" t "\t" x; g=$2; t=""; x="" }
    /^!Sample_title/   { t=$2 }
    /^!Sample_relation = SRA:/ { n=split($2, a, "term="); x=a[n] }
    END                { if (g!="") print g "\t" t "\t" x }' <(tr -d '\r' < metadata/geo_samples.soft) > metadata/geo.tsv
  SRX_LIST=$(cut -f3 metadata/geo.tsv | grep -E '^SRX' | tr '\n' ' ')
  [[ -n "$SRX_LIST" ]] || { log "ERROR: no SRX accessions found in GEO record"; exit 1; }
  log "SRA experiments: $SRX_LIST"

  log "Fetching SRA run info (NCBI)"
  TERM_Q=$(echo $SRX_LIST | sed 's/ /%20OR%20/g')
  ES=$(curl -fsSL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=${TERM_Q}&retmax=100&usehistory=y" || true)
  WEBENV=$(echo "$ES" | sed -n 's:.*<WebEnv>\(.*\)</WebEnv>.*:\1:p')
  QKEY=$(echo "$ES" | sed -n 's:.*<QueryKey>\(.*\)</QueryKey>.*:\1:p')
  rm -f metadata/runs.tsv
  if [[ -n "$WEBENV" ]]; then
    curl -fsSL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&query_key=${QKEY}&WebEnv=${WEBENV}&rettype=runinfo&retmode=text" \
      -o metadata/runinfo.csv || true
    Rscript -e 'ri <- tryCatch(read.csv("metadata/runinfo.csv", stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(ri) && all(c("Run", "Experiment") %in% names(ri))) {
        ri <- ri[ri$Run != "Run" & ri$Run != "", ]
        write.table(data.frame(run = ri$Run, srx = ri$Experiment, layout = ri$LibraryLayout,
                               read_length = ri$avgLength, spots = ri$spots, instrument = ri$Model),
                    "metadata/runs.tsv", sep = "\t", quote = FALSE, row.names = FALSE) }'
  fi
  if [[ ! -s metadata/runs.tsv ]]; then
    log "NCBI run info unavailable, trying ENA"
    printf "run\tsrx\tlayout\tread_length\tspots\tinstrument\n" > metadata/runs.tsv
    for x in $SRX_LIST; do
      curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${x}&result=read_run&fields=run_accession,experiment_accession,library_layout,base_count,read_count,instrument_model&format=tsv" \
        | awk -F'\t' 'NR>1 { rl = ($5>0 ? int($4/$5) : "NA"); print $1 "\t" $2 "\t" $3 "\t" rl "\t" $5 "\t" $6 }' >> metadata/runs.tsv
    done
  fi

  Rscript - <<'REOF'
geo  <- read.delim("metadata/geo.tsv", header = FALSE, col.names = c("gsm", "title", "srx"), stringsAsFactors = FALSE)
runs <- read.delim("metadata/runs.tsv", stringsAsFactors = FALSE)
m <- merge(geo, runs, by = "srx")
if (nrow(m) == 0) { message("*** STOP: could not map GEO samples to SRA runs. Send a screenshot of the metadata folder to Claude."); quit(status = 2) }
m$condition <- ifelse(grepl("caff", m$title, ignore.case = TRUE), "caffeine",
               ifelse(grepl("veh|control|ctrl|nacl|saline", m$title, ignore.case = TRUE), "control", NA))
m <- m[order(m$condition, m$title), ]
m$sample <- ave(m$condition, m$condition, FUN = function(x) paste0(x, "_", seq_along(x)))
out <- m[, c("sample", "condition", "run", "gsm", "title", "srx", "layout", "read_length", "spots", "instrument")]
write.table(out, "metadata/samples.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(out, row.names = FALSE)
ok <- !any(is.na(out$condition)) && sum(out$condition == "caffeine") == 3 &&
      sum(out$condition == "control") == 3 && !any(duplicated(out$sample))
if (!ok) {
  file.rename("metadata/samples.tsv", "metadata/samples_CHECK.tsv")
  message("\n*** STOP: expected exactly 3 caffeine + 3 control runs (one run per sample).\n",
          "*** See metadata/samples_CHECK.tsv and send a screenshot to Claude.\n")
  quit(status = 2)
}
REOF
fi
log "Samples:"; cat metadata/samples.tsv | tee -a logs/pipeline.log

# -----------------------------------------------------------------------------
# 2. Reference: Ensembl GRCm39 cDNA -> kallisto index + transcript-to-gene map
# -----------------------------------------------------------------------------
if [[ ! -s ref/transcriptome.idx ]]; then
  if [[ -z "$REF_FASTA" ]]; then
    REF_FASTA="ref/Mus_musculus.GRCm39.cdna.all.fa.gz"
    log "Downloading Ensembl release ${ENSEMBL_RELEASE} GRCm39 cDNA"
    curl -fsSL "https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz" \
      -o "$REF_FASTA"
  fi
  log "Building transcript-to-gene table"
  { echo -e "target_id\tens_gene\text_gene\tgene_biotype"
    (gzip -cdf "$REF_FASTA") | grep '^>' | sed 's/^>//' | awk '{
        tid=$1; g=""; s=""; b="";
        for (i=2;i<=NF;i++) {
          if ($i ~ /^gene:/)         { g=substr($i,6); sub(/\.[0-9]+$/,"",g) }
          if ($i ~ /^gene_symbol:/)  { s=substr($i,13) }
          if ($i ~ /^gene_biotype:/) { b=substr($i,14) }
        }
        if (s=="") s=g;
        print tid "\t" g "\t" s "\t" b }'
  } > ref/t2g.tsv
  log "Building kallisto index (needs ~4-8 GB RAM, ~10 min)"
  kallisto index -i ref/transcriptome.idx "$REF_FASTA" 2>&1 | tee -a logs/kallisto_index.log
  echo "ensembl_release=${ENSEMBL_RELEASE}" > ref/REFERENCE_INFO.txt
fi

# -----------------------------------------------------------------------------
# 3. Download reads + kallisto quant (one sample at a time; FASTQ deleted after)
# -----------------------------------------------------------------------------
tail -n +2 metadata/samples.tsv | while IFS=$'\t' read -r sample condition run rest; do
  if [[ -s "quant/${sample}/abundance.h5" ]]; then
    log "quant/${sample} exists, skipping"; continue
  fi
  log "=== ${sample} (${run})"
  if [[ -n "$FASTQ_DIR" ]]; then
    fq=$(ls "$FASTQ_DIR/${run}"*.f*q* | head -1)
  else
    prefetch --max-size 50G -O tmp "$run" < /dev/null >> logs/download.log 2>&1
    fasterq-dump --split-3 -e "$THREADS" -O tmp "tmp/${run}" < /dev/null >> logs/download.log 2>&1
    rm -rf "tmp/${run}"
    n=$(ls tmp/${run}*.fastq | wc -l | tr -d ' ')
    if [[ "$n" != "1" ]]; then
      log "ERROR: ${run} produced ${n} FASTQ files (expected single-end). Check metadata/samples.tsv LibraryLayout."; exit 1
    fi
    fq="tmp/${run}.fastq"
  fi
  kallisto quant -i ref/transcriptome.idx -o "quant/${sample}" --single -l "$FRAG_LEN" -s "$FRAG_SD" \
    -b "$BOOTSTRAPS" --seed "$SEED" -t "$THREADS" "$fq" < /dev/null 2>&1 | tee "logs/kallisto_${sample}.log"
  if [[ -z "$FASTQ_DIR" ]]; then rm -f "$fq"; fi
done

# -----------------------------------------------------------------------------
# 4. Statistics, GSEA, figures, tables
# -----------------------------------------------------------------------------
cat > results/run_parameters.txt <<EOF
GEO accession: ${GSE}
Ensembl release (GRCm39 cDNA): ${ENSEMBL_RELEASE}
kallisto: $(kallisto version)
kallisto quant: --single -l ${FRAG_LEN} -s ${FRAG_SD} -b ${BOOTSTRAPS} --seed ${SEED}
sra-tools: $(fasterq-dump --version 2>/dev/null | grep -Eo '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)
EOF
log "Running analysis.R"
Rscript analysis.R 2>&1 | tee logs/analysis.log
log "Done. See results/SUMMARY.txt, figures/, results/."
