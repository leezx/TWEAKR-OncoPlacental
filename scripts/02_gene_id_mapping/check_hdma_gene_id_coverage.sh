#!/bin/bash
# Coverage check: how much of HDMA's Ensembl-fallback gene IDs actually
# resolve via the HGNC mapping table. Run on Argos (needs the argos-codex
# env for both R/SeuratObject and Python).
#
# Not a qsub job — Adrenal is small enough (readRDS ~7-8s) to run
# interactively on the login node.
set -euo pipefail

source /home/zz950/softwares/miniforge3/etc/profile.d/conda.sh
conda activate argos-codex

DATA=/home/zz950/DATA
MAP="$DATA/1.Databases/HGNC_gene_id_mapping"

# 1. Dump HDMA_Adrenal's gene names, split into ENSG-format vs symbol-format
Rscript -e '
suppressMessages(library(SeuratObject))
obj <- readRDS("'"$DATA"'/scRNAseq/HumanDevelopmentMultiomicAtlas/raw/per_organ_RNA_seurat/Adrenal_RNA_obj_clustered_final.rds")
g <- rownames(obj)
n_ensg <- sum(grepl("^ENSG", g))
cat("total genes:", length(g), "\n")
cat("ENSG-format:", n_ensg, "\n")
cat("symbol-format:", length(g) - n_ensg, "\n")
writeLines(g[grepl("^ENSG", g)], "/tmp/adrenal_ensg_ids.txt")
'

# 2. Check how many of those ENSG IDs resolve via the HGNC mapping
python3 -c "
import csv
ensg = set(l.strip() for l in open('/tmp/adrenal_ensg_ids.txt'))

# any-status check against the raw dump (both Ensembl-ID columns)
in_hgnc_at_all = set()
with open('$MAP/raw/hgnc_custom_download.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        for col in ('Ensembl gene ID', 'Ensembl ID(supplied by Ensembl)'):
            e = row.get(col, '')
            if e:
                in_hgnc_at_all.add(e)

found = [e for e in ensg if e in in_hgnc_at_all]
missing = [e for e in ensg if e not in in_hgnc_at_all]
print('total ENSG to resolve:', len(ensg))
print('resolved (any HGNC status):', len(found))
print('genuinely absent from HGNC:', len(missing))
"
