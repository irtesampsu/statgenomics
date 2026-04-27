# Project Structure and Run Guide

## Top-level files and folders

- `encode ids.txt`
  - ENCODE file ID mapping used to identify correct inputs.
- `download_files.R`
  - Utility script to download project data from ENCODE (if re-fetching files).
- `data/`
  - Input data used by analysis scripts.
- `scripts/`
  - Main analysis code (RNA, ATAC, integration, shared utilities).
- `results/`
  - Generated outputs (tables, plots, and report-ready summaries).
  
## Data layout (`data/`)

- `data/rna_seq/`
  - ScriptSeq RNA expression tables (`*.tsv`) used for RNA analyses.
- `data/atac_seq/`
  - ATAC peak files (`*.bigBed`) used for accessibility analyses.
- `data/bam/`
  - BAM files (available for optional stronger counting-based workflows).
- `data/manifests/`
  - Manifest and metadata files documenting sample/file mappings.

## Scripts layout (`scripts/`)

- `scripts/utils/project_metadata.R`
  - Central sample metadata and helper functions (`get_pairwise_*`, `get_all_*`) with input validation.
- `scripts/utils/rna_utils.R`
  - RNA helper functions: count matrix assembly, filtering, ID helpers, heatmap helpers.
- `scripts/utils/atac_utils.R`
  - ATAC helper functions: peak import, consensus construction, quantification helpers.

- `scripts/rna_seq/deseq2_pairwise_cmp_vs_ery.R`
  - Pairwise RNA differential expression (CMP vs Erythroblast) with DESeq2 + QC/plots + GO.
- `scripts/rna_seq/limma_voom_pairwise_cmp_vs_ery.R`
  - Same pairwise RNA analysis using limma-voom + QC/plots + GO.
- `scripts/rna_seq/summarize_pairwise_cmp_vs_ery.R`
  - Compares DESeq2 vs limma outputs and writes RNA summary artifacts.
- `scripts/rna_seq/heatmap_top1000_cmp_vs_ery.R`
  - Generates top-variable/top-feature heatmap artifact for pairwise RNA view.
- `scripts/rna_seq/q4_rna_clustering_all_cells.R`
  - All-cell RNA clustering (Q4): correlation heatmap, PCA, hierarchical clustering, top-variable gene heatmap.
- `scripts/rna_seq/run_pairwise_rna_generic.R`
  - Parameterized RNA pairwise runner for additional cell-line pairs (default method decision support).

- `scripts/atac_seq/atac_pairwise_cmp_vs_ery.R`
  - Pairwise ATAC differential accessibility (CMP vs Erythroblast), annotation, GO, and QC plots.
- `scripts/atac_seq/postprocess_atac_pairwise_cmp_vs_ery.R`
  - Post-processing for stringent DAR subset and directional GO terms.
- `scripts/atac_seq/q8_atac_clustering_all_cells.R`
  - All-cell ATAC clustering (Q8): correlation/PCA/hclust/top-variable peaks.

- `scripts/integration/integrate_pairwise_cmp_vs_ery.R`
  - RNA-ATAC integration for CMP vs Erythroblast (direction-concordant overlap, integration GO, summary text).

## Results layout (`results/`)

- `results/deseq2/`
  - Pairwise RNA DESeq2 outputs (DE tables, GO, PCA/heatmaps/volcano/correlation).
- `results/limma_voom/`
  - Pairwise RNA limma-voom outputs mirroring DESeq2 structure.
- `results/atac_seq/pairwise_cmp_vs_ery/`
  - Pairwise ATAC outputs (DARs, annotations, GO, PCA/volcano/correlation, stringent post-processing outputs).
- `results/clustering/q4_rna_all_cells/`
  - Q4 all-cell RNA clustering deliverables and interpretation README.
- `results/clustering/q8_atac_all_cells/`
  - Q8 all-cell ATAC clustering deliverables and interpretation README.
- `results/integration/pairwise_cmp_vs_ery/`
  - Integrated RNA-ATAC tables and integration summary for current pair.
- `results/summary/`
  - Report-ready summaries/checklists/storyline/method-decision templates and comparison tables.

## Recommended run sequence

Run from repository root:

1. **RNA pairwise analyses (both methods, first pair requirement)**
   - `Rscript scripts/rna_seq/deseq2_pairwise_cmp_vs_ery.R`
   - `Rscript scripts/rna_seq/limma_voom_pairwise_cmp_vs_ery.R`

2. **RNA method comparison summary**
   - `Rscript scripts/rna_seq/summarize_pairwise_cmp_vs_ery.R`
   - Optional visualization: `Rscript scripts/rna_seq/heatmap_top1000_cmp_vs_ery.R`

3. **ATAC pairwise analysis**
   - `Rscript scripts/atac_seq/atac_pairwise_cmp_vs_ery.R`
   - `Rscript scripts/atac_seq/postprocess_atac_pairwise_cmp_vs_ery.R`

4. **All-cell clustering required questions**
   - `Rscript scripts/rna_seq/q4_rna_clustering_all_cells.R`
   - `Rscript scripts/atac_seq/q8_atac_clustering_all_cells.R`

5. **RNA-ATAC integration (current pair)**
   - `Rscript scripts/integration/integrate_pairwise_cmp_vs_ery.R`

6. **Prepare report artifacts**
   - Review `results/summary/`, `results/integration/pairwise_cmp_vs_ery/`, and clustering READMEs.

## Sequence for additional pairwise cell lines

After first-pair method comparison, use the chosen RNA method for additional pairs:

- Generic runner example:
  - `Rscript scripts/rna_seq/run_pairwise_rna_generic.R --cell_a HSC --cell_b CMP --method limma_voom --out_dir results/limma_voom/pairwise_hsc_vs_cmp`

Then apply corresponding ATAC pairwise and integration workflow for that pair, and reuse `results/summary/pairwise_summary_template.txt` for consistent reporting.
