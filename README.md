# SCART spatial transcriptomics (2026)

Analysis code for the study of **spatially fractionated radiotherapy (SCART)** and its
effect on the tumour immune microenvironment in murine H22 hepatocellular carcinoma,
profiled by SeekGene SeekSpace spatial transcriptomics.

> **Replace before publishing:** `<ORG>`, the GEO accession, and the Zenodo DOI.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

---

## What this repository is

Code only. It reproduces every figure and table in the manuscript from the data
deposited at GEO. **No data files are stored here** — see [Data](#data).

## Study design

Four subcutaneous H22 tumours, **one mouse per arm**, harvested on Day 3 immediately
after the final fraction:

| Arm | Mouse | Capture chip | Spots analysed |
|---|---|---|---|
| Control (untreated) | D20 | Chip 2 | 5,401 |
| SBRT (uniform dose) | D31 | Chip 1 | 7,298 |
| SCART ×1 | D21 | Chip 2 | 10,060 |
| SCART ×3 | D56 | Chip 1 | 13,718 |

58,160 spots passed quality control in total; 36,477 were assigned to a treatment arm.

Each spot carries one of three dose-defined zone labels derived from the CT-based dose
simulation: **STV** (core), **TTV** (junction), **Rim** (periphery).

> ⚠️ **n = 1 per arm.** This is an exploratory, descriptive study. Spots within a
> section are not independent replicates, so no between-animal inferential statistics
> are reported and none should be inferred from the code.

## Data

All data are at **NCBI GEO accession [GSEXXXXXX](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSEXXXXXX)**:

- raw sequencing reads (FASTQ)
- per-arm raw count matrices in 10x Matrix Market format, plus barcodes, features and
  spatial coordinates
- `SCART_spatial_annotated.h5ad` — a series-level AnnData with **all 58,160 spots**
  (including the 21,683 not assigned to an arm), raw counts in `layers["counts"]`,
  log-normalised values in `.X`, coordinates in `.obsm["spatial"]`, and both the
  pre- and post-normalisation zone labels in `.obs`

Download it and point `SCART_DATA` at the directory:

```bash
export SCART_DATA=/path/to/downloaded/data      # Windows: $env:SCART_DATA="D:\data"
```

Multiplex immunofluorescence whole-slide images are not deposited here because of
file size; they are available from the corresponding author on request.

## Environment

### Option A — conda

```bash
conda env create -f environment.yml
conda activate scart-spatial
```

### Option B — Docker (recommended for exact reproduction)

```bash
docker pull ghcr.io/<ORG>/scart-spatial:v1.0.0

docker run --rm \
  -v "$SCART_DATA":/data \
  -v "$PWD/results":/results \
  -e SCART_DATA=/data \
  ghcr.io/<ORG>/scart-spatial:v1.0.0 \
  python scripts/regenerate_fig4C.py
```

Versions in `environment.yml` are those actually present on the machine that produced
the manuscript figures, not minimum requirements.

## Running the analysis

```bash
make all          # full pipeline
make figures      # figures only, from cached intermediates
make fig4         # regenerate Figure 4 alone
```

Or step by step:

| Script | Produces |
|---|---|
| `phase2_03_scanpy_load.py` | QC, normalisation, integration → `data_labeled.h5ad` |
| `phase2_04_spatial_celltype.py` | PTPRC+ immune spots, 14-subtype annotation |
| `phase2_05_squidpy_neighborhood.py` | neighbourhood-enrichment z-scores (Figure 3A/5B) |
| `phase2_06_deg_rim.py` | SCART×3 Rim vs SBRT Rim differential expression (Figure 4A) |
| `phase2_07_pathway_decoupler.py` | Hallmark pathway scores (Figure 4C) |
| `phase2_08_morans_i.py` | Moran's I spatial autocorrelation (Figure S1) |
| `phase2_10_journal_figures.py` | Figures 1–4 at journal specification |
| `regenerate_fig4C.py` | Figure 4C and the complete Figure 4 |
| `zone_truth.py` | zone spot counts and Rim enrichment table |
| `normalisation_audit.py` | before/after comparison of the zone-normalisation step |
| `export_geo_matrices.py` | GEO submission files from the h5ad |

### A note on pathway scores

Figure 4C scores are computed with **Scanpy `sc.tl.score_genes`**, which measures
expression of a gene set relative to a matched background set. They are *not* AUCell
scores, despite the historical filename `pathway_aucell_per_cell.csv`. Zero is
meaningful and negative values are legitimate.

### A note on the zone normalisation

The published Figure 2F describes an area-normalisation step intended to give all four
arms identical STV:TTV:Rim ratios. `normalisation_audit.py` compares the zone ratios
before and after that step and reports what it actually achieved, together with the
number of spots it discarded. Both states are preserved in the deposited h5ad so the
comparison can be repeated independently. Conclusions based on per-1,000-Rim-spot
densities are unaffected either way, because those densities do not depend on how large
the Rim is.

## Citation

If you use this code, please cite the paper and the archived release:

```bibtex
@article{scart_spatial_2026,
  title   = {...},
  author  = {...},
  journal = {...},
  year    = {2026},
  doi     = {...}
}
```

See `CITATION.cff` for the software citation.

## License

MIT — see [LICENSE](LICENSE).
