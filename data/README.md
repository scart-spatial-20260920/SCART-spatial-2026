# Data directory

**This directory is intentionally empty.** No data files are tracked in this
repository — see `.gitignore`.

All data for this study are deposited at NCBI GEO under accession
**[GSEXXXXXX](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSEXXXXXX)**.

## What to download

| File | Contents |
|---|---|
| `SCART_spatial_annotated.h5ad` | Series-level AnnData, all 58,160 QC-passing spots. Raw counts in `layers["counts"]`, log-normalised values in `.X`, coordinates in `.obsm["spatial"]`, pre- and post-normalisation zone labels in `.obs`. **This one file is enough to reproduce every figure.** |
| `<ARM>_matrix.mtx.gz` | Per-arm raw integer counts, genes × spots |
| `<ARM>_barcodes.tsv.gz` | Spot barcodes |
| `<ARM>_features.tsv.gz` | Gene symbols (three-column 10x convention; the upstream pipeline did not retain Ensembl identifiers, so the symbol appears in both id columns) |
| `<ARM>_tissue_positions.csv.gz` | barcode, in_tissue, x, y, zone, mouse, chip |
| FASTQ | Raw sequencing reads, four samples |

`<ARM>` is one of `Control`, `SBRT`, `SCARTx1`, `SCARTx3`.

## Where to put it

Anywhere. Then point the environment variable at it:

```bash
export SCART_DATA=/path/to/downloaded/data      # Windows: $env:SCART_DATA="D:\data"
```

Every script and the `Makefile` read that variable. Nothing writes back into the
data directory.

## Verifying the download

Each GEO supplementary file has a published MD5. Check them before running anything:

```bash
md5sum -c md5sums.txt
```

```powershell
# Windows
Get-ChildItem -File | ForEach-Object {
  "{0}  {1}" -f (Get-FileHash $_ -Algorithm MD5).Hash.ToLower(), $_.Name
}
```
