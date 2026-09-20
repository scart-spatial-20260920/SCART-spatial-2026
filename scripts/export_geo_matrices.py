# -*- coding: utf-8 -*-
"""
export_geo_matrices.py

Split data_labeled.h5ad into the per-sample processed files GEO expects, plus an
annotated series-level AnnData, plus an MD5 manifest.

GEO wants processed data as RAW COUNTS, not log-normalised values.  The h5ad
carries log1p-normalised values in .X and the integer counts in .layers['counts'],
so the export reads the layer and refuses to run if it does not look like counts.

Usage
-----
    python export_geo_matrices.py \
        --h5ad "$SCART_DATA/data_labeled.h5ad" \
        --out  "SCART_GEO_submission"

Output (flat, no subdirectories - GEO does not accept nested folders)
    <ARM>_barcodes.tsv.gz
    <ARM>_features.tsv.gz
    <ARM>_matrix.mtx.gz
    <ARM>_tissue_positions.csv.gz
    SCART_spatial_annotated.h5ad
    md5sums.txt
"""
import argparse, gzip, hashlib, shutil, sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite
import scanpy as sc

ARMS = ['control', 'SBRT', 'SCART1', 'SCART3']
ARM_FILE = {'control': 'Control', 'SBRT': 'SBRT',
            'SCART1': 'SCARTx1', 'SCART3': 'SCARTx3'}
# mouse id per arm, as reported in the manuscript
ARM_MOUSE = {'control': 'D20', 'SBRT': 'D31', 'SCART1': 'D21', 'SCART3': 'D56'}
ARM_CHIP = {'control': 'Chip2', 'SBRT': 'Chip1', 'SCART1': 'Chip2', 'SCART3': 'Chip1'}
# annotations worth shipping so a reader can reproduce the figures
KEEP_OBS = ['phase2_group', 'phase2_zone', 'GPT_PTPRC_and_PTPRCAP', 'singleR_Ai',
            'P_group_all', 'P_group_all_HML', 'nCount_RNA', 'nFeature_RNA',
            'spatial_x', 'spatial_y']


def gzip_write(path: Path, lines):
    with gzip.open(path, 'wt', encoding='utf-8', newline='\n') as fh:
        for ln in lines:
            fh.write(ln + '\n')


def md5(path: Path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for blk in iter(lambda: fh.read(chunk), b''):
            h.update(blk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--h5ad', required=True)
    ap.add_argument('--out', default='SCART_GEO_submission')
    ap.add_argument('--layer', default='counts',
                    help="layer holding raw counts (default: counts)")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    print(f'[load] {args.h5ad}')
    ad = sc.read_h5ad(args.h5ad)
    print(f'  {ad.n_obs:,} spots x {ad.n_vars:,} genes')

    # --- pick the counts matrix and verify it really is counts ------------
    if args.layer in ad.layers:
        M = ad.layers[args.layer]
        src = f'layers["{args.layer}"]'
    else:
        M = ad.X
        src = '.X'
        print(f'  [warn] layer "{args.layer}" not found, falling back to .X')

    probe = M[:500]
    probe = probe.toarray() if sp.issparse(probe) else np.asarray(probe)
    if not np.allclose(probe, np.round(probe)) or probe.min() < 0:
        sys.exit(
            f'ABORT: {src} is not integer counts (min={probe.min():.3f}, '
            f'max={probe.max():.3f}).  GEO expects raw counts as processed data. '
            f'Pass --layer with the correct layer name.')
    print(f'  using {src} - verified integer, non-negative')

    if not sp.issparse(M):
        M = sp.csr_matrix(M)

    # --- per-arm export ---------------------------------------------------
    written = []
    for arm in ARMS:
        mask = (ad.obs['phase2_group'].astype(str) == arm).values
        n = int(mask.sum())
        if n == 0:
            print(f'  [skip] {arm}: no spots')
            continue
        tag = ARM_FILE[arm]
        print(f'[{tag}] {n:,} spots  (mouse {ARM_MOUSE[arm]}, {ARM_CHIP[arm]})')

        sub = M[mask]
        bcs = ad.obs_names[mask].astype(str).tolist()

        p = out / f'{tag}_barcodes.tsv.gz'
        gzip_write(p, bcs); written.append(p)

        p = out / f'{tag}_features.tsv.gz'
        # 10x 3-column convention: id, symbol, type.  No separate gene ids here,
        # so the symbol is used for both and that is stated in the README.
        gzip_write(p, [f'{g}\t{g}\tGene Expression' for g in ad.var_names.astype(str)])
        written.append(p)

        # mmwrite wants genes x spots (10x orientation)
        p_mtx = out / f'{tag}_matrix.mtx'
        mmwrite(str(p_mtx), sub.T.tocoo().astype(np.int32), field='integer')
        with open(p_mtx, 'rb') as fi, gzip.open(str(p_mtx) + '.gz', 'wb') as fo:
            shutil.copyfileobj(fi, fo)
        p_mtx.unlink()
        written.append(Path(str(p_mtx) + '.gz'))

        xy = ad.obsm['spatial'][mask]
        pos = pd.DataFrame({'barcode': bcs,
                            'in_tissue': 1,
                            'x': xy[:, 0], 'y': xy[:, 1],
                            'zone': ad.obs['phase2_zone'].astype(str).values[mask],
                            'mouse': ARM_MOUSE[arm], 'chip': ARM_CHIP[arm]})
        p = out / f'{tag}_tissue_positions.csv.gz'
        pos.to_csv(p, index=False, compression='gzip'); written.append(p)

    # --- series-level annotated object ------------------------------------
    print('[series] SCART_spatial_annotated.h5ad')
    keep = [c for c in KEEP_OBS if c in ad.obs.columns]
    slim = sc.AnnData(X=M, obs=ad.obs[keep].copy(), var=pd.DataFrame(index=ad.var_names))
    slim.obsm['spatial'] = ad.obsm['spatial']
    for k in ('X_umap', 'X_pca'):
        if k in ad.obsm:
            slim.obsm[k] = ad.obsm[k]
    slim.layers['counts'] = M
    p = out / 'SCART_spatial_annotated.h5ad'
    slim.write_h5ad(p, compression='gzip'); written.append(p)

    # --- manifest ---------------------------------------------------------
    print('[md5] hashing')
    with open(out / 'md5sums.txt', 'w', encoding='utf-8', newline='\n') as fh:
        for p in written:
            fh.write(f'{md5(p)}  {p.name}\n')

    total = sum(p.stat().st_size for p in written) / 1e9
    print(f'\n[done] {len(written)} files, {total:.2f} GB -> {out.resolve()}')
    print('       md5sums.txt written; paste these into the GEO metadata sheet.')
    print('\nREMINDER: GEO also requires the raw FASTQ from SeekGene. '
          'This script cannot produce them.')


if __name__ == '__main__':
    main()
