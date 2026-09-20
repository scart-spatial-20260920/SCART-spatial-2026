# -*- coding: utf-8 -*-
"""AUTHORITATIVE zone denominators, recomputed directly from data_labeled.h5ad.

This is the source of truth for the STV/TTV/Rim area-normalisation dispute.
The 2026-08-07 memo back-computed these from rim_L_per1000.csv + composition_counts.csv;
here we read them straight off the spot table so there is nothing to argue about.
"""
import os
import numpy as np, pandas as pd, scanpy as sc

from pathlib import Path
DATA    = Path(os.environ.get('SCART_DATA', './data'))
RESULTS = Path(os.environ.get('RESULTS', './results'))
H5AD = DATA / 'SCART_spatial_annotated.h5ad'   # name of the GEO-deposited AnnData
if not H5AD.exists():
    H5AD = DATA / 'data_labeled.h5ad'          # original working-file name
OUT = RESULTS / 'tables'
OUT.mkdir(parents=True, exist_ok=True)

print('[load]', H5AD)
a = sc.read_h5ad(H5AD)
obs = a.obs
print(f'  spots = {a.n_obs:,}')

GROUPS = ['control', 'SBRT', 'SCART1', 'SCART3']
ZONE_NAME = {'H': 'STV (core)', 'M': 'TTV (junction)', 'L': 'Rim (periphery)', 'X': 'unzoned'}
fine = 'GPT_PTPRC_and_PTPRCAP'
imm = obs[fine].astype(str) != 'Undefined'

# ---- 1. spot counts per group x zone ------------------------------------
ct = pd.crosstab(obs['phase2_group'], obs['phase2_zone']).reindex(index=GROUPS)
ct = ct.reindex(columns=[c for c in ['H', 'M', 'L', 'X'] if c in ct.columns], fill_value=0)
print('\n=== TISSUE SPOTS per group x zone (all spots) ===')
print(ct.to_string())

imm_ct = pd.crosstab(obs.loc[imm, 'phase2_group'], obs.loc[imm, 'phase2_zone']).reindex(index=GROUPS)
imm_ct = imm_ct.reindex(columns=ct.columns, fill_value=0).fillna(0).astype(int)
print('\n=== IMMUNE (PTPRC+) SPOTS per group x zone ===')
print(imm_ct.to_string())

# ---- 2. the table that settles the dispute -------------------------------
rows = []
for g in GROUPS:
    tot = int(ct.loc[g].sum())
    zoned = int(ct.loc[g, [c for c in ['H', 'M', 'L'] if c in ct.columns]].sum())
    rim_t = int(ct.loc[g, 'L'])
    rim_i = int(imm_ct.loc[g, 'L'])
    imm_tot = int(imm_ct.loc[g].sum())
    rows.append(dict(
        group=g,
        STV_spots=int(ct.loc[g, 'H']), TTV_spots=int(ct.loc[g, 'M']), Rim_spots=rim_t,
        unzoned_X=int(ct.loc[g, 'X']) if 'X' in ct.columns else 0,
        section_total=tot, zoned_total=zoned,
        Rim_pct_of_section=100 * rim_t / tot,
        Rim_pct_of_zoned=100 * rim_t / zoned,
        Rim_immune=rim_i, immune_total=imm_tot,
        immune_density_per1000_Rim=1000 * rim_i / rim_t,
        pct_immune_in_Rim=100 * rim_i / imm_tot,
        enrichment_vs_section=(rim_i / imm_tot) / (rim_t / tot),
        enrichment_vs_zoned=(rim_i / imm_tot) / (rim_t / zoned),
    ))
T = pd.DataFrame(rows).set_index('group')
pd.set_option('display.width', 250)
print('\n=== AUTHORITATIVE RIM TABLE (recomputed from h5ad) ===')
print(T.to_string(float_format=lambda v: f'{v:.3f}'))

# ---- 3. STV:TTV:Rim ratio, the thing the figure legends claim is equal ----
print('\n=== STV : TTV : Rim  spot-share of ZONED tissue (%)  ===')
print('   figure legends for Fig 2F / 3A claim these are identical across arms')
for g in GROUPS:
    z = ct.loc[g, ['H', 'M', 'L']]
    s = z.sum()
    print(f'  {g:9s}  {100*z["H"]/s:5.1f} : {100*z["M"]/s:5.1f} : {100*z["L"]/s:5.1f}   (zoned n={s:,})')

# ---- 4. compare to the 2026-08-07 memo -----------------------------------
MEMO = {'control': (968, 5401, 48), 'SBRT': (2034, 7298, 230),
        'SCART1': (4477, 10060, 215), 'SCART3': (9310, 13718, 393)}
print('\n=== CROSS-CHECK vs 2026-08-07 memo (back-computed values) ===')
print(f'{"group":9s} {"Rim spots":>20s} {"section total":>24s} {"Rim immune":>20s}')
for g in GROUPS:
    m_rim, m_tot, m_imm = MEMO[g]
    a_rim, a_tot, a_imm = int(ct.loc[g, 'L']), int(ct.loc[g].sum()), int(imm_ct.loc[g, 'L'])
    f = lambda memo, act: f'{act:>7,d} vs {memo:>7,d} {"OK" if act == memo else "DIFF"}'
    print(f'{g:9s} {f(m_rim, a_rim):>20s} {f(m_tot, a_tot):>24s} {f(m_imm, a_imm):>20s}')

T.to_csv(OUT / 'zone_truth_table.csv', encoding='utf-8-sig')
ct.to_csv(OUT / 'zone_spot_counts.csv', encoding='utf-8-sig')
imm_ct.to_csv(OUT / 'zone_immune_counts.csv', encoding='utf-8-sig')
print('\n[wrote] zone_truth_table.csv, zone_spot_counts.csv, zone_immune_counts.csv')
