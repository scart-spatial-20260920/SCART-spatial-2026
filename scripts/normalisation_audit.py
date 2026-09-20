# -*- coding: utf-8 -*-
"""
normalisation_audit.py

What did the Fig 2F 'area normalisation' actually do, and does the 2026-08-07
conclusion survive if we ignore it entirely?

The Fig 2F legend states that the three non-SBRT sections were partitioned using
the SBRT STV:TTV:Rim area proportions, 'producing identical zone area ratios in
all four arms'.  The h5ad keeps BOTH states:
    P_group_all / P_group_all_HML  -> before normalisation
    phase2_group / phase2_zone     -> after normalisation
so the claim is directly testable, and the whole analysis can be re-run on the
raw sections as a robustness check.
"""
import os
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc

DATA    = Path(os.environ.get('SCART_DATA', './data'))
RESULTS = Path(os.environ.get('RESULTS', './results'))
H5AD = DATA / 'SCART_spatial_annotated.h5ad'   # name of the GEO-deposited AnnData
if not H5AD.exists():
    H5AD = DATA / 'data_labeled.h5ad'          # original working-file name
OUT = RESULTS / 'tables'
OUT.mkdir(parents=True, exist_ok=True)

a = sc.read_h5ad(H5AD)
obs = a.obs.copy()
ARMS = ['control', 'SBRT', 'SCART1', 'SCART3']
ZON = ['H', 'M', 'L']
imm = obs['GPT_PTPRC_and_PTPRCAP'].astype(str) != 'Undefined'

# ---- decode the pre-normalisation labels ---------------------------------
pre_arm_map = {'control_D20': 'control', 'SBRT': 'SBRT',
               'SCART_D21': 'SCART1', 'SCART3_D1': 'SCART3'}
obs['pre_arm'] = obs['P_group_all'].astype(str).map(pre_arm_map)
obs['pre_zone'] = obs['P_group_all_HML'].astype(str).str.rsplit('_', n=1).str[-1]
obs = obs[obs['pre_zone'].isin(ZON)]

pre_t = pd.crosstab(obs['pre_arm'], obs['pre_zone']).reindex(index=ARMS, columns=ZON).fillna(0).astype(int)
pre_i = pd.crosstab(obs.loc[imm.reindex(obs.index, fill_value=False), 'pre_arm'],
                    obs.loc[imm.reindex(obs.index, fill_value=False), 'pre_zone'])
pre_i = pre_i.reindex(index=ARMS, columns=ZON).fillna(0).astype(int)

post_t = pd.crosstab(a.obs['phase2_group'], a.obs['phase2_zone']).reindex(index=ARMS, columns=ZON).fillna(0).astype(int)
post_i = pd.crosstab(a.obs.loc[imm, 'phase2_group'], a.obs.loc[imm, 'phase2_zone'])
post_i = post_i.reindex(index=ARMS, columns=ZON).fillna(0).astype(int)


def ratio_line(row):
    s = row.sum()
    return ' : '.join(f'{100*row[z]/s:5.1f}' for z in ZON)


print('=' * 100)
print('1.  DID THE NORMALISATION ACHIEVE WHAT THE Fig 2F LEGEND CLAIMS?')
print('=' * 100)
print('\n   target = the SBRT ratio, which the legend says all four arms were made to match\n')
print(f'{"arm":9s} {"BEFORE  STV:TTV:Rim":>26s} {"AFTER   STV:TTV:Rim":>26s} {"spots kept":>22s}')
for g in ARMS:
    kept = post_t.loc[g].sum(); had = pre_t.loc[g].sum()
    print(f'{g:9s} {ratio_line(pre_t.loc[g]):>26s} {ratio_line(post_t.loc[g]):>26s} '
          f'{kept:>8,d}/{had:>7,d} ({100*kept/had:4.1f}%)')
print(f'\n   TARGET (SBRT)                                {ratio_line(post_t.loc["SBRT"]):>26s}')
print(f'\n   -> discarded overall: {pre_t.values.sum()-post_t.values.sum():,} of '
      f'{pre_t.values.sum():,} spots ({100*(1-post_t.values.sum()/pre_t.values.sum()):.1f}%)')

print('\n   spots discarded per arm x zone:')
print((pre_t - post_t).to_string())

print('\n' + '=' * 100)
print('2.  ROBUSTNESS: the 2026-08-07 conclusion computed BOTH ways')
print('=' * 100)
rows = []
for lbl, T, I in [('AFTER normalisation (as published)', post_t, post_i),
                  ('BEFORE normalisation (raw sections)', pre_t, pre_i)]:
    print(f'\n--- {lbl}')
    print(f'{"arm":9s} {"Rim% of section":>16s} {"immune in Rim%":>15s} {"enrichment":>11s} '
          f'{"Rim immune / 1000 Rim spots":>29s}')
    for g in ARMS:
        rt, tt = T.loc[g, 'L'], T.loc[g].sum()
        ri, ti = I.loc[g, 'L'], I.loc[g].sum()
        enr = (ri / ti) / (rt / tt)
        dens = 1000 * ri / rt
        print(f'{g:9s} {100*rt/tt:15.1f}% {100*ri/ti:14.1f}% {enr:11.2f} {dens:29.1f}')
        rows.append(dict(state=lbl, arm=g, rim_spots=rt, section_spots=tt,
                         rim_immune=ri, immune_total=ti,
                         rim_pct=100*rt/tt, immune_in_rim_pct=100*ri/ti,
                         enrichment=enr, rim_density_per1000=dens))

pd.DataFrame(rows).to_csv(OUT / 'normalisation_audit.csv', index=False,
                          encoding='utf-8-sig', float_format='%.3f')
pre_t.to_csv(OUT / 'zone_spots_BEFORE_normalisation.csv', encoding='utf-8-sig')
print(f'\n[wrote] normalisation_audit.csv, zone_spots_BEFORE_normalisation.csv')

print('\n' + '=' * 100)
print('3.  READING')
print('=' * 100)
print("""
   enrichment = (share of that arm's immune spots sitting in the Rim)
              / (share of that arm's tissue that IS Rim)
   = 1.00 means the immune spots are spread exactly in proportion to how much Rim
   there is, i.e. no preference at all.

   The finding to check is whether SBRT is the only arm with genuine Rim
   over-representation.  If that holds in BOTH panels above, it is a property of
   the tissue and not an artefact of the normalisation step.
""")
