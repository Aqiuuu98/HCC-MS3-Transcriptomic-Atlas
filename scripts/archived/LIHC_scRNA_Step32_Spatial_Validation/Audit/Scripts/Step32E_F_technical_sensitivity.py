#!/usr/bin/env python3
# Step: Step32E and Step32F technical sensitivity
# Conda environment: spatial32_wsl
# Python version: 3.10.21
# Frozen dependencies: /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_environment_freeze.txt (SHA256 6fe6d67eeb0584744c98b7242f6b56375d245b03808b3e19f78257181542c395)
# Relevant package versions: numpy==2.2.6; pandas==2.3.3; scipy==1.15.3
# Executed source: /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32EF_technical_sensitivity.py
# Executed source SHA256: 73866213e4a74cee9c1fe2255103768e1d91a9972314b7c666451f13e82997ee
# Inputs:
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Tables/Step32EF_SpotScores.csv
# Outputs:
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Tables/Step32EF_TechnicalSensitivityBySection.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Tables/Step32EF_TechnicalSensitivityByPatient.csv
# Source body below is the executed version; only this header and final REPORT were added.
"""Step31-style rank-residual technical sensitivity, section by section."""
from hcc_release_paths import project_path as _hcc_path

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr


ROOT = Path(_hcc_path('LIHC_scRNA_Step32_Spatial_Validation'))
TABLES = ROOT / 'Tables'
PROGRAM3 = ['G2M_Independent', 'Mitotic_Entry_Independent',
            'Spindle_Chromosome_Independent']
COVARIATES = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']


def residual(y, design):
    y_rank = rankdata(np.asarray(y, dtype=float), method='average')
    beta = np.linalg.lstsq(design, y_rank, rcond=None)[0]
    return y_rank - design @ beta


def main():
    scores = pd.read_csv(TABLES / 'Step32EF_SpotScores.csv')
    rows = []
    for gsm, group in scores.groupby('gsm', sort=True):
        cov = group[COVARIATES].to_numpy(dtype=float)
        design = np.column_stack([np.ones(len(group))] +
                                 [rankdata(cov[:, j], method='average')
                                  for j in range(cov.shape[1])])
        pairs = [('Core50', program) for program in PROGRAM3]
        pairs.append(('MS3_like', 'MS6_like'))
        for left, right in pairs:
            x_resid = residual(group[left], design)
            y_resid = residual(group[right], design)
            rows.append({'gsm': gsm,
                         'canonical_patient_id': group.canonical_patient_id.iloc[0],
                         'comparison': f'{left}_vs_{right}', 'n_spots': len(group),
                         'raw_rho': float(spearmanr(group[left], group[right]).statistic),
                         'technical_partial_rho': float(spearmanr(x_resid, y_resid).statistic),
                         'covariates': '; '.join(COVARIATES)})
    result = pd.DataFrame(rows)
    result.to_csv(TABLES / 'Step32EF_TechnicalSensitivityBySection.csv', index=False)
    summary = (result.groupby(['canonical_patient_id', 'comparison'])
               .agg(n_sections=('gsm', 'size'),
                    n_positive_partial=('technical_partial_rho', lambda x: int((x > 0).sum())),
                    n_negative_partial=('technical_partial_rho', lambda x: int((x < 0).sum())),
                    median_raw_rho=('raw_rho', 'median'),
                    median_partial_rho=('technical_partial_rho', 'median')).reset_index())
    summary.to_csv(TABLES / 'Step32EF_TechnicalSensitivityByPatient.csv', index=False)
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()

if __name__ == '__main__':
    print('REPORT: Step31-style rank-residual sensitivity; see outputs listed in header.')
