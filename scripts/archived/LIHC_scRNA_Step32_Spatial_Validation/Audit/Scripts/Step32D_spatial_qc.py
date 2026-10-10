#!/usr/bin/env python3
# Step: Step32D raw QC audit
# Conda environment: spatial32_wsl
# Python version: 3.10.21
# Frozen dependencies: /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_environment_freeze.txt (SHA256 6fe6d67eeb0584744c98b7242f6b56375d245b03808b3e19f78257181542c395)
# Relevant package versions: numpy==2.2.6; pandas==2.3.3; scanpy==1.11.5; scipy==1.15.3
# Executed source: /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32D_qc_audit.py
# Executed source SHA256: 1e93b820f70cfd7fb2bef6d3a3d584763ffac50d11a53e8c5bce09e6055789f2
# Inputs:
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32C2_FINAL_LOCK.json
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32C1_GSE281759_SampleStructureAudit.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Data/Standardized/GSE281759/GSM*/**
# Outputs:
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Tables/Step32D_QCSpotMetrics.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Tables/Step32D_QCSectionSummary.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Tables/Step32D_QCCanonicalPatientSummary.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32D_QCAudit.json
# Source body below is the executed version; only this header and final REPORT were added.
"""Unfiltered Step32D raw-count QC audit; no biological scores or spot removal."""
from hcc_release_paths import project_path as _hcc_path

import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


ROOT = Path(_hcc_path('LIHC_scRNA_Step32_Spatial_Validation'))
AUDIT = ROOT / 'Audit'
TABLES = ROOT / 'Tables'
STD = ROOT / 'Data/Standardized/GSE281759'
SAMPLES = AUDIT / 'Step32C1_GSE281759_SampleStructureAudit.csv'
C2_LOCK = AUDIT / 'Step32C2_FINAL_LOCK.json'


def quantiles(values, prefix):
    values = np.asarray(values, dtype=float)
    return {f'{prefix}_{name}': float(np.quantile(values, q)) if len(values) else np.nan
            for name, q in [('min', 0), ('p01', .01), ('p05', .05), ('p25', .25),
                            ('median', .5), ('p75', .75), ('p95', .95),
                            ('p99', .99), ('max', 1)]}


def main():
    if json.loads(C2_LOCK.read_text(encoding='utf-8'))['status'] != 'PASS':
        raise ValueError('Step32C2 is not PASS')
    TABLES.mkdir(parents=True, exist_ok=True)
    samples = list(csv.DictReader(SAMPLES.open(newline='', encoding='utf-8')))
    all_spots = []
    sections = []
    for sample in samples:
        gsm = sample['gsm']
        base = STD / gsm
        adata = sc.read_10x_mtx(base, var_names='gene_ids', make_unique=False,
                                gex_only=True)
        x = adata.X.tocsr()
        symbols = adata.var['gene_symbols'].fillna('').astype(str)
        mt = symbols.str.upper().str.startswith('MT-').to_numpy()
        total = np.asarray(x.sum(axis=1)).ravel().astype(float)
        n_genes = np.asarray(x.getnnz(axis=1)).ravel().astype(int)
        mt_total = np.asarray(x[:, mt].sum(axis=1)).ravel().astype(float)
        pct_mt = np.divide(mt_total * 100, total, out=np.zeros_like(total), where=total > 0)
        pos = pd.read_csv(base / 'spatial/tissue_positions_list.csv.gz',
                          header=None, compression='gzip',
                          names=['barcode', 'in_tissue', 'array_row', 'array_col',
                                 'pxl_row_in_fullres', 'pxl_col_in_fullres'])
        indexed = pos.set_index('barcode').reindex(adata.obs_names)
        if indexed['in_tissue'].isna().any():
            raise ValueError(f'Unmatched expression barcode in {gsm}')
        in_tissue = indexed['in_tissue'].to_numpy(dtype=int)
        if not np.isin(in_tissue, [0, 1]).all():
            raise ValueError(f'Invalid in_tissue flag in {gsm}')
        spots = pd.DataFrame({
            'gsm': gsm, 'barcode': adata.obs_names.astype(str),
            'raw_sample_label': sample['raw_sample_label'],
            'raw_patient_id': sample['raw_patient_id'],
            'canonical_patient_id': sample['canonical_patient_id'],
            'section_code': sample['section_code'],
            'total_counts': total, 'n_genes_by_counts': n_genes,
            'pct_counts_mt': pct_mt, 'in_tissue': in_tissue,
        })
        all_spots.append(spots)
        for scope, mask in [('all_expression', np.ones(len(spots), dtype=bool)),
                            ('in_tissue', in_tissue == 1)]:
            sub = spots.loc[mask]
            row = {'gsm': gsm, 'raw_sample_label': sample['raw_sample_label'],
                   'raw_patient_id': sample['raw_patient_id'],
                   'canonical_patient_id': sample['canonical_patient_id'],
                   'section_code': sample['section_code'], 'scope': scope,
                   'n_spots': len(sub), 'n_mt_genes': int(mt.sum()),
                   'n_zero_count': int((sub['total_counts'] == 0).sum()),
                   'n_under_200_genes': int((sub['n_genes_by_counts'] < 200).sum()),
                   'n_pct_mt_gt_20': int((sub['pct_counts_mt'] > 20).sum())}
            for metric in ['total_counts', 'n_genes_by_counts', 'pct_counts_mt']:
                row.update(quantiles(sub[metric], metric))
            sections.append(row)
        print(f'{gsm}: {len(spots)} expression, {int(in_tissue.sum())} in tissue; '
              f'median genes {int(np.median(n_genes[in_tissue == 1]))}', flush=True)

    spots = pd.concat(all_spots, ignore_index=True)
    section_df = pd.DataFrame(sections)
    patient_rows = []
    for patient, group in section_df.loc[section_df['scope'] == 'in_tissue'].groupby('canonical_patient_id'):
        patient_rows.append({
            'canonical_patient_id': patient, 'n_sections': len(group),
            'n_in_tissue_spots': int(group['n_spots'].sum()),
            'section_codes': '; '.join(group['section_code'].astype(str)),
            'median_section_total_counts': float(group['total_counts_median'].median()),
            'median_section_n_genes': float(group['n_genes_by_counts_median'].median()),
            'median_section_pct_mt': float(group['pct_counts_mt_median'].median()),
            'n_zero_count': int(group['n_zero_count'].sum()),
            'n_under_200_genes': int(group['n_under_200_genes'].sum()),
            'n_pct_mt_gt_20': int(group['n_pct_mt_gt_20'].sum()),
        })
    spots.to_csv(TABLES / 'Step32D_QCSpotMetrics.csv', index=False)
    section_df.to_csv(TABLES / 'Step32D_QCSectionSummary.csv', index=False)
    pd.DataFrame(patient_rows).to_csv(TABLES / 'Step32D_QCCanonicalPatientSummary.csv', index=False)
    report = {
        'stage': 'Step32D', 'status': 'QC_AUDIT_COMPLETE',
        'filter_applied': 'none', 'n_sections': len(samples),
        'n_canonical_patients': len(patient_rows),
        'n_expression_spots': len(spots),
        'n_in_tissue_expression_spots': int(spots['in_tissue'].sum()),
        'n_zero_count_in_tissue': int(((spots['in_tissue'] == 1) & (spots['total_counts'] == 0)).sum()),
        'n_under_200_genes_in_tissue': int(((spots['in_tissue'] == 1) & (spots['n_genes_by_counts'] < 200)).sum()),
        'n_pct_mt_gt_20_in_tissue': int(((spots['in_tissue'] == 1) & (spots['pct_counts_mt'] > 20)).sum()),
    }
    (AUDIT / 'Step32D_QCAudit.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()

if __name__ == '__main__':
    print('REPORT: Step32D unfiltered raw-count QC audit; see outputs listed in header.')
