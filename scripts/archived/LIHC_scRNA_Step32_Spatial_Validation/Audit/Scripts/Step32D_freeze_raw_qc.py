#!/usr/bin/env python3
# Step: Step32D H5AD freeze
# Conda environment: spatial32_wsl
# Python version: 3.10.21
# Frozen dependencies: /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_environment_freeze.txt (SHA256 6fe6d67eeb0584744c98b7242f6b56375d245b03808b3e19f78257181542c395)
# Relevant package versions: pandas==2.3.3; scanpy==1.11.5; anndata==0.11.4; h5py==3.16.0
# Executed source: /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32D_freeze.py
# Executed source SHA256: 73f9720c8cc904065e75490ab5732cd209ac21a8ef383a5bf800114891d7f189
# Inputs:
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32C2_FINAL_LOCK.json
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32D_QCAudit.json
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Tables/Step32D_QCSpotMetrics.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32C1_GSE281759_SampleStructureAudit.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Data/Standardized/GSE281759/GSM*/**
# Outputs:
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Data/QC/GSE281759/*_raw_qc.h5ad
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32D_FINAL_LOCK.json
# Source body below is the executed version; only this header and final REPORT were added.
"""Freeze unfiltered, in-tissue Step32 raw-count universe as one H5AD per GSM."""
from hcc_release_paths import project_path as _hcc_path

import csv
import json
import os
from pathlib import Path

import pandas as pd
import scanpy as sc


ROOT = Path(_hcc_path('LIHC_scRNA_Step32_Spatial_Validation'))
AUDIT = ROOT / 'Audit'
STD = ROOT / 'Data/Standardized/GSE281759'
DEST = ROOT / 'Data/QC/GSE281759'
SPOTS = ROOT / 'Tables/Step32D_QCSpotMetrics.csv'
SAMPLES = AUDIT / 'Step32C1_GSE281759_SampleStructureAudit.csv'


def main():
    c2 = json.loads((AUDIT / 'Step32C2_FINAL_LOCK.json').read_text(encoding='utf-8'))
    qc = json.loads((AUDIT / 'Step32D_QCAudit.json').read_text(encoding='utf-8'))
    if c2['status'] != 'PASS' or qc['status'] != 'QC_AUDIT_COMPLETE':
        raise ValueError('C2 or D audit prerequisite missing')
    samples = list(csv.DictReader(SAMPLES.open(newline='', encoding='utf-8')))
    metrics = pd.read_csv(SPOTS, dtype={'gsm': str, 'barcode': str})
    DEST.mkdir(parents=True, exist_ok=True)
    records = []
    for sample in samples:
        gsm = sample['gsm']
        base = STD / gsm
        target = DEST / f'{gsm}_raw_qc.h5ad'
        sub = metrics.loc[metrics['gsm'] == gsm].set_index('barcode')
        if not (sub['in_tissue'] == 1).all():
            raise ValueError(f'Non-tissue expression spot in {gsm}; freeze rule requires review')
        if target.exists():
            existing = sc.read_h5ad(target, backed='r')
            if existing.n_obs != len(sub) or existing.n_vars != c2.get('n_genes', 36601):
                existing.file.close()
                raise ValueError(f'Existing H5AD shape mismatch; no overwrite: {target}')
            existing.file.close()
            action = 'existing_shape_match'
        else:
            adata = sc.read_10x_mtx(base, var_names='gene_ids', make_unique=False,
                                    gex_only=True)
            if not adata.obs_names.equals(sub.index):
                sub = sub.reindex(adata.obs_names)
            if sub['total_counts'].isna().any():
                raise ValueError(f'QC metrics do not align in {gsm}')
            for col in ['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'in_tissue',
                        'raw_sample_label', 'raw_patient_id', 'canonical_patient_id',
                        'section_code']:
                adata.obs[col] = sub[col].to_numpy()
            pos = pd.read_csv(base / 'spatial/tissue_positions_list.csv.gz',
                              header=None, compression='gzip',
                              names=['barcode', 'in_tissue', 'array_row', 'array_col',
                                     'pxl_row_in_fullres', 'pxl_col_in_fullres'])
            pos = pos.set_index('barcode').reindex(adata.obs_names)
            if pos['pxl_row_in_fullres'].isna().any():
                raise ValueError(f'Spatial positions do not align in {gsm}')
            adata.obs['array_row'] = pos['array_row'].to_numpy()
            adata.obs['array_col'] = pos['array_col'].to_numpy()
            adata.obsm['spatial'] = pos[['pxl_col_in_fullres', 'pxl_row_in_fullres']].to_numpy()
            adata.var['mt'] = adata.var['gene_symbols'].fillna('').str.upper().str.startswith('MT-').to_numpy()
            adata.uns['scalefactors'] = json.loads((base / 'spatial/scalefactors_json.json').read_text(encoding='utf-8'))
            adata.uns['raw_qc_universe'] = ('All filtered-expression matrix spots; all have '
                                            'in_tissue=1; no additional spot filter; X is raw counts')
            adata.uns['hires_image_path'] = str(base / 'spatial/tissue_hires_image.png')
            adata.uns['lowres_image_path'] = str(base / 'spatial/tissue_lowres_image.png')
            temp = target.with_name(target.name + f'.step32d-{os.getpid()}.tmp.h5ad')
            try:
                adata.write_h5ad(temp, compression='lzf')
                if target.exists():
                    raise ValueError(f'Output appeared during write; no overwrite: {target}')
                os.replace(temp, target)
            finally:
                if temp.exists():
                    temp.unlink()
            action = 'created'
        records.append({'gsm': gsm, 'canonical_patient_id': sample['canonical_patient_id'],
                        'path': str(target), 'n_spots': len(sub),
                        'size_bytes': target.stat().st_size, 'action': action})
        print(f'{gsm}: {action} {len(sub)} spots', flush=True)
    lock = {'stage': 'Step32D', 'status': 'PASS', 'filter_rule': 'in_tissue == 1 only',
            'additional_qc_thresholds': None, 'n_sections': len(records),
            'n_canonical_patients': len({r['canonical_patient_id'] for r in records}),
            'n_spots_retained': sum(r['n_spots'] for r in records),
            'fraction_retained_of_expression': sum(r['n_spots'] for r in records) / qc['n_expression_spots'],
            'h5ad_files': records}
    (AUDIT / 'Step32D_FINAL_LOCK.json').write_text(json.dumps(lock, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in lock.items() if k != 'h5ad_files'}, indent=2), flush=True)


if __name__ == '__main__':
    main()

if __name__ == '__main__':
    print('REPORT: Step32D 17-section raw/QC H5AD lock; see outputs listed in header.')
