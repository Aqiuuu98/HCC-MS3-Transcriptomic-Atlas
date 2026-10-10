#!/usr/bin/env python3
# Step: Step32 gene ID harmonization
# Conda environment: spatial32_wsl
# Python version: 3.10.21
# Frozen dependencies: /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_environment_freeze.txt (SHA256 6fe6d67eeb0584744c98b7242f6b56375d245b03808b3e19f78257181542c395)
# Relevant package versions: pandas==2.3.3; scanpy==1.11.5; anndata==0.11.4
# Executed source: /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_GeneCoverageAudit.py
# Executed source SHA256: 9cf770cb547d944df538b8722d1532b72f4be53dbd8ac40a6254e9d6455eceee
# Inputs:
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32D_FINAL_LOCK.json
#   /mnt/e/TCGA/LIHC_scRNA_Step25_CrossPatient_Robustness/Reference/Step25_FROZEN_MS3_Conserved_Core50.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step31_Independent_Validation/Audit/Step31_Core50_FrozenValidationMapping.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step31_Independent_Validation/Audit/Step31E2A_ValidationProgramMapping.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Data/QC/GSE281759/*_raw_qc.h5ad
# Outputs:
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_GeneIdentifierCoverage.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_GeneProgramCoverageSummary.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_GeneCoverageAudit.json
# Source body below is the executed version; only this header and final REPORT were added.
"""Audit frozen Step25/31 gene identifiers against the Step32 Visium universe."""
from hcc_release_paths import project_path as _hcc_path

import json
from pathlib import Path

import pandas as pd
import scanpy as sc


ROOT = Path(_hcc_path('LIHC_scRNA_Step32_Spatial_Validation'))
AUDIT = ROOT / 'Audit'
STEP25 = Path(_hcc_path('LIHC_scRNA_Step25_CrossPatient_Robustness/Reference'))
STEP31 = Path(_hcc_path('LIHC_scRNA_Step31_Independent_Validation/Audit'))


def main():
    if json.loads((AUDIT / 'Step32D_FINAL_LOCK.json').read_text())['status'] != 'PASS':
        raise ValueError('Step32D must be frozen before gene coverage audit')
    frozen = pd.read_csv(STEP25 / 'Step25_FROZEN_MS3_Conserved_Core50.csv')
    core = pd.read_csv(STEP31 / 'Step31_Core50_FrozenValidationMapping.csv')
    programs = pd.read_csv(STEP31 / 'Step31E2A_ValidationProgramMapping.csv')
    if len(frozen) != 50 or len(core) != 50:
        raise ValueError('Frozen Core50 row count mismatch')
    if frozen['Gene'].tolist() != core['canonical_symbol'].tolist():
        raise ValueError('Step25/Step31 frozen Core50 membership/order mismatch')
    ids = []
    for path in sorted((ROOT / 'Data/QC/GSE281759').glob('*_raw_qc.h5ad')):
        adata = sc.read_h5ad(path, backed='r')
        ids.append(set(adata.var_names.astype(str)))
        adata.file.close()
    if len(ids) != 17 or len({len(x) for x in ids}) != 1:
        raise ValueError('Section gene universes differ in size')
    shared = set.intersection(*ids)
    union = set.union(*ids)
    if shared != union:
        raise ValueError('Section gene IDs differ; review before scoring')
    rows = []
    for _, row in core.iterrows():
        gene_id = str(row['ensembl_ids'])
        rows.append({'program': 'Core50', 'frozen_symbol': row['canonical_symbol'],
                     'validation_symbol': row['validation_symbol'],
                     'ensembl_id': gene_id, 'step31_status': 'frozen_core50',
                     'present_in_spatial': gene_id in shared})
    for _, row in programs.loc[programs['mapping_status'] == 'retained'].iterrows():
        gene_id = str(row['selected_ensembl_id'])
        rows.append({'program': row['program'], 'frozen_symbol': row['requested_symbol'],
                     'validation_symbol': row['requested_symbol'],
                     'ensembl_id': gene_id, 'step31_status': 'retained',
                     'present_in_spatial': gene_id in shared})
    coverage = pd.DataFrame(rows)
    summary = (coverage.groupby('program').agg(
        n_frozen_genes=('ensembl_id', 'size'),
        n_present=('present_in_spatial', 'sum')).reset_index())
    summary['coverage_fraction'] = summary['n_present'] / summary['n_frozen_genes']
    coverage.to_csv(AUDIT / 'Step32_GeneIdentifierCoverage.csv', index=False)
    summary.to_csv(AUDIT / 'Step32_GeneProgramCoverageSummary.csv', index=False)
    report = {'n_sections': len(ids), 'n_spatial_gene_ids': len(shared),
              'core50_present': int(summary.loc[summary.program == 'Core50', 'n_present'].iloc[0]),
              'core50_frozen': 50, 'programs': summary.to_dict(orient='records'),
              'source_core50': str(STEP31 / 'Step31_Core50_FrozenValidationMapping.csv'),
              'source_programs': str(STEP31 / 'Step31E2A_ValidationProgramMapping.csv')}
    (AUDIT / 'Step32_GeneCoverageAudit.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()

if __name__ == '__main__':
    print('REPORT: Frozen Core50 and program Ensembl coverage; see outputs listed in header.')
