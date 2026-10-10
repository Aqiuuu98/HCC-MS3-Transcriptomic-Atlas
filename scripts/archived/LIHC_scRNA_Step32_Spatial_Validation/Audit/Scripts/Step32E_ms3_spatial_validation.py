#!/usr/bin/env python3
# Step: Step32E primary MS3 spatial validation
# Conda environment: spatial32_wsl
# Python version: 3.10.21
# Frozen dependencies: /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_environment_freeze.txt (SHA256 6fe6d67eeb0584744c98b7242f6b56375d245b03808b3e19f78257181542c395)
# Relevant package versions: numpy==2.2.6; pandas==2.3.3; scanpy==1.11.5; scipy==1.15.3; scikit-learn==1.7.2; anndata==0.11.4
# Executed source: /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32EF_score.py
# Executed source SHA256: 49b24efcb9788b6a549b881cb3c365b6f6c4bd896ce3f5cff9d8cc86f897dcc1
# Inputs:
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32D_FINAL_LOCK.json
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_GeneCoverageAudit.json
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32C1_GSE281759_SampleStructureAudit.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step31_Independent_Validation/Audit/Step31_Core50_FrozenValidationMapping.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step31_Independent_Validation/Audit/Step31E2A_ValidationProgramMapping.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Data/QC/GSE281759/*_raw_qc.h5ad
# Outputs:
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Tables/Step32EF_SpotScores.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Tables/Step32E_Core50ProgramBySection.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Tables/Step32E_Core50ProgramByPatient.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Tables/Step32F_MS3MS6BySection.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Tables/Step32F_MS3MS6ByPatient.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Tables/Step32EF_SectionMedianScores.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Tables/Step32EF_PatientMedianScores.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Tables/Step32EF_PatientLevelAssociation.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Tables/Step32EF_LOPO.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_SelectedGeneDetection.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32EF_ScoringMethod.json
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32EF_FINAL_LOCK.json
# Source body below is the executed version; only this header and final REPORT were added.
"""Frozen Core50/program scores and section-aware spatial associations."""
from hcc_release_paths import project_path as _hcc_path

import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr
from sklearn.neighbors import NearestNeighbors


ROOT = Path(_hcc_path('LIHC_scRNA_Step32_Spatial_Validation'))
AUDIT = ROOT / 'Audit'
TABLES = ROOT / 'Tables'
STEP31 = Path(_hcc_path('LIHC_scRNA_Step31_Independent_Validation/Audit'))
SAMPLES = AUDIT / 'Step32C1_GSE281759_SampleStructureAudit.csv'
PROGRAM3 = ['G2M_Independent', 'Mitotic_Entry_Independent',
            'Spindle_Chromosome_Independent']
PROGRAM6 = ['Hepatic_Differentiation', 'Bile_Acid_Metabolism']


def rho(x, y):
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return float(spearmanr(x, y).statistic)


def main():
    if json.loads((AUDIT / 'Step32D_FINAL_LOCK.json').read_text())['status'] != 'PASS':
        raise ValueError('Step32D not frozen')
    coverage = json.loads((AUDIT / 'Step32_GeneCoverageAudit.json').read_text())
    if coverage['core50_present'] != 50 or any(
            row['n_present'] != row['n_frozen_genes'] for row in coverage['programs']):
        raise ValueError('Frozen program coverage differs from audit')
    samples = list(csv.DictReader(SAMPLES.open(newline='', encoding='utf-8')))
    core = pd.read_csv(STEP31 / 'Step31_Core50_FrozenValidationMapping.csv')
    mapping = pd.read_csv(STEP31 / 'Step31E2A_ValidationProgramMapping.csv')
    programs = {'Core50': core['ensembl_ids'].astype(str).tolist()}
    for name in PROGRAM3 + PROGRAM6:
        programs[name] = mapping.loc[(mapping.program == name) &
                                     (mapping.mapping_status == 'retained'),
                                     'selected_ensembl_id'].astype(str).tolist()
    if len(programs['Core50']) != 50 or any(len(programs[p]) < 5 for p in PROGRAM3 + PROGRAM6):
        raise ValueError('Frozen gene sets malformed')
    for name in PROGRAM3 + PROGRAM6:
        if set(programs[name]) & set(programs['Core50']):
            raise ValueError(f'Core50 overlap in {name}')
    all_ids = sorted(set().union(*map(set, programs.values())))
    chunks = []
    meta = []
    detection = []
    for sample in samples:
        gsm = sample['gsm']
        adata = sc.read_h5ad(ROOT / f'Data/QC/GSE281759/{gsm}_raw_qc.h5ad')
        if adata.n_obs == 0 or not (adata.obs['in_tissue'] == 1).all():
            raise ValueError(f'Unexpected QC universe in {gsm}')
        x = adata[:, all_ids].X.toarray().astype(np.float64)
        totals = adata.obs['total_counts'].to_numpy(dtype=float)
        if not np.allclose(np.asarray(adata.X.sum(axis=1)).ravel(), totals):
            raise ValueError(f'Raw counts changed in {gsm}')
        logexpr = np.log1p((x / totals[:, None]) * 1e4)
        chunks.append(logexpr)
        for j, gene_id in enumerate(all_ids):
            detection.append({'gsm': gsm, 'ensembl_id': gene_id,
                              'n_spots': adata.n_obs,
                              'n_detected': int((x[:, j] > 0).sum()),
                              'detection_fraction': float((x[:, j] > 0).mean())})
        meta.append(pd.DataFrame({
            'gsm': gsm, 'barcode': adata.obs_names.astype(str),
            'raw_sample_label': sample['raw_sample_label'],
            'raw_patient_id': sample['raw_patient_id'],
            'canonical_patient_id': sample['canonical_patient_id'],
            'section_code': sample['section_code'],
            'total_counts': totals,
            'n_genes_by_counts': adata.obs['n_genes_by_counts'].to_numpy(),
            'pct_counts_mt': adata.obs['pct_counts_mt'].to_numpy(),
            'spatial_x': adata.obsm['spatial'][:, 0],
            'spatial_y': adata.obsm['spatial'][:, 1],
        }))
        print(f'{gsm}: {adata.n_obs} spots; selected genes {len(all_ids)}', flush=True)
    raw = np.vstack(chunks)
    spots = pd.concat(meta, ignore_index=True)
    if len(spots) != 27224 or raw.shape != (27224, len(all_ids)):
        raise ValueError('Scoring universe does not match frozen QC universe')
    mu = raw.mean(axis=0)
    sd = raw.std(axis=0, ddof=0)
    z = np.divide(raw - mu, sd, out=np.zeros_like(raw), where=sd > 0)
    index = {gene: j for j, gene in enumerate(all_ids)}
    for program, ids in programs.items():
        cols = [index[gene] for gene in ids]
        spots[program] = z[:, cols].mean(axis=1)
        spots[program + '_detected_n'] = (raw[:, cols] > 0).sum(axis=1)
    spots['MS3_like'] = spots[PROGRAM3].mean(axis=1)
    spots['MS6_like'] = spots[PROGRAM6].mean(axis=1)
    if not np.isfinite(spots[list(programs) + ['MS3_like', 'MS6_like']].to_numpy()).all():
        raise ValueError('Nonfinite score')

    e_rows = []
    f_rows = []
    for sample in samples:
        sub = spots.loc[spots.gsm == sample['gsm']]
        for program in PROGRAM3:
            e_rows.append({'gsm': sample['gsm'],
                           'canonical_patient_id': sample['canonical_patient_id'],
                           'section_code': sample['section_code'], 'program': program,
                           'n_spots': len(sub),
                           'rho_core50_program': rho(sub.Core50, sub[program]),
                           'median_core50_detected_n': float(sub.Core50_detected_n.median())})
        coords = sub[['spatial_x', 'spatial_y']].to_numpy()
        if not np.isfinite(coords).all() or np.unique(coords, axis=0).shape[0] != len(coords):
            raise ValueError(f'Spatial coordinates invalid for neighbor audit: {sample["gsm"]}')
        neighbors = NearestNeighbors(n_neighbors=7, metric='euclidean').fit(coords)
        neighbor_index = neighbors.kneighbors(return_distance=False)[:, :6]
        adjacent_ms6 = sub.MS6_like.to_numpy()[neighbor_index].mean(axis=1)
        f_rows.append({'gsm': sample['gsm'],
                       'canonical_patient_id': sample['canonical_patient_id'],
                       'section_code': sample['section_code'], 'n_spots': len(sub),
                       'rho_within_spot_ms3_ms6': rho(sub.MS3_like, sub.MS6_like),
                       'rho_ms3_neighbor_ms6': rho(sub.MS3_like, adjacent_ms6),
                       'neighbor_k': 6, 'coordinate_space': 'fullres_pixel'})
    e_section = pd.DataFrame(e_rows)
    f_section = pd.DataFrame(f_rows)
    e_patient = (e_section.groupby(['canonical_patient_id', 'program'])
                 .agg(n_sections=('gsm', 'size'),
                      n_positive_sections=('rho_core50_program', lambda x: int((x > 0).sum())),
                      median_section_rho=('rho_core50_program', 'median'),
                      min_section_rho=('rho_core50_program', 'min'),
                      max_section_rho=('rho_core50_program', 'max')).reset_index())
    f_patient = (f_section.groupby('canonical_patient_id')
                 .agg(n_sections=('gsm', 'size'),
                      n_negative_sections=('rho_within_spot_ms3_ms6', lambda x: int((x < 0).sum())),
                      median_section_rho=('rho_within_spot_ms3_ms6', 'median'),
                      median_neighbor_rho=('rho_ms3_neighbor_ms6', 'median')).reset_index())
    section_medians = spots.groupby(['canonical_patient_id', 'gsm'])[
        list(programs) + ['MS3_like', 'MS6_like']].median().reset_index()
    patient_medians = section_medians.groupby('canonical_patient_id')[
        list(programs) + ['MS3_like', 'MS6_like']].median().reset_index()
    patients = patient_medians['canonical_patient_id'].tolist()
    agg_rows = []
    for name in PROGRAM3:
        agg_rows.append({'comparison': f'Core50_vs_{name}', 'n_patients': 4,
                         'patient_level_rho': rho(patient_medians.Core50, patient_medians[name])})
    agg_rows.append({'comparison': 'MS3_like_vs_MS6_like', 'n_patients': 4,
                     'patient_level_rho': rho(patient_medians.MS3_like, patient_medians.MS6_like)})
    lopo = []
    for leave in patients:
        sub = patient_medians.loc[patient_medians.canonical_patient_id != leave]
        for name in PROGRAM3:
            lopo.append({'held_out_patient': leave, 'comparison': f'Core50_vs_{name}',
                         'n_patients': 3, 'rho': rho(sub.Core50, sub[name])})
        lopo.append({'held_out_patient': leave, 'comparison': 'MS3_like_vs_MS6_like',
                     'n_patients': 3, 'rho': rho(sub.MS3_like, sub.MS6_like)})
    spots.to_csv(TABLES / 'Step32EF_SpotScores.csv', index=False)
    pd.DataFrame(detection).to_csv(AUDIT / 'Step32_SelectedGeneDetection.csv', index=False)
    e_section.to_csv(TABLES / 'Step32E_Core50ProgramBySection.csv', index=False)
    e_patient.to_csv(TABLES / 'Step32E_Core50ProgramByPatient.csv', index=False)
    f_section.to_csv(TABLES / 'Step32F_MS3MS6BySection.csv', index=False)
    f_patient.to_csv(TABLES / 'Step32F_MS3MS6ByPatient.csv', index=False)
    section_medians.to_csv(TABLES / 'Step32EF_SectionMedianScores.csv', index=False)
    patient_medians.to_csv(TABLES / 'Step32EF_PatientMedianScores.csv', index=False)
    pd.DataFrame(agg_rows).to_csv(TABLES / 'Step32EF_PatientLevelAssociation.csv', index=False)
    pd.DataFrame(lopo).to_csv(TABLES / 'Step32EF_LOPO.csv', index=False)
    method = {'normalization': 'log1p(raw count / spot total counts * 10000)',
              'gene_standardization': 'global per-gene z score across 27224 spatial spots, ddof=0',
              'program_score': 'arithmetic mean of frozen gene z scores, zero-SD genes contribute zero',
              'Core50_mapping': str(STEP31 / 'Step31_Core50_FrozenValidationMapping.csv'),
              'program_mapping': str(STEP31 / 'Step31E2A_ValidationProgramMapping.csv'),
              'MS3_like': PROGRAM3, 'MS6_like': PROGRAM6,
              'Core50_in_reciprocal_axis': False,
              'association': 'Spearman rho; spot-level descriptive, patient-level n=4 independent units',
              'patient_aggregation': 'median spot score per section, then median section score per patient',
              'spatial_neighbor': 'six nearest non-self spots in full-resolution pixel coordinates'}
    (AUDIT / 'Step32EF_ScoringMethod.json').write_text(json.dumps(method, indent=2) + '\n')
    lock = {'stage': 'Step32E/F', 'status': 'PASS', 'n_spots': len(spots),
            'n_sections': len(samples), 'n_patients': len(patients),
            'n_core50_genes': len(programs['Core50']),
            'program_gene_counts': {k: len(v) for k, v in programs.items()},
            'n_finite_section_associations': int(e_section.rho_core50_program.notna().sum()),
            'n_finite_axis_associations': int(f_section.rho_within_spot_ms3_ms6.notna().sum()),
            'patient_level_associations': agg_rows}
    (AUDIT / 'Step32EF_FINAL_LOCK.json').write_text(json.dumps(lock, indent=2) + '\n')
    print(json.dumps(lock, indent=2))


if __name__ == '__main__':
    main()

if __name__ == '__main__':
    print('REPORT: Step32E primary Core50/proliferative concordance; see outputs listed in header.')
