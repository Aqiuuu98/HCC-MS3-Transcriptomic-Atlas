#!/usr/bin/env python3
# Step: Step32 final lock
# Conda environment: spatial32_wsl
# Python version: 3.10.21
# Frozen dependencies: /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_environment_freeze.txt (SHA256 6fe6d67eeb0584744c98b7242f6b56375d245b03808b3e19f78257181542c395)
# Relevant package versions:
# Executed source: /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_finalize.py
# Executed source SHA256: f73493283961f6ecdc0663da9ce2bd49a86de0e2047269f4204dc515473c78f4
# Inputs:
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32C2_FINAL_LOCK.json
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32D_FINAL_LOCK.json
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32EF_FINAL_LOCK.json
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_FigureQA.json
#   /mnt/e/TCGA/Step32_Spatial_Validation.ipynb
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Data/**
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Tables/**
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Figures/**
# Outputs:
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_FinalAudit.json
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_OutputManifest.csv
# Source body below is the executed version; only this header and final REPORT were added.
"""Final file inventory and cross-stage lock audit; no source-data mutation."""
from hcc_release_paths import project_path as _hcc_path

import csv
import json
from pathlib import Path


ROOT = Path(_hcc_path('LIHC_scRNA_Step32_Spatial_Validation'))
AUDIT = ROOT / 'Audit'
NOTEBOOK = Path(_hcc_path('Step32_Spatial_Validation.ipynb'))


def main():
    c2 = json.loads((AUDIT / 'Step32C2_FINAL_LOCK.json').read_text())
    d = json.loads((AUDIT / 'Step32D_FINAL_LOCK.json').read_text())
    ef = json.loads((AUDIT / 'Step32EF_FINAL_LOCK.json').read_text())
    qa = json.loads((AUDIT / 'Step32_FigureQA.json').read_text())
    notebook = json.loads(NOTEBOOK.read_text(encoding='utf-8'))
    raw = list((ROOT / 'Data/Raw/GSE281759').glob('*'))
    std = list((ROOT / 'Data/Standardized/GSE281759').rglob('*'))
    std_files = [x for x in std if x.is_file()]
    h5ad = list((ROOT / 'Data/QC/GSE281759').glob('*_raw_qc.h5ad'))
    pdf = list((ROOT / 'Figures').glob('Step32*.pdf'))
    png = list((ROOT / 'Figures').glob('Step32*.png'))
    svg = list((ROOT / 'Figures').glob('Step32*.svg'))
    if not (c2['status'] == d['status'] == ef['status'] == 'PASS'):
        raise ValueError('One or more stage locks failed')
    if (len(raw), len(std_files), len(h5ad), len(pdf), len(png), len(svg)) != (153, 119, 17, 5, 5, 5):
        raise ValueError('Final output file counts differ from expected')
    if len(notebook['cells']) != 51:
        raise ValueError('Notebook appendix missing')
    if qa['source_preflight_exit_code'] != 0 or not all(
            x['pdf_text_audit_exit_code'] == 0 for x in qa['figure_checks']):
        raise ValueError('Available figure QA did not pass')
    report = {
        'status': 'COMPLETE_WITH_SECONDARY_NEGATIVE_RESULT',
        'C2': c2['status'], 'D': d['status'], 'E_F': ef['status'],
        'raw_files': len(raw), 'standardized_files': len(std_files),
        'raw_qc_h5ad_files': len(h5ad), 'figures_per_format': len(pdf),
        'notebook_cells': len(notebook['cells']),
        'canonical_patients': 4, 'sections': 17, 'spots': 27224,
        'primary_core50_section_program_positive': '51/51 raw; 51/51 technically adjusted',
        'secondary_reciprocal_axis': 'unsupported as uniform spatial relationship',
        'figure_collision_audit': 'unavailable: PyMuPDF absent; no installation',
    }
    (AUDIT / 'Step32_FinalAudit.json').write_text(json.dumps(report, indent=2) + '\n')
    rows = []
    for path in sorted(ROOT.rglob('*')):
        if not path.is_file() or path.name == 'Step32_OutputManifest.csv':
            continue
        relative = path.relative_to(ROOT).as_posix()
        if relative.startswith('Data/Raw/'):
            role = 'extracted_raw_member'
        elif relative.startswith('Data/Standardized/'):
            role = 'standardized_input'
        elif relative.startswith('Data/QC/'):
            role = 'frozen_raw_qc_h5ad'
        elif relative.startswith('Figures/'):
            role = 'figure_or_figure_qa'
        elif relative.startswith('Tables/'):
            role = 'analysis_table'
        elif relative.startswith('Audit/Step32A') or relative.startswith('Audit/Step32C1'):
            role = 'preexisting_frozen_audit'
        elif relative.endswith('.py'):
            role = 'reproducible_code'
        else:
            role = 'audit_or_report'
        rows.append({'path': str(path), 'relative_path': relative,
                     'size_bytes': path.stat().st_size, 'role': role})
    rows.append({'path': str(NOTEBOOK), 'relative_path': '../Step32_Spatial_Validation.ipynb',
                 'size_bytes': NOTEBOOK.stat().st_size, 'role': 'notebook'})
    rows.append({'path': str(AUDIT / 'Step32_OutputManifest.csv'),
                 'relative_path': 'Audit/Step32_OutputManifest.csv',
                 'size_bytes': '', 'role': 'self_manifest'})
    with (AUDIT / 'Step32_OutputManifest.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, ['path', 'relative_path', 'size_bytes', 'role'])
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(report, indent=2))
    print('manifest_rows', len(rows))


if __name__ == '__main__':
    main()

if __name__ == '__main__':
    print('REPORT: Original 51-cell final snapshot and output inventory; see outputs listed in header.')
