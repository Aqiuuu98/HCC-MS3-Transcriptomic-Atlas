#!/usr/bin/env python3
"""Freeze article-centered Step32 conclusions from existing audit tables only."""
from hcc_release_paths import project_path as _hcc_path

# Step: Step32 article-centered finalization
# Inputs: existing Step32E/F section tables, Step32EF_FINAL_LOCK.json,
#         Step32_Reproducibility_Manifest.json, Step32_CODE_INDEX.md,
#         Step32_MASTER_REPORT.md, Step32_Spatial_Validation.ipynb
# Outputs: Step32_FINAL_LOCK.json, Step32_FINAL_ArticleCenteredSummary.csv,
#          Step32_MASTER_REPORT.md, Step32_CODE_INDEX.md,
#          Step32_Reproducibility_Manifest.json, Notebook finalization cells
# Conda environment: spatial32_wsl
# Python version: 3.10.21
# Frozen dependencies: Audit/Step32_environment_freeze.txt
# Standard-library only; no spatial expression analysis is run here.

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(_hcc_path('LIHC_scRNA_Step32_Spatial_Validation'))
AUDIT = ROOT / 'Audit'
TABLES = ROOT / 'Tables'
NOTEBOOK = Path(_hcc_path('Step32_Spatial_Validation.ipynb'))
SCRIPT = AUDIT / 'Scripts/Step32_article_finalize.py'
REPORT = ROOT / 'Step32_MASTER_REPORT.md'
SUMMARY = TABLES / 'Step32_FINAL_ArticleCenteredSummary.csv'
LOCK = AUDIT / 'Step32_FINAL_LOCK.json'
INDEX = AUDIT / 'Step32_CODE_INDEX.md'
MANIFEST = AUDIT / 'Step32_Reproducibility_Manifest.json'
ENV = AUDIT / 'Step32_environment_freeze.txt'
E_TABLE = TABLES / 'Step32E_Core50ProgramBySection.csv'
TECH_TABLE = TABLES / 'Step32EF_TechnicalSensitivityBySection.csv'
F_TABLE = TABLES / 'Step32F_MS3MS6BySection.csv'
EF_LOCK = AUDIT / 'Step32EF_FINAL_LOCK.json'
PROGRAMS = ('G2M_Independent', 'Mitotic_Entry_Independent',
            'Spindle_Chromosome_Independent')
TAG = 'step32_article_finalization'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle))


def write_json(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def append_notebook_record():
    nb = json.loads(NOTEBOOK.read_text(encoding='utf-8'))
    if any(TAG in cell.get('metadata', {}).get('tags', []) for cell in nb['cells']):
        assert len(nb['cells']) == 56, 'Unexpected finalization Notebook cell count'
        return len(nb['cells'])
    assert len(nb['cells']) == 54, 'Expected the 54-cell Step32 Notebook before finalization'
    markdown = (
        '## Step32 final article-centered lock\n\n'
        '**Primary.** Frozen MS3 Core50 has positive spatial concordance with G2M, '
        'mitotic-entry, and spindle/chromosome-segregation programs in all 17 Visium '
        'sections from four canonical patients. All 51 section-program associations '
        'remain positive after technical-complexity adjustment.\n\n'
        '**Secondary negative finding.** The MS3-like versus MS6-like transcriptional '
        'relationship does not reproduce as a consistent spatially exclusive pattern: '
        'unadjusted same-spot associations are positive in 14/17 sections; adjusted '
        'directions are 8 positive and 9 negative. Do not call this spatial axis '
        'reciprocal or exclusive.\n\n'
        '**Boundaries.** Sections are nested within four patients; spots are not '
        'independent biological replicates. Spatial concordance does not establish '
        'causality. The negative spatial result does not invalidate the prior '
        'transcriptional reciprocal relationship. No new spatial analysis was run '
        'for this lock.\n'
    )
    code = (
        "# Recreate only the article-centered lock from existing Step32 audit outputs.\n"
        "import runpy\n"
        "runpy.run_path('/mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/"
        "Audit/Scripts/Step32_article_finalize.py', run_name='__main__')\n"
    )
    nb['cells'].extend([
        {'cell_type': 'markdown', 'id': 'step32-final-article-markdown',
         'metadata': {'tags': [TAG]}, 'source': markdown.splitlines(keepends=True)},
        {'cell_type': 'code', 'id': 'step32-final-article-code',
         'metadata': {'tags': [TAG]}, 'execution_count': None, 'outputs': [],
         'source': code.splitlines(keepends=True)},
    ])
    NOTEBOOK.write_text(json.dumps(nb, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    return len(nb['cells'])


def main():
    e, technical, f = map(read_csv, (E_TABLE, TECH_TABLE, F_TABLE))
    ef = json.loads(EF_LOCK.read_text(encoding='utf-8'))
    section_ids = {r['gsm'] for r in f}
    patient_ids = {r['canonical_patient_id'] for r in f}
    assert len(f) == len(section_ids) == 17
    assert patient_ids == {'A009', 'B013', 'C004', 'C008'}
    assert ef['n_sections'] == 17 and ef['n_patients'] == 4
    assert len(e) == 51 and {r['gsm'] for r in e} == section_ids
    assert all({r['program'] for r in e if r['gsm'] == gsm} == set(PROGRAMS)
               for gsm in section_ids)
    assert all(float(r['rho_core50_program']) > 0 for r in e)
    primary_adjusted = [r for r in technical if r['comparison'].startswith('Core50_vs_')]
    axis_adjusted = [r for r in technical if r['comparison'] == 'MS3_like_vs_MS6_like']
    assert len(primary_adjusted) == 51 and {r['gsm'] for r in primary_adjusted} == section_ids
    assert all(float(r['technical_partial_rho']) > 0 for r in primary_adjusted)
    assert len(axis_adjusted) == 17 and {r['gsm'] for r in axis_adjusted} == section_ids
    assert sum(float(r['rho_within_spot_ms3_ms6']) > 0 for r in f) == 14
    assert sum(float(r['technical_partial_rho']) > 0 for r in axis_adjusted) == 8
    assert sum(float(r['technical_partial_rho']) < 0 for r in axis_adjusted) == 9
    assert ENV.exists() and all(s in ENV.read_text(encoding='utf-8')
                                for s in ('python --version', 'pip freeze', 'conda list'))

    notebook_cells = append_notebook_record()
    rows = [
        {'finding_id': 'PRIMARY', 'status': 'SUPPORTED',
         'article_centered_conclusion': 'Frozen MS3 Core50 shows reproducible spatial concordance with G2M, mitotic-entry, and spindle/chromosome-segregation programs across all 17 Visium sections from 4 canonical patients.',
         'observed_result': '51/51 positive unadjusted section-program associations; 51/51 positive after technical-complexity adjustment',
         'unit_and_boundary': '17 sections nested within 4 patients; spots are not independent biological replicates; concordance does not establish causality',
         'evidence_files': f'{E_TABLE}; {TECH_TABLE}'},
        {'finding_id': 'SECONDARY_NEGATIVE', 'status': 'UNSUPPORTED_AS_UNIFORM_SPATIAL_PATTERN',
         'article_centered_conclusion': 'The MS3-like versus MS6-like reciprocal transcriptional relationship does not reproduce as a consistent spatially exclusive pattern.',
         'observed_result': 'Unadjusted same-spot associations positive in 14/17 sections; adjusted directions 8 positive and 9 negative',
         'unit_and_boundary': 'Do not describe the MS3-MS6 axis as spatially reciprocal or spatially exclusive; this does not invalidate the prior transcriptional reciprocal relationship',
         'evidence_files': f'{F_TABLE}; {TECH_TABLE}'},
    ]
    with SUMMARY.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    start, end = '<!-- STEP32_ARTICLE_FINAL_START -->', '<!-- STEP32_ARTICLE_FINAL_END -->'
    report = REPORT.read_text(encoding='utf-8')
    if start in report:
        report = report.split(start)[0].rstrip() + '\n'
    final_block = f'''\n{start}
## Final article-centered lock

**PRIMARY — supported.** Frozen MS3 Core50 shows reproducible spatial concordance with G2M, mitotic-entry, and spindle/chromosome-segregation programs across all 17 Visium sections from 4 canonical patients. All 51 section-program associations are positive in the existing unadjusted results and remain positive after technical-complexity adjustment.

**SECONDARY NEGATIVE FINDING.** The MS3-like versus MS6-like reciprocal transcriptional relationship does not reproduce as a consistent spatially exclusive pattern. Unadjusted same-spot associations are positive in 14/17 sections. After technical-complexity adjustment, 8 directions are positive and 9 are negative. The MS3-MS6 axis is **not** described as spatially reciprocal or spatially exclusive.

**Interpretation boundaries.** The 17 sections are nested within 4 canonical patients (A009, B013, C004, C008). Spots are not independent biological replicates. Spatial concordance does not establish causality. The lack of a consistent spatial reciprocal pattern does not invalidate the transcriptional reciprocal relationship from previous steps. No new signatures, clusters, pathways, or communication analyses were added.

Final lock: `Audit/Step32_FINAL_LOCK.json`. Article summary: `Tables/Step32_FINAL_ArticleCenteredSummary.csv`. Reproducibility files: `Audit/Step32_CODE_INDEX.md`, `Audit/Step32_Reproducibility_Manifest.json`, `Audit/Step32_environment_freeze.txt`. Notebook code and standalone executed-code copies are retained. Step32G was not executed. The historical `Audit/Step32_FinalAudit.json` and `Audit/Step32_OutputManifest.csv` remain prior frozen snapshots.

**STEP32: FROZEN**  
**NEXT: Step33 — virtual gene knockout**
{end}
'''
    REPORT.write_text(report.rstrip() + '\n' + final_block, encoding='utf-8')

    index = INDEX.read_text(encoding='utf-8')
    index_marker = '<!-- STEP32_ARTICLE_INDEX_START -->'
    if index_marker in index:
        index = index.split(index_marker)[0].rstrip() + '\n'
    index += (f'\n{index_marker}\n## Article-centered finalization\n\n'
              '| Stage | Notebook cells | Standalone script | Inputs | Outputs | Frozen status |\n'
              '|---|---:|---|---|---|---|\n'
              '| Step32 final article lock | 55-56 | `Audit/Scripts/Step32_article_finalize.py` | Existing Step32E/F section association tables; Step32EF_FINAL_LOCK.json; Notebook; existing report/index/manifest | `Audit/Step32_FINAL_LOCK.json`; `Tables/Step32_FINAL_ArticleCenteredSummary.csv`; `Step32_MASTER_REPORT.md`; `Audit/Step32_CODE_INDEX.md`; `Audit/Step32_Reproducibility_Manifest.json` | FROZEN; audit-only, no new spatial analysis |\n\n'
              'The Notebook finalization code invokes this standalone script. Its inputs, outputs, and fixed interpretation counts are explicit in the script and final lock.\n')
    INDEX.write_text(index, encoding='utf-8')

    manifest = json.loads(MANIFEST.read_text(encoding='utf-8'))
    manifest['status'] = 'STEP32_FROZEN'
    manifest['notebook_cells_after_addendum'] = notebook_cells
    manifest['notebook_code_retained'] = True
    manifest['article_centered_finalization'] = {
        'script': str(SCRIPT), 'script_sha256': sha(SCRIPT),
        'notebook_cells_1_based': '55-56',
        'outputs': [str(LOCK), str(SUMMARY), str(REPORT), str(INDEX), str(MANIFEST)],
        'primary_positive_raw': '51/51', 'primary_positive_adjusted': '51/51',
        'secondary_raw_positive': '14/17',
        'secondary_adjusted_positive': '8/17', 'secondary_adjusted_negative': '9/17',
        'spatial_reciprocity_supported': False,
        'no_new_spatial_analysis': True,
    }
    manifest['scripts'] = [s for s in manifest['scripts']
                           if s['step'] != 'Step32 article-centered finalization']
    manifest['scripts'].append({
        'step': 'Step32 article-centered finalization',
        'notebook_cells_1_based': '55-56',
        'standalone_script': str(SCRIPT), 'standalone_sha256': sha(SCRIPT),
        'executed_source_script': str(SCRIPT), 'executed_source_sha256': sha(SCRIPT),
        'executed_source_body_retained': True,
        'inputs': [str(E_TABLE), str(TECH_TABLE), str(F_TABLE), str(EF_LOCK),
                   str(NOTEBOOK), str(REPORT), str(INDEX), str(MANIFEST)],
        'outputs': [str(LOCK), str(SUMMARY), str(REPORT), str(INDEX), str(MANIFEST)],
        'frozen_status': 'FROZEN; audit-only', 'dependencies': ['python standard library'],
    })
    manifest['verification']['standalone_script_count'] = len(manifest['scripts'])
    manifest['verification']['article_finalization_evidence_checked'] = True
    write_json(MANIFEST, manifest)

    lock = {
        'status': 'STEP32_FROZEN',
        'scope': 'Article-centered conclusion lock from existing Step32 outputs; no new spatial analysis',
        'dataset': 'GSE281759, 10x Genomics Visium',
        'n_sections': 17, 'n_canonical_patients': 4,
        'canonical_patient_ids': sorted(patient_ids),
        'primary': {
            'conclusion': rows[0]['article_centered_conclusion'],
            'programs': list(PROGRAMS),
            'unadjusted_positive_section_program_associations': '51/51',
            'technically_adjusted_positive_section_program_associations': '51/51',
        },
        'secondary_negative_finding': {
            'conclusion': rows[1]['article_centered_conclusion'],
            'unadjusted_same_spot_positive_sections': '14/17',
            'technically_adjusted_positive_sections': '8/17',
            'technically_adjusted_negative_sections': '9/17',
            'spatial_reciprocity_or_exclusivity_claim_allowed': False,
        },
        'boundaries': [
            '17 sections are nested within 4 patients.',
            'Spots are not independent biological replicates.',
            'Spatial concordance does not establish causality.',
            'The absent consistent spatial reciprocity does not invalidate the prior transcriptional reciprocal relationship.',
            'No new signatures, clusters, pathways, or communication analyses.',
        ],
        'evidence_sha256': {str(p): sha(p) for p in (E_TABLE, TECH_TABLE, F_TABLE, EF_LOCK)},
        'deliverables': {str(p): sha(p) for p in (SUMMARY, REPORT, INDEX, MANIFEST)},
        'notebook': str(NOTEBOOK), 'notebook_cells': notebook_cells,
        'standalone_scripts_retained': len(manifest['scripts']),
        'environment_freeze': str(ENV), 'environment_freeze_sha256': sha(ENV),
        'next_step': 'Step33 — virtual gene knockout',
    }
    write_json(LOCK, lock)
    print('REPORT: STEP32: FROZEN; 17 sections / 4 patients; primary 51/51 adjusted positive; secondary 14/17 raw positive, 8/9 adjusted direction split.')
    print('REPORT: Notebook code retained: YES')
    print('REPORT: Standalone scripts retained: YES')
    print('REPORT: Environment frozen: YES')
    print('REPORT: Reproducibility manifest created: YES')


if __name__ == '__main__':
    main()
