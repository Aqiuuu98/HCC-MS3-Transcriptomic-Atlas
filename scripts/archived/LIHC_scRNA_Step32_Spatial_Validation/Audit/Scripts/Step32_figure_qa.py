#!/usr/bin/env python3
# Step: Step32 figure QA
# Conda environment: spatial32_wsl
# Python version: 3.10.21
# Frozen dependencies: /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_environment_freeze.txt (SHA256 6fe6d67eeb0584744c98b7242f6b56375d245b03808b3e19f78257181542c395)
# Relevant package versions: matplotlib==3.10.9
# Executed source: /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_figure_qa.py
# Executed source SHA256: d7cdd4ef0091cd065822f9af886c52750927ecf1bdd4ea8acd1a7c515189db86
# Inputs:
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_make_figures.py
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Figures/Step32*.pdf
#   /mnt/c/Users/36888/.codex/skills/nature-figure/scripts/validate_figure.py
#   /mnt/c/Users/36888/.codex/skills/nature-figure/scripts/audit_pdf_text.py
#   /mnt/c/Users/36888/.codex/skills/nature-figure/scripts/audit_figure_collisions.py
# Outputs:
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_FigureQA.json
# Source body below is the executed version; only this header and final REPORT were added.
"""Run available source/PDF checks for Step32 figures; record dependency limits."""
from hcc_release_paths import project_path as _hcc_path

import importlib.util
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(_hcc_path('LIHC_scRNA_Step32_Spatial_Validation'))
SKILL = Path('/mnt/c/Users/36888/.codex/skills/nature-figure/scripts')
FIG = ROOT / 'Figures'
AUDIT = ROOT / 'Audit'


def run(args):
    result = subprocess.run([sys.executable, *map(str, args)],
                            capture_output=True, text=True)
    return {'exit_code': result.returncode,
            'output': (result.stdout + result.stderr).strip()}


def main():
    source = run([SKILL / 'validate_figure.py', AUDIT / 'Step32_make_figures.py'])
    rows = []
    for pdf in sorted(FIG.glob('Step32*.pdf')):
        text = run([SKILL / 'audit_pdf_text.py', pdf, '--min-pt', '5'])
        alignment = json.loads(pdf.with_suffix('.alignment.json').read_text())
        rows.append({'figure': pdf.name, 'pdf_text_audit_exit_code': text['exit_code'],
                     'pdf_text_audit': text['output'],
                     'panel_alignment_verdict': alignment.get('verdict'),
                     'collision_audit': 'unavailable: PyMuPDF not installed in spatial32_wsl'})
        print(pdf.name, text['exit_code'], alignment.get('verdict'))
    report = {'source_preflight_exit_code': source['exit_code'],
              'source_preflight': source['output'],
              'pymupdf_available': importlib.util.find_spec('fitz') is not None,
              'figure_checks': rows,
              'visual_inspection': 'Core50 map, Core50 consistency, and MS3-MS6 boundary PNGs inspected at rendered size; no visible clipping after legend correction',
              'limitation': 'Rendered PDF collision audit could not run without PyMuPDF; no package was installed.'}
    (AUDIT / 'Step32_FigureQA.json').write_text(json.dumps(report, indent=2) + '\n')
    print('source_preflight_exit_code', source['exit_code'])


if __name__ == '__main__':
    main()

if __name__ == '__main__':
    print('REPORT: Source and PDF text QA; PyMuPDF absence recorded; see outputs listed in header.')
