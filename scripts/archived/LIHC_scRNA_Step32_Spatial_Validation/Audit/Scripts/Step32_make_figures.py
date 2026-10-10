#!/usr/bin/env python3
# Step: Step32 figures
# Conda environment: spatial32_wsl
# Python version: 3.10.21
# Frozen dependencies: /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_environment_freeze.txt (SHA256 6fe6d67eeb0584744c98b7242f6b56375d245b03808b3e19f78257181542c395)
# Relevant package versions: numpy==2.2.6; pandas==2.3.3; matplotlib==3.10.9
# Executed source: /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_make_figures.py
# Executed source SHA256: 8f0f7298be871cb725c912f5fad8a5b57dfdfa6b4d41fc76bea714423c1b5e42
# Inputs:
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Tables/Step32EF_SpotScores.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Tables/Step32E_Core50ProgramBySection.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Tables/Step32EF_TechnicalSensitivityBySection.csv
#   /mnt/c/Users/36888/.codex/skills/nature-figure/scripts/audit_panel_alignment.py
# Outputs:
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Figures/Step32*.pdf
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Figures/Step32*.svg
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Figures/Step32*.png
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Figures/Step32*.alignment.json
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_FigureProvenance.json
# Source body below is the executed version; only this header and final REPORT were added.
"""Compact Python figures for frozen Step32 primary and secondary spatial results."""
from hcc_release_paths import project_path as _hcc_path

import json
import sys
from pathlib import Path

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, '/mnt/c/Users/36888/.codex/skills/nature-figure/scripts')
from audit_panel_alignment import require_matplotlib_panel_alignment


ROOT = Path(_hcc_path('LIHC_scRNA_Step32_Spatial_Validation'))
FIG = ROOT / 'Figures'
TABLES = ROOT / 'Tables'
AUDIT = ROOT / 'Audit'
FIG.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans', 'sans-serif'],
    'pdf.fonttype': 42,
    'svg.fonttype': 'none',
    'font.size': 7,
    'axes.labelsize': 7,
    'axes.titlesize': 8,
    'xtick.labelsize': 6.5,
    'ytick.labelsize': 6.5,
    'legend.fontsize': 6.5,
    'axes.spines.right': False,
    'axes.spines.top': False,
    'axes.linewidth': .7,
    'savefig.facecolor': 'white',
})


def save(fig, stem):
    fig.canvas.draw()
    report = require_matplotlib_panel_alignment(
        fig, json_out=FIG / f'{stem}.alignment.json', strict=True)
    if report['verdict'] not in {'NOT APPLICABLE', 'PASS'}:
        raise ValueError(f'Figure alignment unresolved: {stem}: {report["verdict"]}')
    fig.savefig(FIG / f'{stem}.svg', bbox_inches='tight')
    fig.savefig(FIG / f'{stem}.pdf', bbox_inches='tight')
    fig.savefig(FIG / f'{stem}.png', dpi=600, bbox_inches='tight')
    plt.close(fig)


def spatial_map(data, metric, title, stem, maxabs):
    fig, ax = plt.subplots(figsize=(3.5, 3.4), layout='constrained')
    points = ax.scatter(data.spatial_x, data.spatial_y, c=data[metric],
                        cmap='coolwarm', vmin=-maxabs, vmax=maxabs,
                        s=8, linewidths=0, rasterized=True)
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.set_xlabel('Full-resolution x (pixels)')
    ax.set_ylabel('Full-resolution y (pixels)')
    ax.set_title(title, loc='left', fontweight='bold')
    bar = fig.colorbar(points, ax=ax, fraction=.045, pad=.03)
    bar.set_label('Mean gene z score')
    fig.suptitle('GSM8627602 · A009T1 · n = 1,762 spots', fontsize=7, y=.99)
    save(fig, stem)


def primary_consistency(data):
    patient_colors = {'A009': '#243b58', 'B013': '#317f8c',
                      'C004': '#b07935', 'C008': '#784b80'}
    names = ['G2M_Independent', 'Mitotic_Entry_Independent',
             'Spindle_Chromosome_Independent']
    labels = ['G2/M', 'Mitotic entry', 'Spindle / chromosome']
    fig, ax = plt.subplots(figsize=(6.5, 3.3), layout='constrained')
    offsets = {'A009': -.18, 'B013': -.06, 'C004': .06, 'C008': .18}
    ax.axvline(0, color='#777777', lw=.7, ls='--', zorder=0)
    for j, (name, label) in enumerate(zip(names, labels)):
        for patient, sub in data.loc[data.program == name].groupby('canonical_patient_id'):
            y = np.full(len(sub), j + offsets[patient])
            ax.scatter(sub.rho_core50_program, y,
                       s=22, color=patient_colors[patient], alpha=.78,
                       edgecolor='white', linewidth=.25, zorder=3,
                       label=patient if j == 0 else None)
        ax.text(.985, j, f'{(data.loc[data.program == name, "rho_core50_program"] > 0).sum()}/17 positive',
                transform=ax.get_yaxis_transform(), va='center', ha='right', fontsize=6.5)
    ax.set_yticks(range(3), labels)
    ax.set_ylim(2.55, -.55)
    ax.set_xlim(-.02, .98)
    ax.set_xlabel('Section-level Spearman ρ with Core50')
    ax.set_title('Core50 concordance across all 17 sections', loc='left', fontweight='bold')
    ax.legend(title='Canonical patient', ncol=4, loc='upper center',
              bbox_to_anchor=(.5, -.24), frameon=False, title_fontsize=6.5)
    save(fig, 'Step32E_Core50_SectionConsistency')


def reciprocal_boundary(data):
    patient_colors = {'A009': '#243b58', 'B013': '#317f8c',
                      'C004': '#b07935', 'C008': '#784b80'}
    data = data.sort_values(['canonical_patient_id', 'gsm']).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(6.5, 4.2), layout='constrained')
    ax.axvline(0, color='#777777', lw=.7, ls='--', zorder=0)
    for j, row in data.iterrows():
        color = patient_colors[row.canonical_patient_id]
        ax.plot([row.raw_rho, row.technical_partial_rho], [j, j],
                color=color, lw=.9, alpha=.65, zorder=1)
        ax.scatter(row.raw_rho, j, s=20, color=color, zorder=2)
        ax.scatter(row.technical_partial_rho, j, s=20, facecolor='white',
                   edgecolor=color, lw=1, zorder=3)
    ax.set_yticks(range(len(data)),
                  [f'{row.canonical_patient_id} · {row.gsm[-3:]}'
                   for _, row in data.iterrows()])
    ax.set_ylim(len(data) - .5, -.5)
    ax.set_xlim(-.52, .65)
    ax.set_xlabel('MS3-like versus MS6-like Spearman ρ')
    ax.set_title('Reciprocal spatial organization is not consistent',
                 loc='left', fontweight='bold')
    ax.text(.98, .02, 'Filled: raw   Open: technical-adjusted',
            transform=ax.transAxes, ha='right', va='bottom', fontsize=6.5)
    save(fig, 'Step32F_MS3_MS6_SectionBoundary')


def main():
    scores = pd.read_csv(TABLES / 'Step32EF_SpotScores.csv')
    selected = scores.loc[scores.gsm == 'GSM8627602']
    if len(selected) != 1762:
        raise ValueError('Prespecified representative section changed')
    maxabs = float(np.quantile(np.abs(scores[['Core50', 'MS3_like', 'MS6_like']].to_numpy()), .98))
    maps = [
        ('Core50', 'Frozen MS3 Core50', 'Step32E_Representative_Core50_Map'),
        ('MS3_like', 'Independent proliferative programs', 'Step32E_Representative_MS3_Map'),
        ('MS6_like', 'Hepatic / metabolic programs', 'Step32F_Representative_MS6_Map'),
    ]
    for metric, title, stem in maps:
        spatial_map(selected, metric, title, stem, maxabs)
    primary_consistency(pd.read_csv(TABLES / 'Step32E_Core50ProgramBySection.csv'))
    sensitivity = pd.read_csv(TABLES / 'Step32EF_TechnicalSensitivityBySection.csv')
    reciprocal_boundary(sensitivity.loc[sensitivity.comparison == 'MS3_like_vs_MS6_like'])
    note = {'backend': 'Python/matplotlib', 'representative_section': 'GSM8627602',
            'representative_selection': 'first T-coded section in frozen GSM order; not selected by score',
            'spatial_coordinates': 'full-resolution pixel coordinates; no histology manipulation',
            'map_color_clip': f'symmetric ±{maxabs:.4f} pooled 98th absolute percentile; scores unchanged',
            'source_spot_scores': str(TABLES / 'Step32EF_SpotScores.csv'),
            'source_associations': [str(TABLES / 'Step32E_Core50ProgramBySection.csv'),
                                    str(TABLES / 'Step32EF_TechnicalSensitivityBySection.csv')],
            'n_sections': 17, 'n_canonical_patients': 4,
            'interpretation': 'Core50 positive concordance; reciprocal axis heterogeneous and unsupported as uniform'}
    (AUDIT / 'Step32_FigureProvenance.json').write_text(json.dumps(note, indent=2) + '\n')
    print(json.dumps(note, indent=2))


if __name__ == '__main__':
    main()

if __name__ == '__main__':
    print('REPORT: Five English-text figures, PDF/SVG/600-dpi PNG; see outputs listed in header.')
