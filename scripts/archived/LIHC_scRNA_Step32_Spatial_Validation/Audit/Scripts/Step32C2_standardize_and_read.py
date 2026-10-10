#!/usr/bin/env python3
# Step: Step32C2
# Conda environment: spatial32_wsl
# Python version: 3.10.21
# Frozen dependencies: /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32_environment_freeze.txt (SHA256 6fe6d67eeb0584744c98b7242f6b56375d245b03808b3e19f78257181542c395)
# Relevant package versions: numpy==2.2.6; pandas==2.3.3; scanpy==1.11.5; scipy==1.15.3; pillow==12.3.0
# Executed source: /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32C2_run.py
# Executed source SHA256: 004f0abcbf3d93b8114828321ab5e05ca83ea68c86516d0a944ebdab7b8bffee
# Inputs:
#   /mnt/e/TCGA/LIHC/_Spatial/_Step32/_Transcriptomics/Data/Primary/GSE281759/GSE281759_RAW.tar
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32C1_GSE281759_SampleStructureAudit.csv
# Outputs:
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Data/Raw/GSE281759/*
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Data/Standardized/GSE281759/GSM*/**
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32C2_StandardizedFileMapping.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32C2_ExpressionReadAudit.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32C2_SpatialReadAudit.csv
#   /mnt/e/TCGA/LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32C2_FINAL_LOCK.json
# Source body below is the executed version; only this header and final REPORT were added.
"""Step32C2: safely stage frozen Visium tar and audit raw section readability."""
from hcc_release_paths import project_path as _hcc_path

import csv
import gc
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import tarfile
from collections import defaultdict
from pathlib import Path, PurePosixPath

import numpy as np
import pandas as pd
import scanpy as sc
from PIL import Image


ROOT = Path(_hcc_path('LIHC_scRNA_Step32_Spatial_Validation'))
AUDIT = ROOT / 'Audit'
RAW = ROOT / 'Data/Raw/GSE281759'
STD = ROOT / 'Data/Standardized/GSE281759'
TAR = Path(_hcc_path('LIHC/_Spatial/_Step32/_Transcriptomics/Data/Primary/GSE281759/GSE281759_RAW.tar'))
SAMPLES = AUDIT / 'Step32C1_GSE281759_SampleStructureAudit.csv'
EXPECTED_SHA256 = '75df5e7e81d0ceb9e1259e614d32caf7966c30107b3b92796d326aee7e59789b'
FILES = {
    'matrix.mtx.gz': '_matrix.mtx.gz',
    'barcodes.tsv.gz': '_barcodes.tsv.gz',
    'features.tsv.gz': '_features.tsv.gz',
    'spatial/tissue_positions_list.csv.gz': '_tissue_positions_list.csv.gz',
    'spatial/scalefactors_json.json': '_scalefactors_json.json.gz',
    'spatial/tissue_hires_image.png': '_tissue_hires_image.png.gz',
    'spatial/tissue_lowres_image.png': '_tissue_lowres_image.png.gz',
}


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_csv(path, rows, fields):
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fields)
        writer.writeheader()
        writer.writerows(rows)


def copy_if_needed(source, target, expected_size, decode_gzip=False):
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        if target.is_file() and target.stat().st_size == expected_size:
            return 'existing_size_match'
        raise ValueError(f'Existing output size mismatch; no overwrite: {target}')
    temp = target.with_name(target.name + f'.step32c2-{os.getpid()}.tmp')
    try:
        if decode_gzip:
            with gzip.open(source, 'rb') as src, temp.open('xb') as dst:
                shutil.copyfileobj(src, dst, length=4 * 1024 * 1024)
        else:
            with source.open('rb') as src, temp.open('xb') as dst:
                shutil.copyfileobj(src, dst, length=4 * 1024 * 1024)
        if temp.stat().st_size != expected_size:
            raise ValueError(f'Written size mismatch: {target}')
        if target.exists():
            raise ValueError(f'Output appeared during write; no overwrite: {target}')
        os.replace(temp, target)
        return 'created'
    finally:
        if temp.exists():
            temp.unlink()


def gzip_size(path):
    with path.open('rb') as stream:
        stream.seek(-4, os.SEEK_END)
        return int.from_bytes(stream.read(4), 'little')


def inspect_positions(path, expr_barcodes):
    # Explicitly check the real CSV form before assigning legacy columns.
    with gzip.open(path, 'rt', newline='') as stream:
        first = next(csv.reader(stream))
    if len(first) != 6 or first[0].lower() in {'barcode', 'barcodes'}:
        return {'format': f'non_legacy_{len(first)}_column', 'readable': False,
                'error': f'Unexpected first row: {first[:6]}'}
    names = ['barcode', 'in_tissue', 'array_row', 'array_col',
             'pxl_row_in_fullres', 'pxl_col_in_fullres']
    pos = pd.read_csv(path, header=None, names=names, compression='gzip',
                      dtype={'barcode': str})
    if len(pos.columns) != 6 or pos['barcode'].duplicated().any():
        return {'format': 'legacy_6_column', 'readable': False,
                'error': 'Malformed or duplicated position barcodes'}
    try:
        in_tissue = pd.to_numeric(pos['in_tissue'], errors='raise')
        if not in_tissue.isin([0, 1]).all():
            raise ValueError('in_tissue values outside {0,1}')
        for col in names[2:]:
            pd.to_numeric(pos[col], errors='raise')
    except Exception as exc:
        return {'format': 'legacy_6_column', 'readable': False, 'error': str(exc)}
    overlap = len(expr_barcodes.intersection(set(pos['barcode'])))
    return {'format': 'legacy_6_column', 'readable': True, 'error': '',
            'n_position_barcodes': len(pos), 'n_overlap': overlap,
            'overlap_fraction': overlap / len(expr_barcodes) if expr_barcodes else 0,
            'n_in_tissue_positions': int(in_tissue.sum()),
            'n_in_tissue_expression': int(pos.loc[
                pos['barcode'].isin(expr_barcodes), 'in_tissue'].sum()),
            'n_unmatched_expression': len(expr_barcodes) - overlap}


def main():
    AUDIT.mkdir(parents=True, exist_ok=True)
    RAW.mkdir(parents=True, exist_ok=True)
    STD.mkdir(parents=True, exist_ok=True)
    sample_rows = list(csv.DictReader(SAMPLES.open(newline='', encoding='utf-8')))
    if len(sample_rows) != 17 or len({r['gsm'] for r in sample_rows}) != 17:
        raise ValueError('Frozen C1 sample table does not contain 17 unique GSMs')
    if {r['canonical_patient_id'] for r in sample_rows} != {'A009', 'B013', 'C004', 'C008'}:
        raise ValueError('Frozen C1 patient hierarchy mismatch')
    if TAR.stat().st_size != 494551040 or sha256(TAR) != EXPECTED_SHA256:
        raise ValueError('Canonical tar differs from frozen provenance')

    tar_rows = {}
    with tarfile.open(TAR, 'r:') as archive:
        for member in archive.getmembers():
            name = PurePosixPath(member.name)
            if not member.isfile() or name.is_absolute() or len(name.parts) != 1 or '..' in name.parts:
                raise ValueError(f'Unsafe or unexpected tar member: {member.name}')
            if member.name in tar_rows:
                raise ValueError(f'Duplicate tar member: {member.name}')
            target = RAW / member.name
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.exists():
                if not target.is_file() or target.stat().st_size != member.size:
                    raise ValueError(f'Existing Raw file size mismatch: {target}')
                action = 'existing_size_match'
            else:
                temp = target.with_name(target.name + f'.step32c2-{os.getpid()}.tmp')
                try:
                    with archive.extractfile(member) as source, temp.open('xb') as sink:
                        shutil.copyfileobj(source, sink, length=4 * 1024 * 1024)
                    if temp.stat().st_size != member.size:
                        raise ValueError(f'Tar member size mismatch: {member.name}')
                    os.replace(temp, target)
                    action = 'created'
                finally:
                    if temp.exists():
                        temp.unlink()
            tar_rows[member.name] = {'size': member.size, 'raw_path': str(target),
                                     'raw_action': action}
    if len(tar_rows) != 153:
        raise ValueError(f'Expected 153 regular tar members, found {len(tar_rows)}')

    mappings = []
    for sample in sample_rows:
        gsm, label = sample['gsm'], sample['raw_sample_label']
        for relative, suffix in FILES.items():
            member = f'{gsm}_{label}{suffix}'
            if member not in tar_rows:
                raise ValueError(f'Missing tar member: {member}')
            source = RAW / member
            target = STD / gsm / relative
            decode = not relative.endswith('.gz')
            expected_size = gzip_size(source) if decode else source.stat().st_size
            action = copy_if_needed(source, target, expected_size, decode)
            mappings.append({'gsm': gsm, 'raw_sample_label': label,
                             'raw_patient_id': sample['raw_patient_id'],
                             'canonical_patient_id': sample['canonical_patient_id'],
                             'archive_member': member, 'raw_path': str(source),
                             'standardized_path': str(target),
                             'raw_size_bytes': source.stat().st_size,
                             'standardized_size_bytes': target.stat().st_size,
                             'action': action})

    expression = []
    spatial = []
    for sample in sample_rows:
        gsm = sample['gsm']
        base = STD / gsm
        er = {'gsm': gsm, 'raw_sample_label': sample['raw_sample_label'],
              'raw_patient_id': sample['raw_patient_id'],
              'canonical_patient_id': sample['canonical_patient_id'],
              'section_code': sample['section_code'], 'readable': False,
              'error': ''}
        sr = dict(er)
        adata = None
        try:
            adata = sc.read_10x_mtx(base, var_names='gene_ids', make_unique=False,
                                    gex_only=True)
            x = adata.X
            values = x.data if hasattr(x, 'data') else np.asarray(x).ravel()
            ids = pd.Index(adata.var_names.astype(str))
            symbols = adata.var['gene_symbols'] if 'gene_symbols' in adata.var else pd.Series([], dtype=str)
            barcodes = set(adata.obs_names.astype(str))
            er.update({'readable': True, 'n_spots': int(adata.n_obs),
                       'n_genes': int(adata.n_vars), 'nnz': int(x.nnz),
                       'total_counts': float(x.sum()),
                       'nonnegative': bool(np.isfinite(values).all() and (values >= 0).all()),
                       'integer_like_raw_counts': bool(np.allclose(values, np.rint(values), rtol=0, atol=1e-8)),
                       'unique_ensembl_gene_ids': bool(ids.is_unique and ids.str.match(r'^ENSG\d+(?:\.\d+)?$').all()),
                       'duplicate_gene_ids': int(ids.duplicated().sum()),
                       'n_gene_symbols_present': int(symbols.notna().sum() if len(symbols) else 0),
                       'gene_symbol_coverage': float((symbols.fillna('').astype(str).str.strip() != '').mean()) if len(symbols) else 0})
            sr['n_expression_barcodes'] = len(barcodes)
            sr.update(inspect_positions(base / 'spatial/tissue_positions_list.csv.gz', barcodes))
            with (base / 'spatial/scalefactors_json.json').open(encoding='utf-8') as stream:
                scales = json.load(stream)
            sr['scalefactors_readable'] = isinstance(scales, dict) and len(scales) > 0
            for kind in ('hires', 'lowres'):
                path = base / f'spatial/tissue_{kind}_image.png'
                with Image.open(path) as image:
                    image.verify()
                sr[f'{kind}_image_readable'] = True
            sr['barcode_alignment_ok'] = bool(sr.get('readable') and sr.get('n_unmatched_expression') == 0)
            sr['readable'] = bool(sr.get('readable') and sr['scalefactors_readable']
                                  and sr['hires_image_readable'] and sr['lowres_image_readable'])
        except Exception as exc:
            if not er['readable']:
                er['error'] = repr(exc)
            sr['error'] = repr(exc)
            sr['readable'] = False
        expression.append(er)
        spatial.append(sr)
        del adata
        gc.collect()
        print(f'{gsm}: expression={er["readable"]} spatial={sr["readable"]} '
              f'overlap={sr.get("overlap_fraction", "NA")}', flush=True)

    mapping_fields = ['gsm', 'raw_sample_label', 'raw_patient_id', 'canonical_patient_id',
                      'archive_member', 'raw_path', 'standardized_path',
                      'raw_size_bytes', 'standardized_size_bytes', 'action']
    expr_fields = ['gsm', 'raw_sample_label', 'raw_patient_id', 'canonical_patient_id',
                   'section_code', 'readable', 'n_spots', 'n_genes', 'nnz',
                   'total_counts', 'nonnegative', 'integer_like_raw_counts',
                   'unique_ensembl_gene_ids', 'duplicate_gene_ids',
                   'n_gene_symbols_present', 'gene_symbol_coverage', 'error']
    spatial_fields = ['gsm', 'raw_sample_label', 'raw_patient_id', 'canonical_patient_id',
                      'section_code', 'readable', 'format', 'n_expression_barcodes',
                      'n_position_barcodes', 'n_overlap', 'overlap_fraction',
                      'n_unmatched_expression', 'n_in_tissue_positions',
                      'n_in_tissue_expression', 'barcode_alignment_ok',
                      'scalefactors_readable', 'hires_image_readable',
                      'lowres_image_readable', 'error']
    write_csv(AUDIT / 'Step32C2_StandardizedFileMapping.csv', mappings, mapping_fields)
    write_csv(AUDIT / 'Step32C2_ExpressionReadAudit.csv', expression, expr_fields)
    write_csv(AUDIT / 'Step32C2_SpatialReadAudit.csv', spatial, spatial_fields)
    passed = (len(expression) == 17 and len(spatial) == 17
              and all(r.get('readable') and r.get('nonnegative')
                      and r.get('integer_like_raw_counts') for r in expression)
              and all(r.get('readable') and r.get('barcode_alignment_ok') for r in spatial))
    lock = {'stage': 'Step32C2', 'status': 'PASS' if passed else 'FAIL',
            'canonical_tar': str(TAR), 'canonical_tar_sha256': EXPECTED_SHA256,
            'n_archive_files': len(tar_rows), 'n_sections': len(sample_rows),
            'n_expression_readable': sum(bool(r.get('readable')) for r in expression),
            'n_nonnegative': sum(bool(r.get('nonnegative')) for r in expression),
            'n_integer_like': sum(bool(r.get('integer_like_raw_counts')) for r in expression),
            'n_spatial_readable': sum(bool(r.get('readable')) for r in spatial),
            'n_barcode_alignment_ok': sum(bool(r.get('barcode_alignment_ok')) for r in spatial),
            'n_canonical_patients': len({r['canonical_patient_id'] for r in sample_rows}),
            'anomalies': [r['gsm'] for r in expression if not r.get('readable') or not r.get('nonnegative') or not r.get('integer_like_raw_counts')]
                         + [r['gsm'] for r in spatial if not r.get('readable') or not r.get('barcode_alignment_ok')]}
    (AUDIT / 'Step32C2_FINAL_LOCK.json').write_text(json.dumps(lock, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(lock, indent=2), flush=True)
    if not passed:
        raise SystemExit(2)


if __name__ == '__main__':
    main()

if __name__ == '__main__':
    print('REPORT: Step32C2 read and barcode-alignment lock; see outputs listed in header.')
