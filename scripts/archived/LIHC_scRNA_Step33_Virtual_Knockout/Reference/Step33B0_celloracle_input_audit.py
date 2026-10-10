#!/usr/bin/env python3
"""Step33B0: audit raw-count provenance, frozen topology, endpoints, and pilot readiness.

Reads only frozen upstream artifacts. Writes Step33B0 audit artifacts under Step33/Audit.
No CellOracle import, installation, coefficient fitting, or perturbation is performed.
"""
from hcc_release_paths import project_path as _hcc_path

# Step: Step33B0 — CellOracle input/topology compatibility audit
# Inputs: exact paths in INPUTS below (frozen Step23/25/26/27/29/31/32/33A)
# Outputs: seven Step33B0 audit files in Audit/, plus Notebook report via companion builder
# Conda environment: spatial32_wsl
# Python version: 3.10.21 (runtime also recorded in EnvironmentPlan)
# Frozen dependencies: existing spatial32_wsl h5py, scipy, numpy; no new packages

import csv
import hashlib
import importlib.metadata as metadata
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse


TCGA = Path(_hcc_path(''))
ROOT = TCGA / 'LIHC_scRNA_Step33_Virtual_Knockout'
AUDIT = ROOT / 'Audit'
STEP23_RAW = TCGA / 'LIHC_scRNA_Step23_Malignant_State_Analysis/Data/GSE149614_STEP23D3_LINEAGE_CLEAN_RAW_COUNTS.h5ad'
STEP23A_RAW = TCGA / 'LIHC_scRNA_Step23_Malignant_State_Analysis/Data/GSE149614_STEP23A_MALIGNANT_RAW_COUNTS.h5ad'
STEP26_LOG = TCGA / 'LIHC_scRNA_Step26_MS3_Biology/Data/GSE149614_STEP26_FINAL_FROZEN_MS3_BIOLOGY.h5ad'
STEP29_RAW = TCGA / 'LIHC_scRNA_Step29_GRN/Data/Step29D_GRN_INPUT_FILTERED_RAW_COUNTS.h5ad'
STEP29_META = TCGA / 'LIHC_scRNA_Step29_GRN/Tables/Step29D_Frozen_Cell_Metadata.csv'
STEP29_REG = TCGA / 'LIHC_scRNA_Step29_GRN/Data/Step29F_Regulon_Targets.json'
STEP29_REG_SUM = TCGA / 'LIHC_scRNA_Step29_GRN/Tables/Step29F_Regulon_Summary.csv'
STEP29_MOTIF = TCGA / 'LIHC_scRNA_Step29_GRN/Tables/Step29F_MotifPruned_Regulons.csv'
STEP25_CORE = TCGA / 'LIHC_scRNA_Step25_CrossPatient_Robustness/Reference/Step25_FROZEN_MS3_Conserved_Core50.csv'
STEP31_PROGRAMS = TCGA / 'LIHC_scRNA_Step31_Independent_Validation/Audit/Step31E2A_ValidationProgramMapping.csv'
STEP32_SCORE = TCGA / 'LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32EF_score.py'
STEP32_METHOD = TCGA / 'LIHC_scRNA_Step32_Spatial_Validation/Audit/Step32EF_ScoringMethod.json'
STEP33A_DESIGN = AUDIT / 'Step33A_FROZEN_DESIGN.json'
STEP23_NOTEBOOK = TCGA / 'Step23A_Malignant_Only_Data_Audit.ipynb'
STEP29_NOTEBOOK = TCGA / 'Step29_GRN_RegulatoryNetwork.ipynb'
STEP33_NOTEBOOK = TCGA / 'Step33_Virtual_Gene_Knockout.ipynb'

INPUTS = (STEP23_RAW, STEP23A_RAW, STEP26_LOG, STEP29_RAW, STEP29_META,
          STEP29_REG, STEP29_REG_SUM, STEP29_MOTIF, STEP25_CORE,
          STEP31_PROGRAMS, STEP32_SCORE, STEP32_METHOD, STEP33A_DESIGN,
          STEP23_NOTEBOOK, STEP29_NOTEBOOK, STEP33_NOTEBOOK)
FORMAL_TFS = ('FOXM1', 'MYBL2', 'E2F8', 'TFDP1', 'E2F1')
FROZEN_CANDIDATES = FORMAL_TFS + ('BRCA1',)
PRIMARY_ENDPOINTS = ('Core50', 'G2M_Independent', 'Mitotic_Entry_Independent',
                     'Spindle_Chromosome_Independent')
SECONDARY_ENDPOINTS = ('Hepatic_Differentiation', 'Bile_Acid_Metabolism')
PROGRAMS = PRIMARY_ENDPOINTS[1:] + SECONDARY_ENDPOINTS


def read_csv(path):
    with path.open(newline='', encoding='utf-8-sig') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, records, columns):
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(records)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def decode_index(dataset):
    return [x.decode('utf-8') if isinstance(x, bytes) else str(x) for x in dataset[...]]


def frame_index(frame):
    key = frame.attrs.get('_index', '_index')
    if isinstance(key, bytes):
        key = key.decode('utf-8')
    return decode_index(frame[key])


def numeric_audit(dataset, block=2_000_000):
    n = len(dataset)
    minimum = float('inf')
    maximum = float('-inf')
    finite = True
    integer_like = True
    for start in range(0, n, block):
        values = np.asarray(dataset[start:start + block], dtype=np.float64)
        if not len(values):
            continue
        finite &= bool(np.isfinite(values).all())
        minimum = min(minimum, float(values.min()))
        maximum = max(maximum, float(values.max()))
        integer_like &= bool(np.all(np.abs(values - np.rint(values)) <= 1e-8))
    return finite, minimum, maximum, integer_like


def raw_candidate(path, frozen_cells, step29_genes, role, full_numeric=True):
    with h5py.File(path, 'r') as h5:
        x = h5['X']
        shape = tuple(map(int, x.attrs['shape'])) if isinstance(x, h5py.Group) else x.shape
        cells = frame_index(h5['obs'])
        genes = frame_index(h5['var'])
        sparse_type = x.attrs.get('encoding-type', 'dense')
        if isinstance(sparse_type, bytes):
            sparse_type = sparse_type.decode('utf-8')
        if full_numeric:
            data = x['data'] if isinstance(x, h5py.Group) else x
            finite, minimum, maximum, integer_like = numeric_audit(data)
        else:
            finite = minimum = maximum = integer_like = 'NOT_CHECKED_NONRAW_REFERENCE'
        overlap = set(cells) & frozen_cells
        counts_layer = 'raw' in h5 or ('layers' in h5 and any('count' in k.lower() for k in h5['layers']))
        gene_system = 'gene_symbol' if any(g in genes for g in ('FOXM1', 'PLK1')) else 'UNDETERMINED'
        result = {
            'source_role': role, 'path': str(path), 'n_cells': shape[0], 'n_genes': shape[1],
            'X_encoding': sparse_type, 'X_dtype': str((x['data'] if isinstance(x, h5py.Group) else x).dtype),
            'separate_raw_or_counts_slot': counts_layer,
            'nonnegative': minimum >= 0 if full_numeric else 'NOT_CHECKED',
            'finite': finite, 'integer_like': integer_like,
            'minimum': minimum, 'maximum': maximum,
            'barcode_unique': len(cells) == len(set(cells)),
            'frozen_barcode_overlap': len(overlap),
            'exact_barcode_match_count': sum(c in frozen_cells for c in cells),
            'missing_frozen_cells': len(frozen_cells - set(cells)),
            'extra_cells': len(set(cells) - frozen_cells),
            'gene_identifier_system': gene_system,
            'gene_unique': len(genes) == len(set(genes)),
            'exact_gene_overlap_with_Step29': len(set(genes) & step29_genes),
            'true_raw_counts': bool(full_numeric and finite and minimum >= 0 and integer_like),
        }
    return result, cells, genes


def load_csr(path):
    with h5py.File(path, 'r') as h5:
        x = h5['X']
        assert x.attrs['encoding-type'] in ('csr_matrix', b'csr_matrix')
        shape = tuple(map(int, x.attrs['shape']))
        matrix = sparse.csr_matrix((x['data'][...], x['indices'][...], x['indptr'][...]), shape=shape)
        return matrix


def raw_audit(frozen_cells):
    with h5py.File(STEP29_RAW, 'r') as h5:
        step29_genes = set(frame_index(h5['var']))
    descriptions = ((STEP23_RAW, 'exact_full_gene_raw_source', True),
                    (STEP29_RAW, 'frozen_filtered_raw_GRN_input', True),
                    (STEP23A_RAW, 'upstream_malignant_raw_ancestor', True),
                    (STEP26_LOG, 'log_normalized_reference_not_raw', True))
    result_rows = []
    indices = {}
    for path, role, numeric in descriptions:
        result, cells, genes = raw_candidate(path, frozen_cells, step29_genes, role, numeric)
        result_rows.append(result)
        indices[role] = (cells, genes)
    full_cells, full_genes = indices['exact_full_gene_raw_source']
    grn_cells, grn_genes = indices['frozen_filtered_raw_GRN_input']
    assert len(frozen_cells) == 12057
    assert len(full_cells) == len(grn_cells) == len(frozen_cells)
    assert full_cells == grn_cells, 'Step23 and Step29 raw barcode order differs'
    assert set(full_cells) == frozen_cells, 'Step23 raw does not cover frozen cells exactly'
    assert set(grn_genes) <= set(full_genes)
    matrix23 = load_csr(STEP23_RAW)
    matrix29 = load_csr(STEP29_RAW)
    gene_index = {g: i for i, g in enumerate(full_genes)}
    subset = matrix23[:, [gene_index[g] for g in grn_genes]]
    n_count_mismatches = int((subset != matrix29).nnz)
    del matrix23, matrix29, subset
    for r in result_rows:
        r['exact_count_match_to_Step23_on_Step29_genes'] = (
            n_count_mismatches == 0 if r['source_role'] in
            ('exact_full_gene_raw_source', 'frozen_filtered_raw_GRN_input') else 'NOT_APPLICABLE')
        r['count_value_mismatches'] = n_count_mismatches if r['source_role'] in (
            'exact_full_gene_raw_source', 'frozen_filtered_raw_GRN_input') else ''
    full, grn = result_rows[0], result_rows[1]
    passed = bool(full['true_raw_counts'] and grn['true_raw_counts']
                  and full['barcode_unique'] and grn['barcode_unique']
                  and full['missing_frozen_cells'] == grn['missing_frozen_cells'] == 0
                  and full['extra_cells'] == grn['extra_cells'] == 0
                  and n_count_mismatches == 0)
    return result_rows, passed, n_count_mismatches, set(full_genes), step29_genes


def topology_audit(expression_genes, step29_genes):
    source = json.loads(STEP29_REG.read_text(encoding='utf-8'))
    summary = {r['Regulon']: r for r in read_csv(STEP29_REG_SUM)}
    assert set(source) == set(summary), 'Regulon JSON and frozen summary disagree'
    edges = []
    pair_counts = Counter()
    for regulon, targets in source.items():
        assert regulon.endswith('(+)'), f'Unexpected regulon name: {regulon}'
        tf = regulon[:-3]
        assert int(summary[regulon]['N_targets']) == len(targets)
        for target in targets:
            pair_counts[(tf, target)] += 1
            edges.append({
                'TF': tf, 'target_gene': target, 'source_regulon': regulon,
                'source_step': 'Step29F motif-filtered regulon',
                'motif_filtered': True, 'edge_in_frozen_topology': True,
                'exact_expression_gene_match': target in expression_genes and tf in expression_genes,
                'exact_Step29_gene_match': target in step29_genes and tf in step29_genes,
            })
    total_unique = len(pair_counts)
    duplicate_count = sum(n - 1 for n in pair_counts.values())
    self_edges = sum(1 for tf, target in pair_counts if tf == target)
    tf_set = {tf for tf, _ in pair_counts}
    mapped = {pair for pair in pair_counts if pair[0] in expression_genes and pair[1] in expression_genes}
    rows = [{
        'scope': 'ALL', 'TF': 'ALL', 'regulons': len(source), 'total_TFs': len(tf_set),
        'source_edge_rows': len(edges), 'unique_TF_target_edges': total_unique,
        'duplicate_edge_rows': duplicate_count, 'self_edges': self_edges,
        'frozen_targets': len({target for _, target in pair_counts}),
        'retained_exact_expression_mapping': len(mapped),
        'retained_percentage': 100 * len(mapped) / total_unique if total_unique else 0,
        'TF_present_in_expression': 'ALL' if tf_set <= expression_genes else False,
        'no_new_edges': True,
    }]
    for tf in FORMAL_TFS:
        pairs = {pair for pair in pair_counts if pair[0] == tf}
        present = tf in tf_set
        mapped_tf = {pair for pair in pairs if pair in mapped}
        rows.append({
            'scope': 'FORMAL_TF', 'TF': tf, 'regulons': sum(r.startswith(tf + '(') for r in source),
            'total_TFs': '', 'source_edge_rows': sum(n for pair, n in pair_counts.items() if pair[0] == tf),
            'unique_TF_target_edges': len(pairs),
            'duplicate_edge_rows': sum(n - 1 for pair, n in pair_counts.items() if pair[0] == tf),
            'self_edges': sum(pair[0] == pair[1] for pair in pairs),
            'frozen_targets': len({pair[1] for pair in pairs}),
            'retained_exact_expression_mapping': len(mapped_tf),
            'retained_percentage': 100 * len(mapped_tf) / len(pairs) if pairs else 0,
            'TF_present_in_expression': present and tf in expression_genes,
            'no_new_edges': True,
        })
    conversion_pass = (len(source) > 0 and len(tf_set) > 0
                       and all(r['TF_present_in_expression'] and r['retained_percentage'] == 100
                               for r in rows[1:]) and len(mapped) == total_unique)
    return rows, edges, conversion_pass


def endpoint_audit():
    core = [r['Gene'] for r in read_csv(STEP25_CORE)]
    mapping = read_csv(STEP31_PROGRAMS)
    sets = {'Core50': core}
    for program in PROGRAMS:
        sets[program] = [r['requested_symbol'] for r in mapping
                         if r['program'] == program and r['mapping_status'] == 'retained']
    expected = {'Core50': 50, 'G2M_Independent': 22,
                'Mitotic_Entry_Independent': 16,
                'Spindle_Chromosome_Independent': 29,
                'Hepatic_Differentiation': 19, 'Bile_Acid_Metabolism': 18}
    for endpoint, genes in sets.items():
        assert len(genes) == expected[endpoint] and len(set(genes)) == len(genes), endpoint
    rows = []
    for tf in FORMAL_TFS:
        for endpoint in PRIMARY_ENDPOINTS + SECONDARY_ENDPOINTS:
            genes = sets[endpoint]
            present = tf in genes
            rows.append({
                'target': tf, 'endpoint': endpoint,
                'target_present_in_endpoint': present,
                'endpoint_size_original': len(genes),
                'endpoint_size_leave_target_out': len(genes) - int(present),
                'primary_perturbation_score_rule': 'leave_target_out_only',
                'secondary_descriptive_score_rule': 'full_frozen_set',
                'gene_set_source': str(STEP25_CORE if endpoint == 'Core50' else STEP31_PROGRAMS),
            })
    return rows, sets


def scoring_audit():
    method = json.loads(STEP32_METHOD.read_text(encoding='utf-8'))
    source = STEP32_SCORE.read_text(encoding='utf-8')
    required = ('np.log1p((x / totals[:, None]) * 1e4)',
                'raw.std(axis=0, ddof=0)', 'z[:, cols].mean(axis=1)')
    assert all(snippet in source for snippet in required)
    return {
        'status': 'ENDPOINT_SCORING_DECISION_REQUIRED',
        'exact_previous_scoring_source': str(STEP32_SCORE),
        'source_sha256': sha256(STEP32_SCORE),
        'method_record': str(STEP32_METHOD),
        'frozen_Core50_source': str(STEP25_CORE),
        'frozen_program_source': str(STEP31_PROGRAMS),
        'previous_definition': {
            'normalization': method['normalization'],
            'gene_standardization': method['gene_standardization'],
            'program_score': method['program_score'],
            'zero_variance_rule': 'contribute zero',
        },
        'CellOracle_documented_input': 'normalized and log-transformed expression; not scaled or centered',
        'CellOracle_documented_simulation_output': 'simulated_count layer from imputed expression and fitted coefficients; not verified as raw integer counts',
        'directly_transferable_to_baseline_and_simulated_layers': False,
        'reason': ('The frozen Step32 formula starts from raw integer counts and uses a fixed '
                   'log1p(count/total*10000) transformation followed by global gene z scores. '
                   'CellOracle simulation returns an imputed simulated_count layer, not verified '
                   'as raw counts. A baseline reference distribution, output scale, and paired '
                   'score application must be explicitly locked before any endpoint calculation.'),
        'primary_KO_endpoint_rule': 'remove only the perturbed TF if present in that frozen gene set',
        'secondary_descriptive_endpoint_rule': 'use full frozen gene set',
        'official_source': 'https://morris-lab.github.io/CellOracle.documentation/modules/celloracle.html',
    }


def environment_plan(raw_pass, topology_pass, scoring_pass):
    def version(name):
        try:
            return metadata.version(name)
        except metadata.PackageNotFoundError:
            return 'NOT_INSTALLED'
    env_dir = Path('/home/LOCAL_USER/miniforge3/envs/ko33_wsl')
    celloracle_status = version('celloracle')
    eligible_to_create = bool(raw_pass and topology_pass and scoring_pass)
    text = [
        'STEP33B0 — ISOLATED CELLORACLE ENVIRONMENT PLAN',
        f'Audit interpreter: {sys.executable}', f'Audit Python: {sys.version.split()[0]}',
        'Audit conda environment: spatial32_wsl (frozen; do not modify)',
        f'Audit scanpy: {version("scanpy")}', f'Audit anndata: {version("anndata")}',
        'Intended new environment: ko33_wsl (Python 3.10)',
        'Intended CellOracle version: 0.18.0',
        f'ko33_wsl exists: {env_dir.is_dir()}',
        f'CellOracle installed in audit environment: {celloracle_status}',
        f'Prerequisites 3–6 all PASS: {eligible_to_create}',
        'Environment created by Step33B0: False',
        'CellOracle installed by Step33B0: False',
        'No Step27–32 environment modified.',
        'Official install documentation: https://morris-lab.github.io/CellOracle.documentation/installation/index.html',
        'Official documentation notes that WSL execution can be very slow.',
    ]
    if not eligible_to_create:
        text.append('Blocking reason: endpoint scoring transfer is not yet validly locked; do not create or install ko33_wsl.')
    return '\n'.join(text) + '\n', env_dir.is_dir(), celloracle_status


def main():
    assert TCGA.is_dir(), 'LOCAL_FILESYSTEM_UNAVAILABLE'
    missing = [str(p) for p in INPUTS if not p.is_file()]
    if missing:
        raise FileNotFoundError('Required frozen inputs missing: ' + '; '.join(missing))
    assert AUDIT.is_dir()
    frozen_meta = read_csv(STEP29_META)
    frozen_cells = {r['Cell_ID'] for r in frozen_meta}
    assert len(frozen_meta) == len(frozen_cells) == 12057
    assert sum(r['final_malignant_state_id'] == 'MS3' for r in frozen_meta) == 720
    assert sum(r['final_malignant_state_id'] == 'MS6' for r in frozen_meta) == 557
    assert all(r['patient'] and r['final_malignant_state_id'] for r in frozen_meta)
    old_design = json.loads(STEP33A_DESIGN.read_text(encoding='utf-8'))
    assert set(old_design['primary_targets']) == set(FROZEN_CANDIDATES)

    raw_rows, raw_pass, mismatches, expression_genes, step29_genes = raw_audit(frozen_cells)
    topology_rows, conversion_rows, topology_pass = topology_audit(expression_genes, step29_genes)
    endpoint_rows, endpoint_sets = endpoint_audit()
    scoring = scoring_audit()
    scoring_pass = scoring['directly_transferable_to_baseline_and_simulated_layers']
    env_text, ko_env_exists, celloracle_status = environment_plan(raw_pass, topology_pass, scoring_pass)

    raw_file = AUDIT / 'Step33B0_RawCountProvenance.csv'
    topo_file = AUDIT / 'Step33B0_FrozenTopologyAudit.csv'
    base_file = AUDIT / 'Step33B0_CellOracleBaseGRN.csv'
    endpoint_file = AUDIT / 'Step33B0_TargetEndpointOverlap.csv'
    score_file = AUDIT / 'Step33B0_ScoringProvenance.json'
    lock_file = AUDIT / 'Step33B0_MethodLock.json'
    env_file = AUDIT / 'Step33B0_EnvironmentPlan.txt'
    write_csv(raw_file, raw_rows, tuple(raw_rows[0]))
    write_csv(topo_file, topology_rows, tuple(topology_rows[0]))
    write_csv(base_file, conversion_rows, tuple(conversion_rows[0]))
    write_csv(endpoint_file, endpoint_rows, tuple(endpoint_rows[0]))
    write_json(score_file, scoring)
    env_file.write_text(env_text, encoding='utf-8')

    fox = next(r for r in topology_rows if r['TF'] == 'FOXM1')
    pilot_ready = bool(raw_pass and topology_pass and scoring_pass
                       and fox['TF_present_in_expression']
                       and fox['retained_percentage'] == 100 and ko_env_exists
                       and celloracle_status != 'NOT_INSTALLED')
    blockers = []
    if not raw_pass:
        blockers.append('CELLORACLE_RAW_INPUT_NOT_VERIFIED')
    if not topology_pass:
        blockers.append('FROZEN_TOPOLOGY_EXACT_MAPPING_NOT_VERIFIED')
    if not scoring_pass:
        blockers.append('ENDPOINT_SCORING_DECISION_REQUIRED')
    if not ko_env_exists or celloracle_status == 'NOT_INSTALLED':
        blockers.append('ko33_wsl/CellOracle 0.18.0 not created or installed because prerequisite scoring decision is unresolved')
    lock = {
        'step': 'Step33B0', 'status': 'COMPLETE_PILOT_NOT_READY',
        'method': 'CellOracle-based in silico TF perturbation',
        'frozen_regulatory_topology_source': str(STEP29_REG),
        'frozen_motif_pruned_regulon_table': str(STEP29_MOTIF),
        'frozen_regulon_summary': str(STEP29_REG_SUM),
        'CellOracleBaseGRN_status': 'UNFITTED_FROZEN_INPUT_CANDIDATE',
        'CellOracle_import_topology_preservation_verified': False,
        'CellOracle_self_edge_handling_note': ('The frozen topology contains 140 self edges, including MYBL2, E2F8, TFDP1, and E2F1. Their handling by CellOracle 0.18.0 has not been tested; none were silently removed.'),
        'raw_GRNBoost2_edges_used_as_perturbation_coefficients': False,
        'frozen_candidate_list': list(FROZEN_CANDIDATES),
        'formal_primary_CellOracle_TF_KO_set': list(FORMAL_TFS),
        'BRCA1_formal_CellOracle_TF_KO_eligible': False,
        'BRCA1_exclusion_reason': 'Computational outgoing edges/regulon-like representation do not establish TF-specific perturbation eligibility.',
        'network_topology_rediscovery_allowed': False,
        'new_regulatory_edges_allowed': False,
        'frozen_topology_modification_allowed': False,
        'simulation_coefficient_fitting_within_frozen_topology_allowed': True,
        'raw_count_provenance_pass': raw_pass,
        'exact_12057_frozen_cell_mapping': raw_pass,
        'Step29_count_value_mismatches_vs_Step23_raw': mismatches,
        'frozen_topology_conversion_pass': topology_pass,
        'frozen_topology_exact_gene_mapping_pass': topology_pass,
        'frozen_regulons': int(topology_rows[0]['regulons']),
        'frozen_unique_TF_target_edges': int(topology_rows[0]['unique_TF_target_edges']),
        'FOXM1_frozen_targets': int(fox['frozen_targets']),
        'FOXM1_exactly_mapped_targets': int(fox['retained_exact_expression_mapping']),
        'endpoint_scoring_transfer_pass': scoring_pass,
        'endpoint_scoring_decision_required': not scoring_pass,
        'primary_perturbation_endpoint_rule': 'leave-target-out only for target-containing frozen endpoints',
        'secondary_descriptive_endpoint_rule': 'full frozen gene sets',
        'KO_condition_prelocked': {'target_expression': 0, 'n_propagation_primary': 3,
                                   'n_propagation_sensitivity_later': [1, 2]},
        'Step34_transition_or_vector_field_analysis_allowed_here': False,
        'FOXM1_PILOT_READY': 'YES' if pilot_ready else 'NO',
        'pilot_blocking_reasons': blockers,
        'ko33_wsl_exists': ko_env_exists,
        'CellOracle_installed_in_audit_environment': celloracle_status,
        'formal_knockout_performed': False,
        'new_topology_discovered': False, 'new_regulatory_edges_added': False,
        'targets_reselected': False, 'original_RiskScore_reconstructed': False,
        'Step34_started': False, 'Step35_started': False,
        'outputs': [str(p) for p in (raw_file, topo_file, base_file, endpoint_file,
                                    score_file, lock_file, env_file)],
    }
    write_json(lock_file, lock)
    print('STEP33B0 — CELLORACLE INPUT / TOPOLOGY AUDIT')
    print('[1] Raw-count provenance')
    print(f'  Source: {STEP23_RAW}')
    print(f'  true raw counts: {raw_pass}; frozen cells matched: {raw_rows[0]["exact_barcode_match_count"]}/12057; barcode mismatches: {raw_rows[0]["missing_frozen_cells"]}; integer-like: {raw_rows[0]["integer_like"]}; count mismatches vs Step29: {mismatches}')
    print('[2] Frozen Step29 topology')
    print(f'  Regulons: {topology_rows[0]["regulons"]}; TFs: {topology_rows[0]["total_TFs"]}; unique edges: {topology_rows[0]["unique_TF_target_edges"]}; duplicate rows: {topology_rows[0]["duplicate_edge_rows"]}; self edges: {topology_rows[0]["self_edges"]}')
    for r in topology_rows[1:]:
        print(f'  {r["TF"]}: {r["frozen_targets"]} frozen targets; {r["retained_exact_expression_mapping"]} exact mapped; {r["retained_percentage"]:.2f}% retained')
    print('[3] Endpoint circularity')
    for tf in FORMAL_TFS:
        entries = [r for r in endpoint_rows if r['target'] == tf]
        print('  ' + tf + ': ' + '; '.join(f'{r["endpoint"]}={r["endpoint_size_leave_target_out"]}/{r["endpoint_size_original"]} (overlap={r["target_present_in_endpoint"]})' for r in entries))
    print('[4] Scoring provenance')
    print(f'  Source: {STEP32_SCORE}; directly transferable: {scoring_pass}; status: {scoring["status"]}')
    print('[5] Method lock: CellOracle; frozen Step29 motif-filtered topology; predictive coefficient fitting within that topology only')
    print(f'[6] Environment plan: ko33_wsl exists={ko_env_exists}; CellOracle installed={celloracle_status != "NOT_INSTALLED"}; no environment created or package installed')
    print('[7] Firewall')
    print('Formal knockout performed: False')
    print('New topology discovered: False')
    print('New regulatory edges added: False')
    print('Targets reselected: False')
    print('Original RiskScore reconstructed: False')
    print('Step34 started: False')
    print('Step35 started: False')
    print('[8] Decision')
    print(f'FOXM1_PILOT_READY = {lock["FOXM1_PILOT_READY"]}')
    if blockers:
        print('Blocking reason: ' + '; '.join(blockers))
    print('STEP33B0: COMPLETE')
    print('NEXT: return report to ChatGPT before formal FOXM1 pilot.')


if __name__ == '__main__':
    main()
