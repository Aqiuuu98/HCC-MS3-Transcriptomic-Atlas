#!/usr/bin/env python3
"""Step33B1: pre-fit topology compatibility gate for the FOXM1 pilot.

Input: Step33B0_CellOracleBaseGRN.csv; Step29D frozen metadata and gene list;
       installed CellOracle 0.18.0 oracle_GRN.py.
Output: Step33B1 environment, topology, input, model, simulation, scoring,
        and method-lock audits under Step33/Audit.
Conda environment: ko33_wsl.
Python version: 3.10.21.
Frozen dependencies: celloracle==0.18.0, numpy==1.26.4,
                     Cython==0.29.37; remaining package versions recorded
                     in Step33B1_EnvironmentFreeze.txt.

This script enforces the pre-specified stop if CellOracle removes self edges.
It does not fit coefficients or simulate a knockout after that stop.
"""
from hcc_release_paths import project_path as _hcc_path

import ast
import csv
import hashlib
import importlib.metadata as md
import json
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(_hcc_path('LIHC_scRNA_Step33_Virtual_Knockout'))
AUDIT = ROOT / 'Audit'
TOPOLOGY = AUDIT / 'Step33B0_CellOracleBaseGRN.csv'
META = Path(_hcc_path('LIHC_scRNA_Step29_GRN/Tables/Step29D_Frozen_Cell_Metadata.csv'))
GENES = Path(_hcc_path('LIHC_scRNA_Step29_GRN/Tables/Step29D_GRN_Genes.csv'))
RAW = Path(_hcc_path('LIHC_scRNA_Step23_Malignant_State_Analysis/Data/GSE149614_STEP23D3_LINEAGE_CLEAN_RAW_COUNTS.h5ad'))
ORACLE_SOURCE = Path(sys.prefix) / 'lib/python3.10/site-packages/celloracle/trajectory/oracle_GRN.py'
EXPECTED_EDGES = 32530
EXPECTED_SELF_EDGES = 140
EXPECTED_CELLS = 12057
EXPECTED_GENES = 12916


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for block in iter(lambda: fh.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_csv(path):
    with path.open(newline='', encoding='utf-8-sig') as fh:
        return list(csv.DictReader(fh))


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')


def version(name):
    try:
        return md.version(name)
    except md.PackageNotFoundError:
        return 'NOT_INSTALLED'


def process_output(args):
    p = subprocess.run(args, text=True, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, check=False)
    return f'$ {" ".join(args)}\nexit_code={p.returncode}\n{p.stdout}'


def environment_freeze():
    packages = ('celloracle', 'scanpy', 'anndata', 'numpy', 'pandas',
                'scipy', 'scikit-learn', 'Cython')
    lines = [
        'STEP33B1 ENVIRONMENT RECORD',
        'Environment: ko33_wsl',
        'Installation scope: isolated ko33_wsl only; frozen Step27-32 environments untouched.',
        'Before installation: Python 3.10.21; celloracle, scanpy, anndata, numpy, pandas, scipy, scikit-learn NOT_INSTALLED.',
        'Before pip list: packaging==26.3, pip==26.2.1, setuptools==84.0.0, wheel==0.48.0.',
        'First pip install celloracle==0.18.0 failed while building velocyto: missing numpy in build isolation.',
        'Dependency-resolving retry was stopped before installation because it pulled unrelated Jupyter/motif packages.',
        'Final installation: numpy==1.26.4, Cython==0.29.37, celloracle==0.18.0 --no-deps.',
        'This is a topology-gate environment, not a complete runnable simulation environment.',
        f'Python runtime: {sys.version.split()[0]}',
        *[f'{p}: {version(p)}' for p in packages],
        process_output([sys.executable, '--version']),
        process_output([sys.executable, '-m', 'pip', 'freeze']),
        process_output(['/home/LOCAL_USER/miniforge3/bin/conda', 'list', '-n', 'ko33_wsl']),
    ]
    (AUDIT / 'Step33B1_EnvironmentFreeze.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return {p: version(p) for p in packages}


def source_self_edge_behavior():
    assert ORACLE_SOURCE.is_file(), ORACLE_SOURCE
    src = ORACLE_SOURCE.read_text(encoding='utf-8')
    tree = ast.parse(src)
    fit = next((x for x in tree.body if isinstance(x, ast.FunctionDef)
                and x.name == '_getCoefMatrix'), None)
    assert fit is not None, 'CellOracle _getCoefMatrix not found'
    matches = []
    for node in ast.walk(fit):
        if not isinstance(node, ast.If):
            continue
        test = node.test
        if not (isinstance(test, ast.Compare) and isinstance(test.left, ast.Name)
                and test.left.id == 'target_gene' and len(test.ops) == 1
                and isinstance(test.ops[0], ast.In) and len(test.comparators) == 1
                and isinstance(test.comparators[0], ast.Name)
                and test.comparators[0].id == 'reggenes'):
            continue
        for child in node.body:
            if not isinstance(child, ast.Expr) or not isinstance(child.value, ast.Call):
                continue
            call = child.value
            if (isinstance(call.func, ast.Attribute) and call.func.attr == 'remove'
                    and isinstance(call.func.value, ast.Name)
                    and call.func.value.id == 'reggenes' and len(call.args) == 1
                    and isinstance(call.args[0], ast.Name)
                    and call.args[0].id == 'target_gene'):
                matches.append({'if_line': node.lineno, 'remove_line': child.lineno,
                                'source': ast.get_source_segment(src, node)})
    return matches


def main():
    assert version('celloracle') == '0.18.0', 'Wrong CellOracle version'
    assert sys.version_info[:2] == (3, 10), 'Wrong Python version'
    versions = environment_freeze()
    rows = read_csv(TOPOLOGY)
    edges = [(r['TF'], r['target_gene']) for r in rows]
    assert len(edges) == EXPECTED_EDGES
    assert len(set(edges)) == EXPECTED_EDGES
    self_edges = sorted((tf, target) for tf, target in edges if tf == target)
    assert len(self_edges) == EXPECTED_SELF_EDGES
    genes = [r['Gene'] for r in read_csv(GENES)]
    assert len(genes) == EXPECTED_GENES and len(set(genes)) == EXPECTED_GENES
    meta = read_csv(META)
    assert len(meta) == EXPECTED_CELLS
    assert len({r['Cell_ID'] for r in meta}) == EXPECTED_CELLS
    assert RAW.is_file(), RAW
    # CellOracle's TFdict is target -> candidate TF list. Preserve all edges here.
    tfdict = defaultdict(list)
    for tf, target in edges:
        tfdict[target].append(tf)
    assert sum(map(len, tfdict.values())) == EXPECTED_EDGES
    assert sum(target in tfs for target, tfs in tfdict.items()) == EXPECTED_SELF_EDGES
    matches = source_self_edge_behavior()
    if len(matches) != 1:
        raise RuntimeError('Self-edge handling uncertain; STOP before fitting')
    patient_counts = Counter(r['patient'] for r in meta)
    state_counts = Counter(r['final_malignant_state_id'] for r in meta)
    canonical = '\n'.join(f'{tf}\t{target}' for tf, target in sorted(edges)) + '\n'
    topology = {
        'step': 'Step33B1', 'source_path': str(TOPOLOGY),
        'source_file_sha256': sha256(TOPOLOGY),
        'canonical_sorted_edge_tsv_sha256': hashlib.sha256(canonical.encode()).hexdigest(),
        'expected_unique_edges': EXPECTED_EDGES,
        'input_unique_edges': len(set(edges)),
        'TFdict_edge_count_before_fitting': sum(map(len, tfdict.values())),
        'input_self_edges': len(self_edges),
        'self_edge_examples': self_edges[:10],
        'celloracle_source_path': str(ORACLE_SOURCE),
        'celloracle_source_sha256': sha256(ORACLE_SOURCE),
        'celloracle_version': version('celloracle'),
        'coefficient_function': '_getCoefMatrix',
        'self_edge_removal_source_line': matches[0]['remove_line'],
        'self_edge_removal_source': matches[0]['source'],
        'CellOracle_import_TF_data_tested': False,
        'TFdict_conversion_pass': True,
        'fit_preserves_all_self_edges': False,
        'reason': 'CellOracle 0.18.0 _getCoefMatrix removes target_gene from reggenes before Ridge fit.',
        'fitted_topology_edge_count': None,
        'status': 'STOP_INTERNAL_SELF_EDGE_REMOVAL',
    }
    write_json(AUDIT / 'Step33B1_TopologyChecksum.json', topology)
    input_rows = [
        {'item': 'Step23 lineage-clean raw counts', 'path': str(RAW), 'observed': 'file exists; B0 validated exact counts', 'status': 'B0_VERIFIED_NOT_REIMPORTED'},
        {'item': 'frozen cell metadata', 'path': str(META), 'observed': len(meta), 'status': 'PASS'},
        {'item': 'frozen Step29 gene universe', 'path': str(GENES), 'observed': len(genes), 'status': 'PASS'},
        {'item': 'frozen patients', 'path': str(META), 'observed': len(patient_counts), 'status': 'PASS'},
        {'item': 'frozen MS3 cells', 'path': str(META), 'observed': state_counts['MS3'], 'status': 'PASS'},
        {'item': 'CellOracle raw-count import', 'path': '', 'observed': 'not run: topology gate stopped', 'status': 'NOT_RUN'},
    ]
    with (AUDIT / 'Step33B1_CellOracleInputAudit.csv').open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=('item', 'path', 'observed', 'status'))
        writer.writeheader(); writer.writerows(input_rows)
    with (AUDIT / 'Step33B1_EndpointCoverageAudit.csv').open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=('check', 'baseline_layer', 'status', 'reason'))
        writer.writeheader()
        writer.writerow({'check': 'all frozen primary and secondary endpoint genes',
                         'baseline_layer': 'imputed_count',
                         'status': 'NOT_RUN_TOPOLOGY_GATE',
                         'reason': topology['reason']})
    write_json(AUDIT / 'Step33B1_ModelFitAudit.json', {
        'status': 'NOT_RUN_TOPOLOGY_GATE', 'reason': topology['reason'],
        'fitted_states': [], 'fitted_genes': None, 'failed_genes': None,
        'coefficient_matrix_dimensions': None, 'FOXM1_coefficient_available_in_MS3': None,
        'fitted_oracle_object_saved': False,
    })
    write_json(AUDIT / 'Step33B1_SimulationQCAudit.json', {
        'status': 'NOT_RUN_TOPOLOGY_GATE', 'FOXM1_simulated': False,
        'distribution_range_evaluated': False, 'distribution_plots_created': False,
    })
    write_json(AUDIT / 'Step33B1_ScoringParameters.json', {
        'status': 'LOCKED_NOT_EXECUTED',
        'name': 'Baseline-Anchored CellOracle Perturbation Score',
        'baseline_layer': 'imputed_count', 'perturbed_layer': 'simulated_count',
        'reference_cells': EXPECTED_CELLS, 'shared_baseline_mean_and_sd': True,
        'primary_Core50_endpoint': 'Core50-LTO-FOXM1',
        'primary_Core50_gene_count': 49, 'direct_target_removed': 'FOXM1',
        'n_propagation': 3, 'alpha': 10,
        'k': int(0.025 * EXPECTED_CELLS), 'balanced': True,
        'b_sight': int(0.025 * EXPECTED_CELLS) * 8,
        'b_maxl': int(0.025 * EXPECTED_CELLS) * 4,
        'n_comps': None, 'n_comps_status': 'NOT_COMPUTED_TOPOLOGY_GATE',
    })
    write_json(AUDIT / 'Step33B1_MethodLock.json', {
        'step': 'Step33B1', 'status': 'STOPPED_PRE_FIT',
        'decision': 'FOXM1_PILOT_FAILED',
        'blocking_condition': 'CellOracle 0.18.0 internally removes frozen self-edges during coefficient fitting',
        'frozen_input_edges': EXPECTED_EDGES, 'frozen_input_self_edges': EXPECTED_SELF_EDGES,
        'input_TFdict_edges': EXPECTED_EDGES, 'model_fitted': False,
        'only_FOXM1_simulated': False, 'other_formal_TFs_simulated': False,
        'new_topology_discovered': False, 'new_edges_added': False,
        'target_reselection': False, 'original_RiskScore_reconstructed': False,
        'state_transition_analysis': False, 'Step34_started': False, 'Step35_started': False,
        'environment_versions': versions,
        'patient_counts': dict(patient_counts), 'state_counts': dict(state_counts),
        'manual_decision_required': True,
    })
    print('STEP33B1 - FOXM1 VIRTUAL KO PILOT REPORT')
    print('[1] Environment: ko33_wsl; Python', sys.version.split()[0],
          '; CellOracle', version('celloracle'), '; full freeze: Step33B1_EnvironmentFreeze.txt')
    print('[2] Input: frozen metadata', len(meta), 'cells;', len(genes), 'Step29 genes;',
          len(patient_counts), 'patients; MS3', state_counts['MS3'], 'cells; Step23 raw provenance verified in B0.')
    print('[3] Topology:', len(edges), 'frozen edges,', len(self_edges), 'self-edges;',
          'TFdict retained all', sum(map(len, tfdict.values())), 'edges before fitting.')
    print('[3] STOP: installed CellOracle _getCoefMatrix line', matches[0]['remove_line'],
          'removes each target_gene from its reggenes before coefficient fit.')
    print('[4] Technical preprocessing: NOT RUN. Locked k=301, b_sight=2408, b_maxl=1204; n_comps not determined.')
    print('[5] Model fitting: NOT RUN; no fitted object, coefficients, or failed-gene list.')
    print('[6] Simulation QC: NOT RUN; no FOXM1 direct perturbation check.')
    print('[7] MS3 endpoint effects: NOT CALCULATED; patient-level effects unavailable.')
    print('[8] Firewall: no TF simulated; no new topology, edges, target selection, RiskScore reconstruction, state transition, Step34, or Step35.')
    print('[9] Decision: FOXM1_PILOT_FAILED - explicit self-edge preservation requirement is incompatible with CellOracle 0.18.0 fitting.')
    print('STEP33B1: STOPPED_PRE_FIT')
    print('NEXT: manual method decision on 140 frozen self-edges before any multi-TF simulation.')


if __name__ == '__main__':
    main()
