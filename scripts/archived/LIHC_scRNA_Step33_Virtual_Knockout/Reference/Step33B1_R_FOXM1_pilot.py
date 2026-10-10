#!/usr/bin/env python3
"""Step33B1-R: CellOracle FOXM1 pilot after the locked self-loop decision.

Input: frozen Step33B0 base GRN, Step23 lineage-clean counts, Step29 cells/genes,
       Step25 Core50 and Step31 frozen independent programs.
Output: Step33B1-R derived topology, model, QC, endpoint and notebook audits.
Conda environment: ko33_wsl; Python: 3.10.21.
Frozen dependencies: CellOracle 0.18.0; exact runtime versions are recorded
                     in Step33B1_R_EnvironmentFreeze.txt.

Run phases with --phase topology, --phase fit, or --phase simulate. The original
32,530-edge Step33B0 file and prior Step33B1 stopped record are never edited.
"""
from hcc_release_paths import project_path as _hcc_path

import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
import subprocess
import sys
import warnings
from collections import Counter, defaultdict
from pathlib import Path

os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
os.environ.setdefault('OMP_NUM_THREADS', '2')
os.environ.setdefault('MKL_NUM_THREADS', '2')
os.environ.setdefault('MPLBACKEND', 'Agg')


TCGA = Path(_hcc_path(''))
ROOT = TCGA / 'LIHC_scRNA_Step33_Virtual_Knockout'
AUDIT = ROOT / 'Audit'
FROZEN = AUDIT / 'Step33B0_CellOracleBaseGRN.csv'
DERIVED = AUDIT / 'Step33B1_CellOracleCompatibleTopology.csv'
SELF_DECISION = AUDIT / 'Step33B1_SelfLoopDecision.csv'
DERIVATION = AUDIT / 'Step33B1_TopologyDerivation.json'
FORMAL = AUDIT / 'Step33B1_FormalTF_SelfLoopAudit.csv'
FORMAL_TFS = ('FOXM1', 'MYBL2', 'E2F8', 'TFDP1', 'E2F1')
STEP23_RAW = TCGA / 'LIHC_scRNA_Step23_Malignant_State_Analysis/Data/GSE149614_STEP23D3_LINEAGE_CLEAN_RAW_COUNTS.h5ad'
STEP29_META = TCGA / 'LIHC_scRNA_Step29_GRN/Tables/Step29D_Frozen_Cell_Metadata.csv'
STEP29_GENES = TCGA / 'LIHC_scRNA_Step29_GRN/Tables/Step29D_GRN_Genes.csv'
STEP26 = TCGA / 'LIHC_scRNA_Step26_MS3_Biology/Data/GSE149614_STEP26_FINAL_FROZEN_MS3_BIOLOGY.h5ad'
PREFIT = ROOT / 'Data/Step33B1_R_CellOracle_PreFit.celloracle.oracle'
FITTED = ROOT / 'Data/Step33B1_CellOracle_FittedModel.celloracle.oracle'
FIT_AUDIT = AUDIT / 'Step33B1_ModelFitAudit.json'
INPUT_AUDIT = AUDIT / 'Step33B1_CellOracleInputAudit.csv'
SCORING = AUDIT / 'Step33B1_ScoringParameters.json'
SIM_QC = AUDIT / 'Step33B1_SimulationQCAudit.json'
ENDPOINT_AUDIT = AUDIT / 'Step33B1_EndpointCoverageAudit.csv'
CORE50 = TCGA / 'LIHC_scRNA_Step25_CrossPatient_Robustness/Reference/Step25_FROZEN_MS3_Conserved_Core50.csv'
PROGRAMS = TCGA / 'LIHC_scRNA_Step31_Independent_Validation/Audit/Step31E2A_ValidationProgramMapping.csv'
PROGRAM_DELTA = ROOT / 'Tables/Step33B1_FOXM1_MS3_ProgramDelta.csv'
PATIENT_SUMMARY = ROOT / 'Tables/Step33B1_FOXM1_PatientSummary.csv'
ENDPOINT_SUMMARY = ROOT / 'Tables/Step33B1_FOXM1_EndpointSummary.csv'
FIGURE_PNG = ROOT / 'Figures/Step33B1_FOXM1_SimulationDistribution.png'
FIGURE_PDF = ROOT / 'Figures/Step33B1_FOXM1_SimulationDistribution.pdf'
EXPECTED_CELLS = 12057
EXPECTED_GENES = 12916
PCA_EVALUATION_COMPONENTS = 100
KNN_K = int(0.025 * EXPECTED_CELLS)
KNN_B_SIGHT = KNN_K * 8
KNN_B_MAXL = KNN_K * 4
ALPHA = 10


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read_csv(path):
    with path.open(newline='', encoding='utf-8-sig') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows, fieldnames):
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')


def canonical_hash(edges):
    payload = ''.join(f'{tf}\t{target}\n' for tf, target in sorted(edges))
    return hashlib.sha256(payload.encode('utf-8')).hexdigest()


def derive_topology():
    original = read_csv(FROZEN)
    original_edges = [(r['TF'], r['target_gene']) for r in original]
    assert len(original) == 32530 and len(set(original_edges)) == 32530
    removed = [edge for edge in original_edges if edge[0] == edge[1]]
    retained = [edge for edge in original_edges if edge[0] != edge[1]]
    assert len(removed) == 140 and len(retained) == 32390
    assert set(original_edges) == set(retained) | set(removed)
    assert not set(retained) & set(removed)
    decision_rows = [
        {'tf': tf, 'target_gene': target,
         'is_self_loop': str(tf == target),
         'original_frozen_topology': 'True',
         'retained_in_celloracle_topology': str(tf != target),
         'removal_reason': 'CELLORACLE_SELF_REGULATION_UNSUPPORTED' if tf == target else ''}
        for tf, target in original_edges]
    write_csv(SELF_DECISION, decision_rows,
              ('tf', 'target_gene', 'is_self_loop', 'original_frozen_topology',
               'retained_in_celloracle_topology', 'removal_reason'))
    derived_rows = [{'tf': tf, 'target_gene': target} for tf, target in retained]
    write_csv(DERIVED, derived_rows, ('tf', 'target_gene'))
    reloaded = [(r['tf'], r['target_gene']) for r in read_csv(DERIVED)]
    assert set(reloaded) == set(retained) and len(reloaded) == 32390
    original_set, simulation_set = set(original_edges), set(reloaded)
    assert original_set - simulation_set == set(removed)
    assert not simulation_set - original_set
    summary = {
        'step': 'Step33B1-R',
        'decision': 'METHOD_REQUIRED_COMPATIBILITY_TRANSFORMATION',
        'original_frozen_topology_path': str(FROZEN),
        'simulation_topology_path': str(DERIVED),
        'original_edges': len(original_edges),
        'self_loops': len(removed),
        'simulation_edges': len(reloaded),
        'added_edges': len(simulation_set - original_set),
        'removed_nonself_edges': len((original_set - simulation_set) - set(removed)),
        'original_file_sha256': sha256(FROZEN),
        'simulation_file_sha256': sha256(DERIVED),
        'original_canonical_edge_sha256': canonical_hash(original_edges),
        'simulation_canonical_edge_sha256': canonical_hash(reloaded),
        'self_loop_decision_file_sha256': sha256(SELF_DECISION),
        'biological_topology_modified': False,
        'network_rediscovery_performed': False,
    }
    write_json(DERIVATION, summary)
    out = Counter(tf for tf, target in retained)
    self_tfs = {tf for tf, target in removed}
    formal_rows = [
        {'target': tf,
         'self_loop_present_in_original': str(tf in self_tfs),
         'self_loop_removed_for_celloracle': str(tf in self_tfs),
         'nonself_targets_remaining': out[tf]}
        for tf in FORMAL_TFS]
    write_csv(FORMAL, formal_rows,
              ('target', 'self_loop_present_in_original',
               'self_loop_removed_for_celloracle', 'nonself_targets_remaining'))
    print('REPORT TOPOLOGY: original', len(original_edges), 'self-loops', len(removed),
          'simulation', len(reloaded), 'added', summary['added_edges'],
          'removed nonself', summary['removed_nonself_edges'])
    for row in formal_rows:
        print('REPORT FORMAL TF:', row)
    return summary


def freeze_environment():
    assert sys.version_info[:2] == (3, 10)
    names = ('celloracle', 'scanpy', 'anndata', 'numpy', 'pandas', 'scipy',
             'scikit-learn', 'velocyto', 'h5py')
    version_lines = [f'{name}=={importlib.metadata.version(name)}' for name in names]
    pip = subprocess.run([sys.executable, '-m', 'pip', 'freeze'], text=True,
                         capture_output=True, check=True).stdout
    conda = subprocess.run(['/home/LOCAL_USER/miniforge3/bin/conda', 'list', '-n', 'ko33_wsl'],
                           text=True, capture_output=True, check=True).stdout
    content = ('STEP33B1-R ko33_wsl runtime freeze\n'
               f'Python: {sys.version.split()[0]}\n'
               + '\n'.join(version_lines) + '\n\nPIP FREEZE\n' + pip
               + '\nCONDA LIST\n' + conda)
    (AUDIT / 'Step33B1_R_EnvironmentFreeze.txt').write_text(content, encoding='utf-8')


def load_frozen_adata():
    import h5py
    import numpy as np
    import pandas as pd
    import scanpy as sc

    meta = pd.read_csv(STEP29_META, dtype={'Cell_ID': str, 'patient': str,
                                            'final_malignant_state_id': str})
    genes = pd.read_csv(STEP29_GENES)['Gene'].astype(str).tolist()
    ids = meta['Cell_ID'].tolist()
    assert len(ids) == EXPECTED_CELLS and len(set(ids)) == EXPECTED_CELLS
    assert len(genes) == EXPECTED_GENES and len(set(genes)) == EXPECTED_GENES
    adata = sc.read_h5ad(STEP23_RAW)
    assert adata.n_obs == EXPECTED_CELLS
    assert set(adata.obs_names) == set(ids)
    assert set(genes).issubset(set(adata.var_names))
    adata = adata[ids, genes].copy()
    assert list(adata.obs_names) == ids and list(adata.var_names) == genes
    adata.obs['patient'] = pd.Categorical(meta['patient'].to_numpy())
    adata.obs['final_malignant_state_id'] = pd.Categorical(
        meta['final_malignant_state_id'].to_numpy())
    with h5py.File(STEP26, 'r') as handle:
        names = [x.decode() if isinstance(x, bytes) else str(x)
                 for x in handle['obs']['_index'][...]]
        embedding = np.asarray(handle['obsm']['X_umap_lineage_clean'][...])
    assert embedding.shape == (EXPECTED_CELLS, 2)
    assert len(set(names)) == EXPECTED_CELLS and set(names) == set(ids)
    position = {name: index for index, name in enumerate(names)}
    adata.obsm['X_umap_lineage_clean'] = embedding[[position[x] for x in ids], :]
    assert np.isfinite(adata.obsm['X_umap_lineage_clean']).all()
    assert adata.X.min() >= 0
    write_csv(INPUT_AUDIT, [
        {'item': 'Step23 raw source', 'path': str(STEP23_RAW), 'observed': adata.shape, 'status': 'PASS'},
        {'item': 'frozen cells in metadata order', 'path': str(STEP29_META), 'observed': len(ids), 'status': 'PASS'},
        {'item': 'frozen Step29 genes in frozen order', 'path': str(STEP29_GENES), 'observed': len(genes), 'status': 'PASS'},
        {'item': 'patient labels', 'path': str(STEP29_META), 'observed': adata.obs['patient'].nunique(), 'status': 'PASS'},
        {'item': 'MS3 cells', 'path': str(STEP29_META), 'observed': int((adata.obs['final_malignant_state_id'] == 'MS3').sum()), 'status': 'PASS'},
        {'item': 'frozen display embedding', 'path': str(STEP26), 'observed': 'X_umap_lineage_clean 12057x2', 'status': 'PASS'},
    ], ('item', 'path', 'observed', 'status'))
    return adata


def fit_model():
    import numpy as np
    import celloracle as co

    topology = json.loads(DERIVATION.read_text(encoding='utf-8'))
    assert topology['original_edges'] == 32530
    assert topology['simulation_edges'] == 32390
    assert topology['added_edges'] == topology['removed_nonself_edges'] == 0
    assert sha256(FROZEN) == topology['original_file_sha256']
    assert sha256(DERIVED) == topology['simulation_file_sha256']
    freeze_environment()
    fit_audit = {'step': 'Step33B1-R', 'status': 'STARTED',
                 'alpha': ALPHA, 'GRN_unit': 'cluster',
                 'cluster_column': 'final_malignant_state_id',
                 'shared_candidate_topology_edges': 32390,
                 'source_raw': str(STEP23_RAW), 'n_cells': EXPECTED_CELLS,
                 'n_genes': EXPECTED_GENES, 'PCA_evaluation_components': PCA_EVALUATION_COMPONENTS,
                 'k': KNN_K, 'balanced': True, 'b_sight': KNN_B_SIGHT,
                 'b_maxl': KNN_B_MAXL, 'n_jobs': 4}
    write_json(FIT_AUDIT, fit_audit)
    adata = load_frozen_adata()
    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(
        adata=adata, cluster_column_name='final_malignant_state_id',
        embedding_name='X_umap_lineage_clean', transform='natural_log')
    assert oracle.adata.shape == (EXPECTED_CELLS, EXPECTED_GENES)
    tfdict = defaultdict(list)
    for edge in read_csv(DERIVED):
        tfdict[edge['target_gene']].append(edge['tf'])
    assert sum(map(len, tfdict.values())) == 32390
    assert all(target not in tfs for target, tfs in tfdict.items())
    oracle.import_TF_data(TFdict=dict(tfdict))
    assert sum(map(len, oracle.TFdict.values())) == 32390
    print('REPORT IMPORT: raw 12057x12916; imported TFdict edges 32390', flush=True)
    oracle.perform_PCA(n_components=PCA_EVALUATION_COMPONENTS)
    variance = np.cumsum(oracle.pca.explained_variance_ratio_)
    # Exact CellOracle documented rule; 100 computed PCs cover its plotted range.
    pc_candidates = np.where(np.diff(np.diff(variance) > 0.002))[0]
    if len(pc_candidates) == 0 or int(pc_candidates[0]) < 1:
        fit_audit.update(status='STOP_DOCUMENTED_PC_RULE_FAILED',
                         PC_rule_candidates=pc_candidates.tolist())
        write_json(FIT_AUDIT, fit_audit)
        raise RuntimeError('CellOracle documented PC-selection rule failed; no alternative rule used')
    n_comps = min(int(pc_candidates[0]), 50)
    fit_audit.update(status='PCA_COMPLETE', n_comps=n_comps,
                     CellOracle_documented_PC_rule_index=int(pc_candidates[0]))
    write_json(FIT_AUDIT, fit_audit)
    print('REPORT PCA: documented rule index', int(pc_candidates[0]),
          'n_comps', n_comps, flush=True)
    oracle.knn_imputation(n_pca_dims=n_comps, k=KNN_K, balanced=True,
                          b_sight=KNN_B_SIGHT, b_maxl=KNN_B_MAXL, n_jobs=4)
    assert oracle.adata.layers['imputed_count'].shape == (EXPECTED_CELLS, EXPECTED_GENES)
    fit_audit['status'] = 'IMPUTATION_COMPLETE'
    write_json(FIT_AUDIT, fit_audit)
    write_json(SCORING, {
        'step': 'Step33B1-R', 'name': 'Baseline-Anchored CellOracle Perturbation Score',
        'status': 'LOCKED_NOT_YET_SCORED', 'baseline_layer': 'imputed_count',
        'perturbed_layer': 'simulated_count', 'reference_cells': EXPECTED_CELLS,
        'baseline_mean_and_SD_reused_for_KO': True, 'SD_ddof': 0,
        'primary_Core50': 'Core50-LTO-FOXM1', 'Core50_primary_gene_count': 49,
        'alpha': ALPHA, 'n_propagation': 3, 'GRN_unit': 'cluster',
        'n_comps': n_comps, 'k': KNN_K, 'balanced': True,
        'b_sight': KNN_B_SIGHT, 'b_maxl': KNN_B_MAXL,
        'source_step31_step32_absolute_score_reused': False})
    oracle.to_hdf5(str(PREFIT))
    fit_audit['prefit_oracle_checkpoint'] = str(PREFIT)
    fit_audit['status'] = 'PREFIT_SAVED'
    write_json(FIT_AUDIT, fit_audit)
    print('REPORT PREFIT: saved', PREFIT, flush=True)
    oracle.fit_GRN_for_simulation(GRN_unit='cluster', alpha=ALPHA,
                                  use_cluster_specific_TFdict=False,
                                  verbose_level=1)
    states = sorted(oracle.coef_matrix_per_cluster)
    expected_states = sorted(adata.obs['final_malignant_state_id'].cat.categories.tolist())
    assert states == expected_states and len(states) == 8
    assert all(x.shape == (EXPECTED_GENES, EXPECTED_GENES)
               for x in oracle.coef_matrix_per_cluster.values())
    ms3_coef = oracle.coef_matrix_per_cluster['MS3']
    foxm1_targets = [e['target_gene'] for e in read_csv(DERIVED) if e['tf'] == 'FOXM1']
    assert len(foxm1_targets) == 58
    foxm1_values = ms3_coef.loc['FOXM1', foxm1_targets]
    nonzero = int(np.count_nonzero(foxm1_values.to_numpy()))
    fit_audit.update(status='FITTED', fitted_states=states,
                     fitted_target_genes=len(tfdict),
                     failed_genes=[],
                     coefficient_matrix_dimensions=[EXPECTED_GENES, EXPECTED_GENES],
                     FOXM1_coefficient_available_in_MS3=True,
                     FOXM1_nonzero_outgoing_coefficients_in_MS3=nonzero,
                     FOXM1_frozen_nonself_targets=58)
    write_json(FIT_AUDIT, fit_audit)
    oracle.to_hdf5(str(FITTED))
    fit_audit.update(status='FITTED_SAVED', fitted_oracle_object=str(FITTED),
                     fitted_oracle_file_sha256=sha256(FITTED))
    write_json(FIT_AUDIT, fit_audit)
    print('REPORT FIT: states', states, 'target genes', len(tfdict),
          'FOXM1 MS3 nonzero outgoing', nonzero, 'saved', FITTED, flush=True)


def frozen_endpoints():
    core = [r['Gene'] for r in read_csv(CORE50)]
    assert len(core) == 50 and len(set(core)) == 50 and core.count('FOXM1') == 1
    mappings = read_csv(PROGRAMS)
    specs = [
        ('Core50-LTO-FOXM1', [g for g in core if g != 'FOXM1'], 'primary', 49),
        ('G2M_Independent', None, 'primary', 22),
        ('Mitotic_Entry_Independent', None, 'primary', 16),
        ('Spindle_Chromosome_Independent', None, 'primary', 29),
        ('Hepatic_Differentiation', None, 'secondary', 19),
        ('Bile_Acid_Metabolism', None, 'secondary', 18),
    ]
    result = {}
    for name, genes, role, expected in specs:
        if genes is None:
            genes = [r['requested_symbol'] for r in mappings
                     if r['program'] == name and r['mapping_status'] == 'retained']
        assert len(genes) == expected and len(set(genes)) == expected, name
        result[name] = {'genes': genes, 'role': role}
    return result


def audit_endpoint_coverage(oracle):
    import numpy as np

    layer = oracle.adata.layers['imputed_count']
    assert layer.shape == (EXPECTED_CELLS, EXPECTED_GENES)
    endpoints = frozen_endpoints()
    gene_index = {g: i for i, g in enumerate(oracle.adata.var_names)}
    cache = {}
    rows = []
    failed_primary = []
    for endpoint, spec in endpoints.items():
        for gene in spec['genes']:
            if gene in gene_index:
                col = np.asarray(layer[:, gene_index[gene]]).reshape(-1)
                finite = bool(np.isfinite(col).all())
                mean = float(np.mean(col)) if finite else None
                sd = float(np.std(col, ddof=0)) if finite else None
                okay = finite and sd > 0
            else:
                finite, mean, sd, okay = False, None, None, False
            rows.append({'endpoint': endpoint, 'role': spec['role'], 'gene': gene,
                         'gene_present': gene in gene_index,
                         'baseline_finite': finite, 'baseline_mean': mean,
                         'baseline_sd_ddof0': sd,
                         'status': 'PASS' if okay else 'FAILED'})
            if okay:
                cache[gene] = (gene_index[gene], mean, sd)
            elif spec['role'] == 'primary':
                failed_primary.append((endpoint, gene))
    write_csv(ENDPOINT_AUDIT, rows,
              ('endpoint', 'role', 'gene', 'gene_present', 'baseline_finite',
               'baseline_mean', 'baseline_sd_ddof0', 'status'))
    if failed_primary:
        raise RuntimeError('PERTURBATION_ENDPOINT_COVERAGE_FAILED: ' + repr(failed_primary))
    if any(r['status'] != 'PASS' for r in rows):
        raise RuntimeError('Secondary endpoint coverage failed; do not silently drop genes')
    return endpoints, cache


def simulate_foxm1():
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import celloracle as co

    fit = json.loads(FIT_AUDIT.read_text(encoding='utf-8'))
    assert fit['status'] == 'FITTED_SAVED'
    assert FITTED.is_file() and sha256(FITTED) == fit['fitted_oracle_file_sha256']
    oracle = co.load_hdf5(str(FITTED))
    assert oracle.adata.shape == (EXPECTED_CELLS, EXPECTED_GENES)
    assert oracle.GRN_unit == 'cluster' and oracle.alpha_for_trajectory_GRN == ALPHA
    endpoints, gene_stats = audit_endpoint_coverage(oracle)
    write_json(SIM_QC, {'step': 'Step33B1-R', 'status': 'STARTED',
                        'target': 'FOXM1', 'perturb_condition': {'FOXM1': 0.0},
                        'n_propagation': 3, 'GRN_unit': 'cluster',
                        'use_randomized_GRN': False, 'clip_delta_X': False})
    with warnings.catch_warnings(record=True) as observed:
        warnings.simplefilter('always')
        oracle.simulate_shift(perturb_condition={'FOXM1': 0.0},
                              GRN_unit='cluster', n_propagation=3,
                              use_randomized_GRN=False, clip_delta_X=False)
    warnings_text = [str(w.message) for w in observed]
    imputed = oracle.adata.layers['imputed_count']
    simulated = oracle.adata.layers['simulated_count']
    assert imputed.shape == simulated.shape == (EXPECTED_CELLS, EXPECTED_GENES)
    foxm1_col = list(oracle.adata.var_names).index('FOXM1')
    base_foxm1 = np.asarray(imputed[:, foxm1_col]).reshape(-1)
    sim_foxm1 = np.asarray(simulated[:, foxm1_col]).reshape(-1)
    foxm1_check = {
        'baseline_median_all_cells': float(np.median(base_foxm1)),
        'simulated_median_all_cells': float(np.median(sim_foxm1)),
        'simulated_max_all_cells': float(np.max(sim_foxm1)),
        'fraction_cells_suppressed': float(np.mean(sim_foxm1 < base_foxm1)),
        'simulated_nonzero_cells': int(np.count_nonzero(sim_foxm1)),
    }
    ood = oracle.evaluate_simulated_gene_distribution_range()
    ood.to_csv(AUDIT / 'Step33B1_R_OODGeneStats.csv')
    oracle.evaluate_and_plot_simulation_value_distribution(n_genes=1, n_bins=20,
                                                            alpha=0.5, save=None)
    fig = plt.gcf()
    fig.savefig(FIGURE_PNG, dpi=200, bbox_inches='tight')
    fig.savefig(FIGURE_PDF, bbox_inches='tight')
    plt.close('all')
    from celloracle.visualizations.config import CONFIG
    threshold = float(CONFIG['OOD_WARNING_EXCEEDING_PERCENTAGE']) / 100
    warning_genes = int((ood['Max_exceeding_ratio'] > threshold).sum())
    qc_review = warning_genes > 0 or any('out of distribution' in x.lower() for x in warnings_text)
    qc = {'step': 'Step33B1-R', 'status': 'FOXM1_SIMULATION_QC_REVIEW_REQUIRED' if qc_review else 'QC_PASS',
          'CellOracle_warnings': warnings_text,
          'CellOracle_OOD_warning_threshold_ratio': threshold,
          'genes_above_CellOracle_OOD_warning_threshold': warning_genes,
          'max_gene_OOD_exceeding_ratio': float(ood['Max_exceeding_ratio'].max()),
          'max_gene_OOD_cell_ratio': float(ood['OOD_cell_ratio'].max()),
          'FOXM1_direct_manipulation_check': foxm1_check,
          'distribution_diagnostic_png': str(FIGURE_PNG),
          'distribution_diagnostic_pdf': str(FIGURE_PDF),
          'OOD_gene_table': str(AUDIT / 'Step33B1_R_OODGeneStats.csv')}
    write_json(SIM_QC, qc)
    print('REPORT QC:', qc['status'], 'warning genes', warning_genes,
          'FOXM1 simulated nonzero cells', foxm1_check['simulated_nonzero_cells'], flush=True)
    if qc_review:
        print('DECISION: FOXM1_SIMULATION_QC_REVIEW_REQUIRED', flush=True)
        return
    obs = oracle.adata.obs
    ms3_mask = np.asarray(obs['final_malignant_state_id'] == 'MS3')
    assert int(ms3_mask.sum()) == 720
    patient = obs['patient'].astype(str).to_numpy()[ms3_mask]
    cell_ids = np.asarray(oracle.adata.obs_names)[ms3_mask]
    patient_counts = Counter(patient)
    delta_rows, patient_rows, endpoint_rows = [], [], []
    for endpoint, spec in endpoints.items():
        columns = [gene_stats[g][0] for g in spec['genes']]
        means = np.array([gene_stats[g][1] for g in spec['genes']])
        sds = np.array([gene_stats[g][2] for g in spec['genes']])
        baseline_score = ((np.asarray(imputed[ms3_mask][:, columns]) - means) / sds).mean(axis=1)
        knockout_score = ((np.asarray(simulated[ms3_mask][:, columns]) - means) / sds).mean(axis=1)
        delta = knockout_score - baseline_score
        assert len(delta) == 720 and np.isfinite(delta).all()
        for cell, pat, base, ko, change in zip(cell_ids, patient,
                                                baseline_score, knockout_score, delta):
            delta_rows.append({'Cell_ID': cell, 'patient': pat, 'state': 'MS3',
                               'endpoint': endpoint, 'role': spec['role'],
                               'n_genes': len(spec['genes']),
                               'ProgramScore_base': float(base),
                               'ProgramScore_KO': float(ko),
                               'DeltaProgram': float(change)})
        medians = {}
        for pat in sorted(patient_counts):
            pat_values = delta[patient == pat]
            pat_median = float(np.median(pat_values))
            medians[pat] = pat_median
            patient_rows.append({'patient': pat, 'n_MS3_cells': patient_counts[pat],
                                 'endpoint': endpoint, 'role': spec['role'],
                                 'median_DeltaProgram': pat_median,
                                 'Q1_DeltaProgram': float(np.quantile(pat_values, .25)),
                                 'Q3_DeltaProgram': float(np.quantile(pat_values, .75))})
        endpoint_rows.append({'endpoint': endpoint, 'role': spec['role'],
                              'n_genes': len(spec['genes']), 'n_MS3_cells': len(delta),
                              'n_MS3_patients': len(patient_counts),
                              'median_DeltaProgram': float(np.median(delta)),
                              'Q1_DeltaProgram': float(np.quantile(delta, .25)),
                              'Q3_DeltaProgram': float(np.quantile(delta, .75)),
                              'negative_patients': sum(x < 0 for x in medians.values()),
                              'positive_patients': sum(x > 0 for x in medians.values()),
                              'zero_patients': sum(x == 0 for x in medians.values()),
                              'patient_medians_json': json.dumps(medians, sort_keys=True)})
    write_csv(PROGRAM_DELTA, delta_rows,
              ('Cell_ID', 'patient', 'state', 'endpoint', 'role', 'n_genes',
               'ProgramScore_base', 'ProgramScore_KO', 'DeltaProgram'))
    write_csv(PATIENT_SUMMARY, patient_rows,
              ('patient', 'n_MS3_cells', 'endpoint', 'role', 'median_DeltaProgram',
               'Q1_DeltaProgram', 'Q3_DeltaProgram'))
    write_csv(ENDPOINT_SUMMARY, endpoint_rows,
              ('endpoint', 'role', 'n_genes', 'n_MS3_cells', 'n_MS3_patients',
               'median_DeltaProgram', 'Q1_DeltaProgram', 'Q3_DeltaProgram',
               'negative_patients', 'positive_patients', 'zero_patients',
               'patient_medians_json'))
    parameters = json.loads(SCORING.read_text(encoding='utf-8'))
    parameters['status'] = 'EXECUTED'
    parameters['MS3_patient_cell_counts'] = dict(patient_counts)
    write_json(SCORING, parameters)
    print('REPORT MS3 PATIENTS:', dict(patient_counts), flush=True)
    for row in endpoint_rows:
        print('REPORT ENDPOINT:', row, flush=True)
    print('DECISION: FOXM1_PILOT_PASS', flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--phase', choices=('topology', 'fit', 'simulate'), required=True)
    args = parser.parse_args()
    if args.phase == 'topology':
        derive_topology()
    elif args.phase == 'fit':
        fit_model()
    else:
        simulate_foxm1()


if __name__ == '__main__':
    main()
