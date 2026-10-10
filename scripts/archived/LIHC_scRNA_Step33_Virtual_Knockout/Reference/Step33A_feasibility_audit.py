#!/usr/bin/env python3
"""Step33A read-only frozen-input, target, and method feasibility audit."""
from hcc_release_paths import project_path as _hcc_path

# Step: Step33A — Frozen-input / Target / Method Feasibility Audit
# Input: frozen Step25/26/27/29/31/32 files listed in INPUT_SPECS below
# Output: Step33A_FrozenTargetEligibility.csv,
#         Step33A_SecondaryVulnerabilityInventory.csv,
#         Step33A_FrozenInputInventory.csv,
#         Step33A_FrozenEndpointAvailability.csv,
#         Step33A_MethodFeasibility.csv,
#         Step33A_EnvironmentAudit.txt,
#         Step33A_FROZEN_DESIGN.json
# Conda environment: spatial32_wsl
# Python version: captured at runtime in Step33A_EnvironmentAudit.txt
# Frozen dependencies: existing spatial32_wsl; package versions captured at runtime
# No knockout, training, network inference, scoring, or package installation.

import csv
import hashlib
import importlib.metadata as metadata
import json
import os
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import h5py


TCGA = Path(_hcc_path(''))
ROOT = TCGA / 'LIHC_scRNA_Step33_Virtual_Knockout'
AUDIT = ROOT / 'Audit'
STEP25 = TCGA / 'LIHC_scRNA_Step25_CrossPatient_Robustness'
STEP26 = TCGA / 'LIHC_scRNA_Step26_MS3_Biology'
STEP27 = TCGA / 'LIHC_scRNA_Step27_Trajectory'
STEP29 = TCGA / 'LIHC_scRNA_Step29_GRN'
STEP31 = TCGA / 'LIHC_scRNA_Step31_Independent_Validation'
STEP32 = TCGA / 'LIHC_scRNA_Step32_Spatial_Validation'

PRIMARY = ('FOXM1', 'MYBL2', 'E2F8', 'TFDP1', 'E2F1', 'BRCA1')
SECONDARY = ('PLK1', 'AURKA', 'AURKB', 'TTK', 'BUB1', 'BUB1B', 'CENPE',
             'KIF2C', 'BIRC5', 'RRM2', 'FANCD2', 'UBE2T', 'CHEK1')
PRIMARY_ENDPOINTS = ('Core50', 'G2M_Independent', 'Mitotic_Entry_Independent',
                     'Spindle_Chromosome_Independent')
SECONDARY_ENDPOINTS = ('Hepatic_Differentiation', 'Bile_Acid_Metabolism')
ENDPOINT_COLUMNS = {
    'Core50': 'MS3_conserved_Core50_score',
    'G2M_Independent': 'Program__Proliferation_G2M',
    'Mitotic_Entry_Independent': 'DeepProgram__Mitotic_Entry',
    'Spindle_Chromosome_Independent': 'DeepProgram__Spindle_Chromosome_Segregation',
    'Hepatic_Differentiation': 'DeepProgram__Hepatic_Differentiation',
    'Bile_Acid_Metabolism': 'DeepProgram__Bile_Acid_Metabolism',
}

INPUT_SPECS = (
    ('frozen_Core50_membership', 'Step25', STEP25 / 'Reference/Step25_FROZEN_MS3_Conserved_Core50.csv', 'endpoint_source'),
    ('frozen_MS3_biology_object', 'Step26', STEP26 / 'Data/GSE149614_STEP26_FINAL_FROZEN_MS3_BIOLOGY.h5ad', 'provenance_only'),
    ('frozen_vulnerability_candidates', 'Step26', STEP26 / 'Reference/Step26_FROZEN_MS3_Mechanistic_Candidates.csv', 'provenance_only'),
    ('frozen_trajectory_object', 'Step27', STEP27 / 'Data/GSE149614_STEP27_FINAL_TRAJECTORY.h5ad', 'provenance_only'),
    ('frozen_raw_malignant_expression', 'Step29', STEP29 / 'Data/Step29D_GRN_INPUT_FILTERED_RAW_COUNTS.h5ad', 'primary_discovery_input'),
    ('frozen_log2_malignant_expression', 'Step29', STEP29 / 'Data/Step29E_GRN_INPUT_LOG2_COUNTS.h5ad', 'network_input_provenance'),
    ('frozen_cell_metadata', 'Step29', STEP29 / 'Tables/Step29D_Frozen_Cell_Metadata.csv', 'primary_discovery_metadata'),
    ('frozen_gene_detection_audit', 'Step29', STEP29 / 'Tables/Step29D_Gene_Filtering_Audit.csv', 'target_detection_source'),
    ('frozen_GRN_gene_list', 'Step29', STEP29 / 'Tables/Step29D_GRN_Genes.csv', 'network_gene_universe'),
    ('frozen_GRN_TF_list', 'Step29', STEP29 / 'Tables/Step29D_GRN_TFs.csv', 'network_TF_universe'),
    ('directed_GRN_adjacencies', 'Step29', STEP29 / 'Tables/Step29E_GRNBoost2_Adjacencies.tsv', 'network_edges_and_importance'),
    ('GRN_QC', 'Step29', STEP29 / 'Tables/Step29E_GRNBoost2_QC.csv', 'network_count_check'),
    ('motif_filtered_regulon_targets', 'Step29', STEP29 / 'Data/Step29F_Regulon_Targets.json', 'regulon_membership'),
    ('motif_filtered_regulon_object', 'Step29', STEP29 / 'Data/Step29F_Regulons.dat', 'regulon_provenance_only'),
    ('motif_filtered_regulon_summary', 'Step29', STEP29 / 'Tables/Step29F_Regulon_Summary.csv', 'regulon_target_counts'),
    ('AUCell_activity', 'Step29', STEP29 / 'Tables/Step29G_AUCell_RegulonActivity.csv', 'activity_availability'),
    ('frozen_MS3_regulator_shortlist', 'Step29', STEP29 / 'Tables/Step29H_MS3_CIN_Regulon_Shortlist.csv', 'target_provenance'),
    ('frozen_regulatory_architecture', 'Step29', STEP29 / 'Audit/Step29_FINAL_RegulatoryArchitecture.json', 'target_provenance'),
    ('external_Core50_mapping', 'Step31', STEP31 / 'Audit/Step31_Core50_FrozenValidationMapping.csv', 'validation_provenance_only'),
    ('external_program_mapping', 'Step31', STEP31 / 'Audit/Step31E2A_ValidationProgramMapping.csv', 'validation_provenance_only'),
    ('external_spatial_lock', 'Step32', STEP32 / 'Audit/Step32_FINAL_LOCK.json', 'validation_provenance_only'),
)

METHOD_SOURCES = {
    'frozen Step29 directed GRN/regulon based perturbation': '',
    'CellOracle': 'https://morris-lab.github.io/CellOracle.documentation/modules/celloracle.html',
    'scTenifoldKnk': 'https://github.com/cailab-tamu/scTenifoldKnk',
    'GEARS': 'https://github.com/snap-stanford/GEARS',
    'scGen': 'https://scgen.readthedocs.io/en/stable/',
    'dynamo': 'https://github.com/aristoteleo/dynamo-release',
}


def rows(path, delimiter=','):
    with path.open(newline='', encoding='utf-8') as stream:
        yield from csv.DictReader(stream, delimiter=delimiter)


def write_csv(path, records, columns):
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(records)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def package_version(names):
    for name in names:
        try:
            return metadata.version(name)
        except metadata.PackageNotFoundError:
            pass
    return None


def h5ad_metadata(path):
    with h5py.File(path, 'r') as data:
        if 'X' not in data:
            return '', '', '', ''
        x = data['X']
        shape = tuple(x.attrs['shape']) if isinstance(x, h5py.Group) else x.shape
        obs = [k for k in data['obs'].keys() if k != '_index'] if 'obs' in data else []
        uns = list(data['uns'].keys()) if 'uns' in data else []
        return shape[0], shape[1], ';'.join(obs), ';'.join(uns)


def inventory():
    output = []
    for role, step, path, use in INPUT_SPECS:
        record = {'input_role': role, 'source_step': step, 'path': str(path),
                  'exists': path.is_file(), 'size': path.stat().st_size if path.is_file() else '',
                  'object_type': path.suffix.lower().lstrip('.'), 'cells': '', 'genes': '',
                  'key_obs_fields': '', 'key_uns_fields': '',
                  'frozen_status': 'FROZEN_READ_ONLY' if path.is_file() else 'MISSING',
                  'used_for_step33': use}
        if path.is_file() and path.suffix.lower() == '.h5ad':
            record['cells'], record['genes'], record['key_obs_fields'], record['key_uns_fields'] = h5ad_metadata(path)
        elif role == 'frozen_cell_metadata' and path.is_file():
            record['key_obs_fields'] = ';'.join(next(rows(path)).keys())
        output.append(record)
    return output


def environment_audit():
    executable = Path(sys.executable)
    env_name = executable.parent.parent.name
    requested_packages = {
        'h5py': ('h5py',), 'anndata': ('anndata',), 'scanpy': ('scanpy',),
        'pandas': ('pandas',), 'numpy': ('numpy',),
        'CellOracle': ('celloracle',), 'scTenifoldKnk_R': (),
        'scTenifoldpy': ('scTenifoldpy',), 'GEARS': ('cell-gears', 'gears'),
        'scGen': ('scgen',), 'dynamo': ('dynamo-release', 'dynamo'),
    }
    versions = {name: package_version(candidates) if candidates else None
                for name, candidates in requested_packages.items()}
    rscript = shutil.which('Rscript')
    if rscript:
        try:
            probe = subprocess.run([rscript, '-e', 'cat(as.character(utils::packageVersion("scTenifoldKnk")))'],
                                   capture_output=True, text=True, timeout=20, check=False)
            if probe.returncode == 0:
                versions['scTenifoldKnk_R'] = probe.stdout.strip()
        except (OSError, subprocess.TimeoutExpired):
            pass
    return env_name, versions, rscript


def main():
    assert TCGA.is_dir(), 'LOCAL_FILESYSTEM_UNAVAILABLE'
    missing = [str(p) for _, _, p, _ in INPUT_SPECS if not p.is_file()]
    if missing:
        raise FileNotFoundError('Required frozen input missing: ' + '; '.join(missing))
    assert ROOT.is_dir() and AUDIT.is_dir(), 'Step33 output directories missing'

    inv = inventory()
    metadata_rows = list(rows(STEP29 / 'Tables/Step29D_Frozen_Cell_Metadata.csv'))
    n_cells = len(metadata_rows)
    patients = sorted({r['patient'] for r in metadata_rows})
    state_counts = Counter(r['final_malignant_state_id'] for r in metadata_rows)
    assert n_cells == 12057 and set(state_counts) == {f'MS{i}' for i in range(1, 9)}
    assert len({r['Cell_ID'] for r in metadata_rows}) == n_cells
    assert len(patients) == 5

    gene_audit = {r['Gene']: r for r in rows(STEP29 / 'Tables/Step29D_Gene_Filtering_Audit.csv')}
    grn_genes = {r['Gene'] for r in rows(STEP29 / 'Tables/Step29D_GRN_Genes.csv')}
    grn_tfs = {r['TF'] for r in rows(STEP29 / 'Tables/Step29D_GRN_TFs.csv')}
    with h5py.File(STEP29 / 'Data/Step29D_GRN_INPUT_FILTERED_RAW_COUNTS.h5ad', 'r') as raw_h5:
        matrix_genes = {g.decode('utf-8') if isinstance(g, bytes) else str(g)
                        for g in raw_h5['var']['_index'][...]}
    assert matrix_genes == grn_genes, 'Frozen H5AD genes differ from Step29 GRN gene list'
    focus = set(PRIMARY) | set(SECONDARY)
    outgoing = Counter()
    incoming = Counter()
    weights_valid = Counter()
    total_edges = 0
    for edge in rows(STEP29 / 'Tables/Step29E_GRNBoost2_Adjacencies.tsv', '\t'):
        total_edges += 1
        tf, target = edge['TF'], edge['target']
        if tf in focus:
            outgoing[tf] += 1
            try:
                weights_valid[tf] += float(edge['importance']) >= 0
            except ValueError:
                pass
        if target in focus:
            incoming[target] += 1
    qc = next(rows(STEP29 / 'Tables/Step29E_GRNBoost2_QC.csv'))
    assert total_edges == int(qc['N_edges'])
    assert len(grn_genes) == 12916 and len(grn_tfs) == 1329

    regulon_targets = json.loads((STEP29 / 'Data/Step29F_Regulon_Targets.json').read_text(encoding='utf-8'))
    regulon_summary = list(rows(STEP29 / 'Tables/Step29F_Regulon_Summary.csv'))
    regulons_by_tf = defaultdict(list)
    for r in regulon_summary:
        regulons_by_tf[r['TF']].append(r)
    with (STEP29 / 'Tables/Step29G_AUCell_RegulonActivity.csv').open(newline='', encoding='utf-8') as stream:
        aucell_header = next(csv.reader(stream))
        n_aucell_cells = sum(1 for _ in stream)
    assert n_aucell_cells == n_cells
    assert set(regulon_targets) == set(aucell_header[1:])

    def target_row(target, primary):
        ga = gene_audit.get(target)
        regs = regulons_by_tf[target]
        names = [r['Regulon'] for r in regs]
        aucell_available = bool(names) and all(name in aucell_header for name in names)
        expr_present = target in matrix_genes and ga is not None and int(ga['Detected_cells']) > 0
        edge_weights = outgoing[target] > 0 and weights_valid[target] == outgoing[target]
        if primary:
            if not expr_present or not outgoing[target]:
                eligibility = 'METHOD_INELIGIBLE'
                reason = 'No detected expression or no frozen outgoing GRN edges; cannot support target-specific frozen-network perturbation.'
            elif target == 'BRCA1':
                eligibility = 'METHOD_INELIGIBLE'
                reason = ('BRCA1 has a computational Step29 regulon label and GRN edges, but is a noncanonical DNA-repair regulator; TF-specific perturbation eligibility is not established. No approved gene-agnostic frozen-network method exists. Reassess only after an explicit method decision.')
            else:
                eligibility = 'METHOD_DECISION_REQUIRED'
                reason = ('Detected and represented in frozen directed GRN, motif-filtered regulon, and AUCell. GRNBoost2 importance is unsigned, not a perturbation-response coefficient; no approved frozen-network KO method is locked.')
        else:
            eligibility = 'INVENTORY_ONLY'
            reason = ('Secondary frozen vulnerability pool only. Expression and Step29 network representation are inventory facts; no KO was tested. An incoming target edge is not an outgoing perturbation route.')
        if primary and target != 'BRCA1':
            gene_specific_support = 'PENDING_METHOD_LOCK; frozen outgoing TF edges exist'
        elif primary:
            gene_specific_support = 'TF_SPECIFIC_INELIGIBLE; gene-agnostic method would need separate approval'
        else:
            gene_specific_support = ('NO_DIRECT_FROZEN_TF_OUTGOING_ROUTE; gene-agnostic method would need separate approval'
                                     if outgoing[target] == 0 else 'PENDING_METHOD_LOCK')
        result = {
            'target': target,
            'target_class': 'MS3_frozen_regulator_candidate' if primary and target != 'BRCA1' else ('MS3_frozen_noncanonical_regulator_candidate' if primary else 'secondary_frozen_vulnerability_pool'),
            'source_step': 'Step29_FINAL_RegulatoryArchitecture' if primary else 'User-specified Step33A frozen pool; Step26 context',
            'expression_present': expr_present,
            'expression_detection_fraction': ga['Detected_fraction'] if ga else '',
            'present_in_grn': target in grn_genes,
            'present_as_GRN_TF': target in grn_tfs,
            'n_outgoing_edges': outgoing[target],
            'n_incoming_edges': incoming[target],
            'present_as_regulon': bool(regs),
            'regulon_target_count': sum(int(r['N_targets']) for r in regs),
            'aucell_available': aucell_available,
            'network_weight_available': edge_weights,
            'gene_specific_perturbation_compatibility': gene_specific_support,
            'methodologically_perturbable': eligibility,
            'eligibility_reason': reason,
        }
        return result

    primary_rows = [target_row(t, True) for t in PRIMARY]
    secondary_rows = [target_row(t, False) for t in SECONDARY]
    assert all(r['expression_present'] for r in primary_rows)
    assert all(r['present_in_grn'] for r in primary_rows)
    assert all(r['n_outgoing_edges'] > 0 for r in primary_rows)
    assert all(r['n_outgoing_edges'] == 0 for r in secondary_rows)

    with h5py.File(STEP29 / 'Data/Step29D_GRN_INPUT_FILTERED_RAW_COUNTS.h5ad', 'r') as raw_h5:
        raw_layers = list(raw_h5['layers'].keys()) if 'layers' in raw_h5 else []
    with h5py.File(STEP27 / 'Data/GSE149614_STEP27_FINAL_TRAJECTORY.h5ad', 'r') as trajectory_h5:
        trajectory_layers = list(trajectory_h5['layers'].keys()) if 'layers' in trajectory_h5 else []
    assert not any('spliced' in layer.lower() or 'velocity' in layer.lower()
                   for layer in raw_layers + trajectory_layers)

    endpoint_rows = []
    for endpoint in PRIMARY_ENDPOINTS + SECONDARY_ENDPOINTS:
        column = ENDPOINT_COLUMNS[endpoint]
        nonmissing = sum(r.get(column, '') not in ('', 'nan', 'NaN') for r in metadata_rows)
        endpoint_rows.append({
            'endpoint': endpoint, 'role': 'PRIMARY' if endpoint in PRIMARY_ENDPOINTS else 'SECONDARY_DIRECTIONAL',
            'source_step': 'Step25/Step26 frozen program; Step29D cell metadata',
            'frozen_score_column': column, 'cells_with_score': nonmissing,
            'available': nonmissing == n_cells,
            'expected_direction_if_future_KO': 'decrease' if endpoint in PRIMARY_ENDPOINTS else 'increase_supportive_only',
            'source_path': str(STEP29 / 'Tables/Step29D_Frozen_Cell_Metadata.csv'),
        })
    assert all(r['available'] for r in endpoint_rows)

    env_name, versions, rscript = environment_audit()
    assert env_name == 'spatial32_wsl', f'Unexpected Python environment: {env_name}'
    installed = {
        'CellOracle': versions['CellOracle'] is not None,
        'scTenifoldKnk': versions['scTenifoldKnk_R'] is not None or versions['scTenifoldpy'] is not None,
        'GEARS': versions['GEARS'] is not None,
        'scGen': versions['scGen'] is not None,
        'dynamo': versions['dynamo'] is not None,
    }
    method_rows = [
        {'method': 'frozen Step29 directed GRN/regulon based perturbation', 'installed': 'N/A: frozen data, no implementation',
         'requires_new_environment': False, 'requires_perturbation_training_data': False,
         'requires_velocity': False, 'requires_new_GRN_reconstruction': False,
         'can_use_frozen_Step29_GRN': True, 'supports_gene_KO': 'UNDETERMINED: no approved implementation',
         'supports_non_TF_targets': 'UNDETERMINED',
         'main_limitations': 'Step29 GRNBoost2 importance is nonnegative/unsigned and not a causal response coefficient. Regulon membership/AUCell are descriptive. No published compatible frozen-network KO implementation or parameters have been locked; do not invent propagation/scaling.',
         'recommendation_status': 'METHOD_DECISION_REQUIRED'},
        {'method': 'CellOracle', 'installed': installed['CellOracle'],
         'requires_new_environment': not installed['CellOracle'], 'requires_perturbation_training_data': False,
         'requires_velocity': False, 'requires_new_GRN_reconstruction': True,
         'can_use_frozen_Step29_GRN': 'PRIOR_TF_TARGET_LIST_ONLY; not direct simulation coefficients',
         'supports_gene_KO': 'TF-focused', 'supports_non_TF_targets': 'NOT_ESTABLISHED',
         'main_limitations': 'Can import a supplied TF-target base GRN, but official simulation fits new GRN regression coefficients. Step29 importance weights cannot be substituted directly; new GRN inference violates frozen-network design. BRCA1 TF suitability unverified.',
         'recommendation_status': 'INCOMPATIBLE_WITH_FROZEN_DESIGN'},
        {'method': 'scTenifoldKnk', 'installed': installed['scTenifoldKnk'],
         'requires_new_environment': not installed['scTenifoldKnk'], 'requires_perturbation_training_data': False,
         'requires_velocity': False, 'requires_new_GRN_reconstruction': True,
         'can_use_frozen_Step29_GRN': False, 'supports_gene_KO': True,
         'supports_non_TF_targets': True,
         'main_limitations': 'Official workflow constructs its own denoised scGRN from wild-type scRNA-seq, then zeros target outdegree edges. Rebuilding a GRN conflicts with Step29 frozen network lock.',
         'recommendation_status': 'INCOMPATIBLE_WITH_FROZEN_DESIGN'},
        {'method': 'GEARS', 'installed': installed['GEARS'],
         'requires_new_environment': not installed['GEARS'], 'requires_perturbation_training_data': True,
         'requires_velocity': False, 'requires_new_GRN_reconstruction': False,
         'can_use_frozen_Step29_GRN': False, 'supports_gene_KO': True,
         'supports_non_TF_targets': True,
         'main_limitations': 'Needs measured perturbation-screen training examples; frozen observational GSE149614 and Step29 GRN alone are insufficient. External generic training would change the design and require a separate decision.',
         'recommendation_status': 'INCOMPATIBLE_WITH_FROZEN_DESIGN'},
        {'method': 'scGen', 'installed': installed['scGen'],
         'requires_new_environment': not installed['scGen'], 'requires_perturbation_training_data': True,
         'requires_velocity': False, 'requires_new_GRN_reconstruction': False,
         'can_use_frozen_Step29_GRN': False, 'supports_gene_KO': 'CONDITION_RESPONSE_MODEL; not gene-specific KO from GRN',
         'supports_non_TF_targets': 'NOT_APPLICABLE_WITHOUT_TRAINING_CONDITIONS',
         'main_limitations': 'Requires control/perturbed training conditions; no measured target-specific KO condition is frozen here. Does not directly consume Step29 GRN.',
         'recommendation_status': 'INCOMPATIBLE_WITH_FROZEN_DESIGN'},
        {'method': 'dynamo', 'installed': installed['dynamo'],
         'requires_new_environment': not installed['dynamo'], 'requires_perturbation_training_data': False,
         'requires_velocity': True, 'requires_new_GRN_reconstruction': False,
         'can_use_frozen_Step29_GRN': False, 'supports_gene_KO': 'GENE_PERTURBATION_WITH_FITTED_DYNAMICS',
         'supports_non_TF_targets': 'POTENTIALLY_WITH_VALID_DYNAMICS',
         'main_limitations': 'Perturbation inference depends on a learned velocity/vector field. Frozen GSE149614 inputs contain no spliced/unspliced or metabolic-labeling velocity layers; cannot force a velocity analysis.',
         'recommendation_status': 'INCOMPATIBLE_WITH_FROZEN_DESIGN'},
    ]
    for method in method_rows:
        method['primary_source_url'] = METHOD_SOURCES[method['method']]
    allowed = {'POTENTIALLY_COMPATIBLE', 'INCOMPATIBLE_WITH_FROZEN_DESIGN', 'METHOD_DECISION_REQUIRED'}
    assert all(m['recommendation_status'] in allowed for m in method_rows)

    input_cols = ('input_role', 'source_step', 'path', 'exists', 'size', 'object_type',
                  'cells', 'genes', 'key_obs_fields', 'key_uns_fields', 'frozen_status', 'used_for_step33')
    target_cols = ('target', 'target_class', 'source_step', 'expression_present',
                   'expression_detection_fraction', 'present_in_grn', 'present_as_GRN_TF',
                   'n_outgoing_edges', 'n_incoming_edges', 'present_as_regulon',
                   'regulon_target_count', 'aucell_available', 'network_weight_available',
                   'gene_specific_perturbation_compatibility',
                   'methodologically_perturbable', 'eligibility_reason')
    method_cols = ('method', 'installed', 'requires_new_environment',
                   'requires_perturbation_training_data', 'requires_velocity',
                   'requires_new_GRN_reconstruction', 'can_use_frozen_Step29_GRN',
                   'supports_gene_KO', 'supports_non_TF_targets', 'main_limitations',
                   'recommendation_status', 'primary_source_url')
    write_csv(AUDIT / 'Step33A_FrozenInputInventory.csv', inv, input_cols)
    write_csv(AUDIT / 'Step33A_FrozenTargetEligibility.csv', primary_rows, target_cols)
    write_csv(AUDIT / 'Step33A_SecondaryVulnerabilityInventory.csv', secondary_rows, target_cols)
    write_csv(AUDIT / 'Step33A_FrozenEndpointAvailability.csv', endpoint_rows, tuple(endpoint_rows[0]))
    write_csv(AUDIT / 'Step33A_MethodFeasibility.csv', method_rows, method_cols)

    source_paths = [p for _, _, p, _ in INPUT_SPECS]
    environment_lines = [
        'STEP33A ENVIRONMENT AUDIT (read-only method check)',
        f'Python version: {sys.version.split()[0]}', f'Python executable: {sys.executable}',
        f'Active conda environment (interpreter prefix): {env_name}',
        f'CONDA_DEFAULT_ENV process variable: {os.environ.get("CONDA_DEFAULT_ENV", "<not set; interpreter selected by absolute path>")}',
        f'Rscript available: {bool(rscript)}', f'Rscript path: {rscript or "<none>"}',
        'Package versions in the spatial32_wsl interpreter:',
    ]
    environment_lines += [f'  {name}: {value or "NOT_INSTALLED"}' for name, value in versions.items()]
    environment_lines += ['', 'Exact frozen input paths (read-only):']
    environment_lines += [f'  {p}' for p in source_paths]
    environment_lines += ['', 'No package installed. No new environment created. No model trained.']
    (AUDIT / 'Step33A_EnvironmentAudit.txt').write_text('\n'.join(environment_lines) + '\n', encoding='utf-8')

    result = {
        'step': 'Step33A', 'status': 'COMPLETE_METHOD_DECISION_REQUIRED',
        'discovery_dataset': 'GSE149614',
        'primary_targets': list(PRIMARY), 'secondary_vulnerability_pool': list(SECONDARY),
        'primary_endpoints': list(PRIMARY_ENDPOINTS),
        'secondary_directional_endpoints': list(SECONDARY_ENDPOINTS),
        'expected_future_direction': {**{k: 'down' for k in PRIMARY_ENDPOINTS},
                                      **{k: 'up_supportive_not_required' for k in SECONDARY_ENDPOINTS}},
        'original_RiskScore_reconstructed': False,
        'target_reselection_allowed': False, 'network_rediscovery_allowed': False,
        'formal_knockout_performed': False, 'Step34_started': False, 'Step35_started': False,
        'discovery_universe': {'cells': n_cells, 'genes_in_Step29_GRN_input': len(grn_genes),
                               'patients': patients, 'n_patients': len(patients),
                               'frozen_states': dict(sorted(state_counts.items())),
                               'MS3_cells': state_counts['MS3'], 'MS6_cells': state_counts['MS6'],
                               'metadata_path': str(STEP29 / 'Tables/Step29D_Frozen_Cell_Metadata.csv'),
                               'expression_path': str(STEP29 / 'Data/Step29D_GRN_INPUT_FILTERED_RAW_COUNTS.h5ad')},
        'network': {'GRN_edges': total_edges, 'TFs_with_edges': int(qc['N_TFs_with_edges']),
                    'genes': len(grn_genes), 'motif_filtered_regulons': len(regulon_targets),
                    'AUCell_cells': n_aucell_cells, 'AUCell_regulons': len(aucell_header) - 1,
                    'importance_weights_available': True,
                    'weights_are_signed_perturbation_coefficients': False},
        'method_decision_required': True,
        'method_installation_scope': 'spatial32_wsl interpreter',
        'compatible_method_locked': False,
        'velocity_input_audit': {'Step29_raw_h5ad_layers': raw_layers,
                                 'Step27_final_h5ad_layers': trajectory_layers,
                                 'spliced_unspliced_velocity_available': False},
        'BRCA1_caveat': 'Computational regulon and edges do not establish suitability for direct TF-specific knockout simulation.',
        'external_validation_firewall': 'Step31 and Step32 are provenance only; no training or network inference from GSE326201 or GSE281759.',
        'primary_method_sources': METHOD_SOURCES,
        'input_inventory_path': str(AUDIT / 'Step33A_FrozenInputInventory.csv'),
        'target_eligibility_path': str(AUDIT / 'Step33A_FrozenTargetEligibility.csv'),
        'method_feasibility_path': str(AUDIT / 'Step33A_MethodFeasibility.csv'),
        'environment_audit_path': str(AUDIT / 'Step33A_EnvironmentAudit.txt'),
        'script_path': str(Path(__file__)),
    }
    (AUDIT / 'Step33A_FROZEN_DESIGN.json').write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')

    print('STEP33A — VIRTUAL KNOCKOUT FEASIBILITY AUDIT REPORT')
    print(f'[1] Environment: {env_name}; Python {sys.version.split()[0]}; Rscript={bool(rscript)}; method packages installed={sum(installed.values())}/5')
    print(f'[2] Frozen discovery universe: {n_cells} cells; {len(grn_genes)} GRN genes; {len(patients)} patients; MS3={state_counts["MS3"]}; MS6={state_counts["MS6"]}')
    print('[3] Primary target eligibility:')
    for r in primary_rows:
        print(f'  {r["target"]}: detection={float(r["expression_detection_fraction"]):.4f}; outgoing={r["n_outgoing_edges"]}; regulon_targets={r["regulon_target_count"]}; {r["methodologically_perturbable"]}')
    print('[4] Secondary vulnerability inventory: ' + ', '.join(SECONDARY) + ' (inventory only; no KO)')
    print('[5] Frozen endpoint availability: ' + ', '.join(f'{r["endpoint"]}={r["cells_with_score"]}/{n_cells}' for r in endpoint_rows))
    print(f'[6] Step29 network: {total_edges} GRN edges; {len(regulon_targets)} motif-filtered regulons; AUCell {n_aucell_cells} cells x {len(aucell_header)-1} regulons; unsigned importance weights available')
    print('[7] Method feasibility:')
    for m in method_rows:
        print(f'  {m["method"]}: {m["recommendation_status"]}')
    print('[8] Compatible options: none locked; frozen-network method requires an explicit published method/parameter decision. CellOracle/scTenifoldKnk need new GRN fitting; GEARS/scGen need perturbation training; dynamo needs velocity.')
    print('[9] Firewall:')
    print('Formal knockout performed: False')
    print('New GRN constructed: False')
    print('Targets reselected: False')
    print('Original RiskScore reconstructed: False')
    print('Step34 started: False')
    print('Step35 started: False')
    print('METHOD_DECISION_REQUIRED')
    print('STEP33A: COMPLETE')
    print('NEXT: return feasibility report to ChatGPT for method lock.')


if __name__ == '__main__':
    main()
