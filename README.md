# HCC-MS3 Transcriptomic Atlas: code and frozen-source release candidate

**Private repository; v1.0.0 candidate, not a public release.** This repository preserves the analysis sources and selected frozen evidence for the HCC-MS3 study. The [PRIVATE GitHub repository](https://github.com/Aqiuuu98/HCC-MS3-Transcriptomic-Atlas) has been created on `main`. Its initial commit (`d91af36c257de082757229c3262e9836387f9b0f`) uploaded 354 files; the remote file set and all 354 SHA256 values were verified. The repository is not publicly accessible. Author metadata, code/data licensing and public-release approval remain pending. No GitHub public release or Zenodo DOI exists.

## Study overview

The study connects TCGA-LIHC bulk transcriptomic and clinical analyses with malignant-cell programs in the GSE149614 discovery atlas. It includes patient-aware conservation, trajectory context, regulatory-network and communication inference; independent single-cell validation in GSE326201; spatial validation in GSE281759; virtual suppression of five frozen regulators; and frozen-query drug-induced transcriptional reversal using L1000FWD/LINCS.

This is an evidence-preserving archive, with **35 steps documented and 46 manuscript Figure/Table identifiers traced**. These counts describe documentation coverage. They do not establish successful independent end-to-end reproduction. No scientific analysis was rerun during release preparation.

## Repository structure

| Path | Contents |
|---|---|
| `notebooks/` | 39 selected source notebooks; original code retained apart from documented path adaptation and prose-cell repair; saved outputs removed |
| `scripts/archived/` | 16 preserved scientific scripts, including compatibility gates and publication helpers |
| `scripts/validate_release.py`, `scripts/inspect_notebook.py` | Static integrity and source-catalog inspection; no analysis execution |
| `hcc_release_paths.py` | Configured path adapter requiring an isolated restored analysis workspace |
| `environment/` | Context-specific historical environment evidence and partial dependency information |
| `frozen_objects/` | Lightweight gene sets, model CSV/JSON specifications, states and frozen result definitions |
| `source_data/` | Frozen plot inputs and all 2 main/25 supplementary Table DATA copies |
| `manifests/` | Candidate files, checksums, source aliases, dependencies, parameters and external-object plans |
| `docs/` | Study design, step route, access instructions, limitations and draft availability statements |

The local working directory also contains `Reports/`, `local_audit/`, nonselected historical copies and preparation scripts. They are excluded from Git. **Do not upload the entire working directory or its earlier draft ZIP.** The allowlist in `manifests/github_candidate_paths.json` defines the intended repository content. Duplicate sources are mapped in `manifests/source_aliases.csv`.

## Analysis route and source order

The complete step-to-notebook map is [step1_35_route.csv](docs/step1_35_route.csv), also available as [READ_ORDER.csv](notebooks/READ_ORDER.csv). The [notebook catalog](docs/notebook_catalog.csv) records 41 audited notebook records, of which 39 selected notebook files are tracked in Git. The excluded records are `notebooks/bulk/LIHC_Step15B_Stability_Validation.ipynb`, a source-cell duplicate of the retained `notebooks/bulk/Step15B_Stability_Validation.ipynb`, and the zero-cell `notebooks/trajectory/Untitled.ipynb`. These exclusions explain the 41-versus-39 counts. The catalog also records the original kernel metadata. These are dependency guides, not an executable Run-All driver.

| Steps | Scope and order |
|---|---|
| 1–2 | TCGA expression/cohort preparation and single-gene survival screening |
| 3–4 | Historical exploratory model/risk branches; not formal M01 producers; originals retained privately |
| 5 → 6/7 | Tumor-normal differential expression and functional interpretation |
| 8 → 9/10 → 11–14 | Original Top10/risk, clinical adjustment, immune context, WGCNA, frozen subtyping and pathways; Step10 code resides in Step9 |
| 15 → 16 → 17 | Model comparison/stability, conditional internal holdout, external-validation specification and bulk freeze |
| 18 → 19 → 20 → 21 → 22 → 23 | Discovery single-cell import, QC, representation, lineage annotation, malignancy support and lineage-clean states |
| 24 → 25 → 26 | Frozen bulk-anchor projection, conserved MS3 programs and biological context |
| 27 → 28 | Pseudotime and directionality/terminal-fate feasibility |
| 29 → 30 | Regulatory architectures and patient-aware communication hypotheses |
| 31 / 32 | Independent single-cell and spatial validation using frozen upstream programs |
| 33 → 34 | Five-TF virtual suppression and accepted native CellOracle displacement |
| 35 | Frozen queries, saved drug signatures, convergence and liver-context reversal |

**Step21 and Step25 original producer code has not been recovered.** Their saved final objects and definitions are documented. They can support downstream reuse after access approval, but do not replace the missing original generation logic. Raw-to-final regeneration of these steps is blocked by that code gap.

## Installation and static checks

Static inspection requires Python 3.10 or newer and its standard library:

```sh
python scripts/validate_release.py
python scripts/inspect_notebook.py --catalog
```

These commands parse JSON, Python syntax and checksums. They do not execute notebook cells, fit models, query APIs, deserialize joblib/CellOracle objects, or load expression matrices.

Scientific environments are separate from the static checker. See [environment/README.md](environment/README.md) and [historical software evidence](environment/historical_software_evidence.csv). `environment/environment.yml` and `environment/requirements.txt` are recorded **Step33 CellOracle** exports, not a universal Step1–35 specification. Their dependency resolution has not been retested. Some bulk, CellRank and regulatory-network version locks remain incomplete.

For a future authorized reproduction in the historical Linux/WSL context, the Step33 export may be restored separately:

```sh
conda env create -f environment/environment.yml
```

This installation was not performed during preparation. Do not install the latest versions of all imported packages as a substitute for the missing historical locks.

## Configuring a restored workspace

Obtain approved inputs, reproduce the original relative directory layout in an **isolated working copy**, and keep the original frozen project read-only. The source-to-repository mapping is in `manifests/files_manifest.csv`; excluded duplicates are listed in `manifests/source_aliases.csv`. Many heavy dependencies are intentionally absent from Git.

Start Jupyter from the cloned release, exposing its root on `PYTHONPATH`. Example for Linux/WSL:

```sh
export PYTHONPATH="$PWD${PYTHONPATH:+:$PYTHONPATH}"
export HCC_PROJECT_ROOT=/scratch/HCC_MS3_WORKSPACE
jupyter lab
```

PowerShell equivalent:

```powershell
$env:PYTHONPATH = (Get-Location).Path
$env:HCC_PROJECT_ROOT = 'D:\ResearchWorkspaces\HCC_MS3_WORKSPACE'
jupyter lab
```

The configured directory must already exist. The path adapter rejects the original frozen project and the release archive as execution workspaces. It does not download or restore inputs. [Portability notes](docs/portability.md) explain the remaining embedded-script/path and setup limitations. No automatic scientific execution is provided.

## Public data and large frozen objects

- TCGA-LIHC: obtain the matching STAR expression and clinical inputs through the [NCI GDC portal](https://portal.gdc.cancer.gov/projects/TCGA-LIHC), observing each file's open/controlled access designation. The saved Step1 audit records 424 mapped STAR files; this archive does not infer that count or redownload them.
- Discovery scRNA-seq: [GSE149614](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614).
- Independent scRNA-seq: [GSE326201](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE326201).
- Spatial validation: [GSE281759](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281759).
- Drug-induced signatures: [L1000FWD](https://maayanlab.cloud/L1000FWD/) and the documented frozen query/signature identities. Saved results are authoritative; requerying a live service is not guaranteed to reproduce its historical responses.

See [data_access.md](docs/data_access.md) for access boundaries and [external_objects.csv](manifests/external_objects.csv) for the 109 previously hashed external objects (17,148,067,078 bytes). This inventory includes intermediates and rejected-branch objects; it is **not** an automatic deposit list. [external_archiving_plan.csv](manifests/external_archiving_plan.csv) distinguishes proposed core derived objects, optional intermediates, acquisition-only matrices and excluded branches. Archive URLs and access terms remain pending. No raw TCGA/GEO download, h5ad/CellOracle object, full cache or draft ZIP is included in Git.

## Frozen model and biological definitions

Portable model CSV/JSON files preserve gene order, coefficients, imputation/scaling values and cutoffs. The M01 model was fitted on 238 development patients and evaluated on 102 internal holdout patients; its historical feature discovery was not entirely isolated from that holdout. The original RiskScore and M01 are distinct. EXT01 is a separate external-validation refit/specification, not the same M01 fitted object. Serialized joblib models are withheld pending privacy review; they were not deserialized.

Original Final Top10, MS1–MS8, Core50/Core100, patient-aware programs, five formal TFs (FOXM1, MYBL2, E2F8, TFDP1, E2F1), drug queries and final evidence tiers remain frozen. Step25 leave-one-patient-out results are discovery-cohort internal portability evidence. Step34 accepts native CellOracle transition geometry; the rejected expression-to-Harmony surrogate is not included as formal analysis evidence. Step35 query/context agreement is transcriptomic evidence within the same resource, not an independent drug-efficacy experiment.

## Figure/Table provenance and limits

[figure_table_provenance.csv](manifests/figure_table_provenance.csv) maps all 46 global identifiers (Figures1–8/S1–S11, Tables1–2/S1–S25) and 78 panel/table records. Plot sources, code locators, claim boundaries and explicitly withheld dependencies are recorded. All 27 final Table DATA copies retain their frozen bytes. Archived Step-local table/figure filenames are not substitutes for manuscript numbering.

Figure S7A has a confirmed artwork–caption mismatch: the final plot displays 39 candidates, comprising 28 primary robust (≥5% relative state-skew margin) and 11 secondary moderate (3–5%) candidates. The frozen legend describes only 28 and network edge widths. Saved tables and main/supplementary results support both tiers. No figure, frozen data or manuscript was changed here. An English replacement for panel A is proposed in [figure_scope_conflicts.md](docs/figure_scope_conflicts.md). The frozen Word legend has not been changed; author approval and manuscript synchronization remain necessary before a public manuscript-matched release.

Known reproducibility limits include missing Step21/25 producers, incomplete environment locks, external-object access, mixed historical notebook branches, and Step34 excerpt setup dependencies. Syntax/structure checks passed after non-scientific repairs, but no notebook was certified by independent Run All. Pseudotime does not establish lineage; inferred regulatory/communication networks do not establish causality; counterfactual displacement does not establish cell conversion; transcriptomic reversal does not establish therapeutic efficacy.

## Citation, licensing and availability

`CITATION.cff` is an explicitly incomplete draft awaiting confirmed authors/order/ORCIDs. `LICENSE` is a pending-rights notice, not a selected code license. Code, derived objects and third-party resources require separate rights decisions. [Draft Code Availability](docs/manuscript_code_availability.md) identifies the current PRIVATE repository. Manuscript wording and any archival identifiers require review after separate public-release approval.

The initial private commit and push are complete. This repository remains PRIVATE; no public GitHub release or Zenodo deposit/DOI has been created. Public release requires separate approval. Documentation synchronization does not certify independent scientific reproduction. See [known limitations](docs/known_limitations.md) before reuse.
