# Reproducibility guide

## Static archive inspection

From the repository root with Python >=3.10:

```sh
python scripts/validate_release.py
python scripts/inspect_notebook.py --catalog
```

These are nonanalytical checks: hashes, JSON/CSV, notebook structure and Python syntax. They do not load matrices or deserialize models. A passing result is not an analytical reproduction claim.

## Future scientific reuse

Read README, docs/step1_35_route.csv, notebook_catalog.csv, manifests/analysis_dependency_map.csv and figure_table_provenance.csv. Select the intended final branch rather than Run All. Step10 lives in Step9; historical Step3/4 are not final M01; Step21/25 lack original producers; Step34 is the accepted native excerpt.

Retrieve matching public inputs through docs/data_access.md. Obtain approved derived objects through the planned separate archive (URLs pending); verify their historical SHA256 values against manifests/external_objects.csv. Raw third-party matrices should be reacquired from their source accessions, not supplied by a Git code clone.

Restore the recorded relative paths in an isolated workspace. files_manifest.csv maps each retained release object to its original scientific relative path; source_aliases.csv resolves removed duplicates. Keep original frozen inputs read-only and route new outputs to the writable restored copy. Set HCC_PROJECT_ROOT to that existing workspace and PYTHONPATH to this repository root. See README examples and docs/portability.md.

Restore context-specific historical environments separately. The Step33 conda export can be a starting point for that context only; installation/resolution was not tested during this audit. Incomplete version locks and generated-script/setup dependencies must be resolved transparently before execution.

Use frozen feature lists, coefficients, scaling, cutoff, state map, Core50/100, TF topology and drug queries. Do not reselect signatures or reconstruct missing original producers with newly written algorithms. Compare future outputs against saved sources only after a separately authorized reproduction. Preserve failed sensitivities and biological claim boundaries.

## What has and has not been verified

Release-copy prose repair, direct-path adaptation, scientific AST equivalence, source-table byte preservation, payload selection and static integrity were checked. No Step1–35 scientific analysis, full notebook execution, package installation or independent from-scratch reproduction was performed. Missing Step21/25 producer logic is blocking for the complete raw-to-final route; frozen-output reuse remains documented rather than claimed as full reproduction.
