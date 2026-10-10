# Historical environment evidence

The original analysis used several kernels, not one environment. See notebook_kernel_versions.csv and historical_software_evidence.csv (versions remain UNVERIFIED where only an import was saved).

- Bulk/TCGA: notebook Python 3.9.7 or 3.10.20; exact complete bulk package lock not recovered.
- CellRank Step28: root kernel metadata Python 3.12.10; separate WSL GPCCA notebooks must be considered separately.
- GRN/communication/independent validation/spatial: notebook Python 3.10.21 (WSL). Step29 saved environment JSON and Step32 full saved freeze are archived separately.
- CellOracle Step33: saved Linux/WSL conda export, Python 3.10.21 and CellOracle 0.18.0; environment.yml is this historical export only. Scanpy/anndata versions may differ from Step32. requirements.txt is only its pip subsection.
- Step34/35 used the ko33_celloracle_wsl kernel according to saved metadata, not an independently verified environment recreation.

No current pip freeze was substituted for the analytical environment. Dependency resolution/installation was not performed. Do not merge different environment pins to create a fictional universal environment.

The local archival runtime is recorded separately in Reports; it is not an analysis environment.


Additional saved environment JSON fields recovered: 12 context-specific software/version records. For example Step30 FINAL records Python3.10.21, LIANA1.7.3, CORNETO0.9.1a6, Scanpy1.11.5, AnnData0.11.4 and SciPy1.15.3. Import-only UNVERIFIED records remain evidence of imports rather than contradicting these scoped saved records.
