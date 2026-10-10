# Historical environments and dependency scope

The static release checker uses only Python >=3.10 standard library. It does not require the scientific environments. No environment was installed or overwritten during final preparation.

environment.yml is the saved Step33 clean CellOracle conda export for Linux/WSL; requirements.txt is its recorded pip subsection. They are not a universal Step1–35 lock. The historical folder preserves additional Step29/32/33 context-specific exports. Do not infer unrecorded versions from the preparation runtime.

historical_software_evidence.csv records the source and scope of each saved package/version. import_inventory.csv lists observed imports, not a tested dependency resolver or exact package-name mapping. notebook_kernel_versions.csv records saved notebook metadata, not successful execution.

Examples of scoped evidence: Step30 final environment records Python3.10.21, LIANA1.7.3, CORNETO0.9.1a6, Scanpy1.11.5, AnnData0.11.4 and SciPy1.15.3; Step33 compatibility/pilot records CellOracle0.18.0. Those versions apply to their recorded contexts. Bulk/CellRank/full GRN locks remain incomplete.

Review the preserved script headers and step dependency route before restoring each environment. Third-party motif, gene-set, prior-network, reference annotation and L1000 resources require their own acquisition instructions and terms. Missing libraries, solver backends and frozen inputs are execution prerequisites; no latest-version substitutions or dependency-resolution success are claimed here.
