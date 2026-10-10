# Data access and external dependencies

| Resource | Access route | Release treatment |
|---|---|---|
| TCGA-LIHC | https://portal.gdc.cancer.gov/projects/TCGA-LIHC | Retrieve the original matching STAR and clinical files; respect per-file access designation. No raw download or clinical row bundle in Git. |
| GSE149614 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614 | Discovery scRNA inputs; acquire from GEO, retain submitted sample labels. |
| GSE326201 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE326201 | Independent scRNA inputs; original 18 filtered 10x matrices are acquisition-only. |
| GSE281759 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281759 | Spatial inputs; derived section-QC objects are separate proposed archive material. |
| L1000FWD/LINCS | https://maayanlab.cloud/L1000FWD/ | Retain frozen query/signature identities and derived result tables. Raw service response/metadata redistribution requires rights review. |
| Motif, pathway, prior-network and gene references | Refer to preserved Step29/30/31/32 code, import and environment evidence | Reacquire original matching resources with version and license evidence. Missing detailed pins are documented; no invented retrieval version. |

The external object inventory includes109 previously hashed files totaling17,148,067,078 bytes. See external_objects.csv and external_archiving_plan.csv. SHA256 values are reused from the completed stable-read audit; final preparation checked availability/size without repeating the large-object hash or loading matrices. The planned derived-data repository is not yet created; URLs/DOIs remain pending.

GDC open and controlled data have distinct access processes; users must not attempt participant reidentification. GEO makes submitted data accessible, but third-party rights may still apply. L1000FWD states academic non-profit use terms for its tools/services, which do not establish a blanket redistribution grant for all saved signatures. The original frozen clinical audit and unsafe serialization content remain local pending review.

Policy sources checked during preparation: [GDC access](https://gdc.cancer.gov/access-data/data-access-processes-and-tools), [GDC analysis policy](https://gdc.cancer.gov/analyze-data/data-analysis-policies), [GEO disclaimer](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html), [L1000FWD terms](https://github.com/MaayanLab/L1000FWD). Source endpoints and service versions may change; saved frozen evidence remains distinct from a new query.
