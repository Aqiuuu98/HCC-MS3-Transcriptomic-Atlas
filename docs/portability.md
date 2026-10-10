# Publication-copy portability and repairs

Only release copies were edited. Original notebook/script files, frozen numerical data, figures and manuscripts were not rewritten. The private baseline and per-file unified diffs are retained under local_audit/FINAL_PREPARATION_BASELINE_20261010 and Reports/FINAL_PREPARATION_DIFFS (excluded from Git).

- 53 code cells contained only Markdown prose and saved SyntaxError outputs. They were converted to Markdown with the source text unchanged; saved errors remain in the private baseline.
- 152 literal project-path expressions across notebook/script sources were adapted to hcc_release_paths.project_path. Adjacent Python string literals were handled as a single AST constant. After normalizing only path constants and adapter imports, each changed executable cell/script has an identical scientific AST.
- All notebook execution outputs and counts were removed from the Git copies. Saved outputs remain available in the private source snapshots. Cell IDs were added where absent; widget output metadata were removed.
- Random seeds, modeling parameters, gene membership/order, feature selection, thresholds, equations and scientific code branches were not changed. No analysis code was executed.

## Use

Set HCC_PROJECT_ROOT to an existing isolated restored workspace and put this repository root on PYTHONPATH before starting Jupyter or an archived script. The adapter returns a native filesystem string. It refuses the original project, release directory and traversal outside the configured root. It never copies data or writes a file.

Legacy notebooks retain generated-script strings, environment-specific commands and mixed historical/setup sequences. Embedded code may still mention archival Windows/WSL paths; those are explicitly listed by the final safety audit. The adapter fixes direct executable literal paths, not every generated external program or historical saved lock string. Consult the step route and inspect a specific module before execution. No Run-All certification or from-scratch driver is implied.

Step34's accepted native excerpt intentionally omits rejected surrogate training and original filesystem-migration cells. Its baseline inputs and variable setup still require reconstruction in an isolated workspace. Do not run the full historical notebook merely to fill those setup gaps.
