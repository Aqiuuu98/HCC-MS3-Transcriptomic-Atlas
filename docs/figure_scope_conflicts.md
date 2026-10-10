# Figure S7A: confirmed displayed-scope and frozen-legend mismatch

The saved Step30 robust table has28 rows: MS6=27, MS3=1. The ≥3% multi-resource provenance table has39 rows: MS6=34, MS5=4, MS3=1. The additional11 are secondary moderate candidates (MS6=7, MS5=4). Margins describe relative state skew, not percentages of supporting patients. Formal primary comparisons use HCC09T/HCC10T; HCC10L is a same-patient site sensitivity.

Original Step30 cell185 reads the 39-row provenance table; cell189 groups it by target and evidence tier to draw the count plot. The selected source PDF, final FigureS7 PNG and native PPTX labels agree with that two-tier display. The final visible bars are MS6 27+7, MS5 0+4 and MS3 1+0.

The frozen English Supplementary document's Figure S7 legend (paragraph57, zero-based) describes all28 robust candidates and says edge widths encode counts. The displayed panel is a two-tier count bar plot, so both scope and visual-description wording are mismatched. Main results and Supplementary methods/results correctly describe both28 and39.

Classification: CONFIRMED_FIGURE_LEGEND_VERSION_MISMATCH. This blocks a public release represented as fully synchronized with the current manuscript, but does not invalidate the saved tier counts or prevent a clearly labeled private code/source review. No figure, Word document or frozen scientific result was modified. The current documentation update proposes an accurate caption for the existing 39-candidate display. The proposal has not been applied to the frozen manuscript; its synchronization still requires author approval.

## Proposed replacement for Figure S7A, awaiting author approval

A, Counts of all 39 frozen state-skewed ligand–receptor candidates, grouped by malignant receiver state and evidence tier. Blue bars show 28 primary robust candidates (27 MS6-associated and one MS3-associated), with cross-patient minimum relative state-skew margins ≥5%. Orange bars show 11 secondary moderate candidates (seven MS6-associated and four MS5-associated), with margins ≥3% and <5%. Bar heights indicate retained candidate counts, rather than measured signaling strength. These descriptive effect-size tiers are not percentages of supporting patients or statistical significance thresholds. Candidates were inferred from HCC09T/HCC10T primary-sample comparisons; HCC10L served only as a same-patient sensitivity context, not a third independent patient. Inferred communication does not establish physical signaling or causal pathway activation.

Replace only the panel A description in the frozen Figure S7 legend. Blue/orange encode primary robust/secondary moderate evidence tiers; bar heights encode counts. Panel A contains no network nodes or edges, so its former edge-width sentence is inapplicable. The panel B text is outside this correction and remains unchanged.

The frozen Supplementary methods/results (paragraphs20/41) and main Results (paragraph87; zero-based) already describe39 at ≥3% and28 at ≥5%; no associated numerical Results change is needed. Verified source: Step30 notebook cells59/185/189, the 39-row Step30E1 provenance table, the 28-row Step30D1 robust table, final FigureS7.png and frozen Supplementary legend paragraph57. No analyses were rerun, no candidates were changed, and no figure or Word file was written.
