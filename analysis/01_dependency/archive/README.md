# Historical material — not the current analysis

The ZIP archives preserve the old summary-level organoid inference (including the
superseded 37-gene bounds and the original, incorrect 1,653-gene input). They are
retained only for provenance. Do not combine their results with the current
1,720-gene, individual-score analysis.

`previous_figures.zip` contains the byte-identical figures supplied in the conversation.
The old cell-line/integration figures contained overlapping labels. Their caption
also used 2/63 as a common frequency boundary despite gene-specific missingness,
and quoted the full 466-gene union although only 462 candidates were jointly
evaluable. These graphics are not the canonical figures under `figures/svg/`.
The underlying statistical tables were not changed for this presentation fix.

`original_code/reproduce_analysis_1165_original.py` is the unmodified statistical
script recovered from the prior reproducibility archive. The executable repository
port is `scripts/reproduce_cell_lines.py`; only I/O paths were changed. The older
250-model workbooks in `data/inputs/` are audit/template inputs required by that
script, not the primary dataset. All comparative statistics are recalculated from
the 1,165-model matrix.
