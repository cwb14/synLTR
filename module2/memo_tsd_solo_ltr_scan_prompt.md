# Memo: prompt for AI to add solo-LTR TSD scan to LTR-RT validation

**Date**: 2026-04-22
**Working dir**: `/data/chris/temp/`
**Purpose**: Hand-off prompt for an AI task. Add a *solo-LTR TSD scan*
alongside the existing full-length TSD check, so we can flag SCN LTR-RT
candidates where the presence of TSDs around either individual LTR (not
around the whole element) signals (1) a solo-LTR misclassified as one
end of an intact LTR-RT, or (2) microsatellite sequence being called as
LTR-RT.

## Why this is needed

A real full-length LTR-RT has TSDs flanking the *entire element* (just
outside `lLTR.start` and just after `rLTR.end`). A **solo-LTR** — the
remnant left after intra-element LTR-LTR recombination excises the
internal region and one LTR — still carries the original TSDs, but they
now flank the single remaining LTR on its own. The current pipeline only
looks for full-length TSDs, so it cannot distinguish:

- **Intact element**: TSD at outer boundaries only.
- **Solo-LTR dressed as one end of a bogus pair**: TSD hugging just the
  real LTR. The annotator paired it with a spurious second LTR (often
  another solo-LTR or a microsatellite-driven match) to fabricate a
  full-length call. Scanning solo-LTR-style TSDs around each LTR exposes
  the pattern.
- **Microsatellite artifact**: TSD-like matches are found at *every*
  boundary because the whole region is a repeat tract. This is the
  degenerate case where both solo-LTR probes and the full-length probe
  all return "hits."

Evidence from `NC_135853.1:8957380-8970410` (*Wolffia australiana*,
`/home/chris/data/temp/Waustraliana.fa`) and its sibling
`8957562-8976405`: these two calls share a near-identical left LTR (~1 kb
of overlap) but have right LTRs ~5 kb apart. A solo-LTR TSD probe at the
shared left-LTR position returns matches (consistent with either a
solo-LTR or a microsatellite), while the "full-length" element calls
lack convincing domain support and sit in CT-microsatellite-rich
flanking DNA. A solo-LTR scan run alongside the usual TSD check would
have flagged both as suspect.

## Input context for the AI

| Item | Value |
|---|---|
| Repo | `/data/chris/temp/synLTR/` |
| File to edit | `synLTR/module2/ltrharvest5.py` |
| TSD core helper | `_find_tsd_seq(left, right, min_len=5)` @ line 2476 |
| TSD wrapper | `_has_exact_tsd` @ line 2503 |
| Existing full-length TSD scans | `wfa_guided_tsd_names` @ 2507, `rescue_nonautonomous_by_tsd_from_scn` @ 2589, `tsd_positive_full_length_from_scn` @ 2677 |
| SCN LTR-boundary loader | `load_scn_ltr_boundaries` @ 1858 |
| SCN file path pattern | `{prefix}.work/{prefix}.ltrtools.stitched.scn` |
| Output TSV schema | `#name LTR_len aln_len subs ti tv raw_d raw_T JC69_d JC69_T K2P_d K2P_T left_trim right_trim tsd domains nest_status` (see `dedup_kmer2ltr_tsv` @ 1952, header line written @ 2310) |
| Genome for verification | `/home/chris/data/temp/Waustraliana.fa` |
| Failing sibling pair | `NC_135853.1:8957380-8970410` and `NC_135853.1:8957562-8976405` |
| Intact-looking contrast | `NC_135853.1:8939955-8981795` (still should not trigger solo flags) |

## The prompt (hand this to the AI verbatim)

> **Task**: Add a *solo-LTR TSD scan* to `synLTR/module2/ltrharvest5.py` that, for each SCN LTR-RT candidate, probes for TSDs hugging each individual LTR (not only the full element). Use the result to flag candidates whose TSD topology is inconsistent with an intact full-length insertion. Do **not** remove or change the existing full-length TSD detection; add the new scan alongside it.
>
> **Biological rationale**: A real full-length LTR-RT has TSDs flanking the whole element: a 4-6 bp duplicated target sequence sits immediately 5' of `lLTR.start` and immediately 3' of `rLTR.end`. A *solo-LTR* (post-recombination remnant) carries TSDs hugging just the one remaining LTR, i.e., identical 4-6 bp motifs immediately 5' of the LTR and immediately 3' of the LTR. If a candidate reported by ltrharvest/LTRfinder contains a *solo-LTR signature around one of its two LTRs*, the call is suspect: either the tool paired a genuine solo-LTR with a spurious partner to fabricate a full-length element, or the whole region is a microsatellite where "TSDs" appear at every boundary.
>
> **What to implement**:
>
> 1. **Helper to extract boundary flanks from the genome** (or reuse existing helpers). Given an SCN record with 1-based inclusive coordinates (`lLTR_s, lLTR_e, rLTR_s, rLTR_e`) and a chromosome sequence:
>
>    - `full_left_flank`  = seq[`lLTR_s-1-W` : `lLTR_s-1`]                 # bp before lLTR start
>    - `full_right_flank` = seq[`rLTR_e`    : `rLTR_e + W`]                 # bp after rLTR end
>    - `solo_L_inner_flank` = seq[`lLTR_e`  : `lLTR_e + W`]                 # bp just after lLTR end
>    - `solo_R_inner_flank` = seq[`rLTR_s-1-W` : `rLTR_s-1`]                # bp just before rLTR start
>
>    where `W` is the flank-scan window (keep consistent with the existing scan; expose as a parameter, default 20).
>
> 2. **Scan three TSD probes per candidate** using the existing `_find_tsd_seq` (with whatever microsatellite-aware filters are currently active — this prompt is agnostic to that):
>
>    - `tsd_full`   = `_find_tsd_seq(full_left_flank,  full_right_flank)`   # intact-element TSD
>    - `tsd_soloL`  = `_find_tsd_seq(full_left_flank,  solo_L_inner_flank)` # TSD hugging left LTR alone
>    - `tsd_soloR`  = `_find_tsd_seq(solo_R_inner_flank, full_right_flank)` # TSD hugging right LTR alone
>
> 3. **Classify each candidate**'s TSD topology into one of these tags (string literal):
>
>    - `intact`         : `tsd_full` present AND both `tsd_soloL` and `tsd_soloR` absent. Good.
>    - `solo-left`      : `tsd_soloL` present AND `tsd_full` absent AND `tsd_soloR` absent. Likely a solo-LTR at the left boundary with a spurious right partner.
>    - `solo-right`     : `tsd_soloR` present AND `tsd_full` absent AND `tsd_soloL` absent. Likely a solo-LTR at the right boundary with a spurious left partner.
>    - `solo-both`      : both `tsd_soloL` and `tsd_soloR` present AND `tsd_full` absent. Two separate solo-LTRs paired by the tool, or a microsatellite hotspot.
>    - `ambiguous`      : `tsd_full` present AND one or both solo probes also present. Cannot resolve from TSDs alone (most commonly microsatellite-driven).
>    - `none`           : no TSD in any of the three probes.
>
> 4. **Expose the scan and the tag**:
>
>    - Add a new function, e.g. `scan_tsd_topology_from_scn(stitched_scn, genome_fa, flank=20, min_len=5) -> Dict[te_key, dict]` where `te_key` is `"chrom:s-e"` (consistent with the rest of the code) and the value is a dict `{"tsd_full", "tsd_soloL", "tsd_soloR", "topology"}`.
>    - Wire it into the same stage of the pipeline where the existing full-length TSD scan runs (see `main()` @ line 2772; step ~9b is where `dedup_kmer2ltr_tsv` is called). Results should be threaded into `dedup_kmer2ltr_tsv` and written to the output TSV.
>
> 5. **Extend the output TSV schema** by appending two columns *after* the existing `tsd` column (and updating the header @ line 2310 and any downstream parsers):
>
>    - `tsd_topology` : one of the six tag strings above.
>    - `tsd_solo_motifs` : concise summary like `L=CTCTC;R=.` (empty side shown as `.`), or `.` if neither solo probe returned a motif.
>
>    `dedup_kmer2ltr_tsv` (line 1952) currently writes `tsd`, `domains`, `nest_status`. Add the two new columns between `tsd` and `domains` so the stable tail (`domains`, `nest_status`) stays at end-of-row.
>
> 6. **Do not auto-discard candidates** based on this tag. The tag is diagnostic; downstream filtering is the user's call. (If you add any filtering it must be behind a default-off CLI flag, e.g. `--drop-non-intact-tsd-topology`.)
>
> **Preserve**:
> - The existing `_find_tsd_seq` / `_has_exact_tsd` signatures and default behavior.
> - All existing output columns in their current order, with the two new columns inserted as specified above.
> - The existing full-length TSD scan path (`tsd_positive_full_length_from_scn`, `wfa_guided_tsd_names`, `rescue_nonautonomous_by_tsd_from_scn`); do not duplicate their work — reuse their results if convenient.
>
> **Tests** (new file `synLTR/module2/tests/test_tsd_topology.py`):
>
> - `intact` case: synthetic 10-kb region with a distinct 5-bp TSD (e.g. `ATGCA`) flanking the whole element, and random (non-repeating) interior sequence + LTR-adjacent bases that do NOT match the TSD. Expect `topology == "intact"`.
> - `solo-left` case: the 5-bp motif hugs the left LTR only; right-LTR flanks are random and do not match each other.
> - `solo-right` case: symmetric to `solo-left`.
> - `solo-both` case: identical 5-bp motif hugs each LTR individually, left and right match each other too but also match the solo probes.
> - `ambiguous` case: same motif appears at all four boundaries (i.e. microsatellite-like context). Expect `topology == "ambiguous"`.
> - `none` case: all flanks random; no motif shared anywhere. Expect `topology == "none"`.
>
> **Verification on real data** (after implementation):
>
> - Extract the genome chromosome `NC_135853.1` from `/home/chris/data/temp/Waustraliana.fa`.
> - Run the new scan on `NC_135853.1:8939955-8981795` (the intact-looking candidate). Expect `topology == "intact"` OR `ambiguous` (acceptable because the pipeline currently reports `tsd=.` for this element, suggesting the full-length motif is subtle).
> - Run on `NC_135853.1:8957380-8970410` and `NC_135853.1:8957562-8976405` (the sibling pair with shared left LTR). Expect `topology ∈ {solo-left, solo-both, ambiguous}` for at least one of them. A tag of `intact` on both is a bug.
> - Print a summary table showing the topology call for every element in the failing stack so the user can spot-check.

## Suggested adjustments before sending

- If you already plumbed a microsatellite filter into `_find_tsd_seq`
  (see `memo_tsd_microsat_filter_prompt.md`), this scan benefits from it
  automatically — a microsatellite flank will no longer produce TSDs at
  every boundary, so `ambiguous` becomes much rarer and the solo-L/R
  distinctions sharpen.
- Consider whether the two new TSV columns should sit between `tsd` and
  `domains` (as specified) or at end-of-row. End-of-row minimizes risk
  of breaking existing parsers but makes grouped reading less
  ergonomic.
- Flank window `W` default 20 bp matches the scripts already used for
  diagnostic prints. If you raise it, keep an eye on solo-probe false
  positives — a larger window gives the k-mer search more chances.

## Notes

- This work composes cleanly with the microsatellite filter in
  `memo_tsd_microsat_filter_prompt.md`. A combined implementation is
  described in `memo_tsd_microsat_and_solo_combined_prompt.md`.
- The 6-element stack that motivated this is plotted in
  `/data/chris/temp/nested_ltrs_NC_135853_8.94-8.98Mb.pdf` and further
  discussed in the session notes.
