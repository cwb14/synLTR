# Memo: combined prompt — microsatellite-aware TSD filter + solo-LTR TSD scan

**Date**: 2026-04-22
**Working dir**: `/data/chris/temp/`
**Purpose**: One hand-off prompt that addresses both TSD-detection
issues in a single AI conversation:

1. **Microsatellite filter** — the existing TSD detector is fooled by
   simple-repeat flanks (CT/TC/CAG/homopolymer), reporting "TSDs" that
   are just slices of the local microsatellite.
2. **Solo-LTR TSD scan** — the pipeline only checks for TSDs flanking
   the whole element; it does not notice when TSDs hug each LTR
   individually, which is the signature of a solo-LTR misclassified as
   one end of an intact element (or of a microsatellite-driven
   artifact).

The two enhancements compose: the microsatellite filter makes every TSD
probe (full-length and solo) more trustworthy; the solo scan turns
those trustworthy probes into a structural diagnostic that catches
false-positive full-length calls.

## Why this is needed

### (1) Microsatellite false-positive TSDs

On `NC_135853.1` ~8.94-8.98 Mb in *Wolffia australiana* every inner
LTR-RT call in a 6-element stack had a "TSD" reported but every outer
flank had dinucleotide entropy <= 1.95 bits (random DNA ~3.8-4.0 bits).
The reported motifs (`CTCTC`, `CTCTCTC`, `TCTCTC`, ...) are simple
dinucleotide tandems that happen to appear on both sides of any
boundary drawn through the CT-microsatellite tract. The existing
>=2-distinct-bases rule in `_find_tsd_seq` passes them.

### (2) Solo-LTR / paired-solo false-positive full-length calls

In the same stack, two candidates (`NC_135853.1:8957380-8970410` and
`NC_135853.1:8957562-8976405`) share a near-identical left LTR (~1 kb of
overlap) but have right LTRs ~5 kb apart. Neither has TEsorter domain
hits. K2P divergences are biologically inconsistent (inner older than
outer). The most parsimonious interpretation: a genuine solo-LTR at the
shared-left-LTR position got paired by the annotator with two different
spurious right-LTR candidates (both within a microsatellite tract),
yielding two bogus "full-length" elements. A solo-LTR TSD scan would
have flagged both.

## Input context for the AI

| Item | Value |
|---|---|
| Repo | `/data/chris/temp/synLTR/` |
| File to edit | `synLTR/module2/ltrharvest5.py` |
| Core TSD helper | `_find_tsd_seq(left, right, min_len=5)` @ line 2476 |
| TSD wrapper | `_has_exact_tsd` @ line 2503 |
| Existing full-length TSD scans | `wfa_guided_tsd_names` @ 2507, `rescue_nonautonomous_by_tsd_from_scn` @ 2589, `tsd_positive_full_length_from_scn` @ 2677 |
| SCN LTR-boundary loader | `load_scn_ltr_boundaries` @ 1858 |
| Dedup / output writer | `dedup_kmer2ltr_tsv` @ 1952 (TSV header written @ 2310) |
| CLI arg | `--tsd-min-len` (used @ 3344, 3446) |
| Genome for verification | `/home/chris/data/temp/Waustraliana.fa` |
| Example microsat-fail coords | `NC_135853.1:8946922-8980902`, `:8957380-8970410`, `:8958661-8969155` |
| Example solo-pair | `NC_135853.1:8957380-8970410` and `NC_135853.1:8957562-8976405` |
| Intact-looking contrast | `NC_135853.1:8939955-8981795` |

## The combined prompt (hand this to the AI verbatim)

> **Task**: Make two composable enhancements to TSD detection in `synLTR/module2/ltrharvest5.py`:
>
> **(A) Microsatellite-aware filter** inside the existing TSD detector, so that "TSDs" which are just slices of a simple-repeat tract are rejected.
>
> **(B) Solo-LTR TSD scan** alongside the existing full-length TSD check, producing a per-candidate topology tag (`intact` / `solo-left` / `solo-right` / `solo-both` / `ambiguous` / `none`) that is emitted in the output TSV.
>
> Do not alter the stringent core approach (longest exact shared k-mer, no mismatches, longest-first search). Do not remove any existing output columns; the two new columns from (B) are added next to `tsd`. Do not auto-discard candidates in (B) — the tag is diagnostic.
>
> ---
>
> ### Background on TSDs
>
> A real LTR-RT insertion duplicates a short (typically 4-6 bp) target site, so identical short sequences flank the element on either side immediately 5' of `lLTR.start` and immediately 3' of `rLTR.end`. A *solo-LTR* (post-recombination remnant) still carries the original TSDs but they now hug a single LTR: identical 4-6 bp motifs sit immediately 5' of the LTR and immediately 3' of the same LTR. Microsatellite tracts (CT/TC/CAG/homopolymer) produce spurious "TSDs" at every boundary because the k-mer search is swamped by the tract's inherent self-similarity.
>
> ### Failing evidence (*Wolffia australiana*, `NC_135853.1`)
>
> Microsatellite false positives (checked with 20-bp flank windows):
>
> | Element coords | Reported TSD | Outer-left flank | Outer-right flank | Dinuc entropy L / R (bits) |
> |---|---|---|---|---|
> | 8946922-8980902 | `CTCTC`   | `TCTCTCTCTCTCTCTCTCTC` | `CCTCTCTCTCTCTCTCTCTC` | 1.00 / 1.24 |
> | 8957380-8970410 | `CTCTCTC` | `CTCTCTCTCTCTCTCTCTCT` | `TCTCTCCCCCTCGCTCGCTC` | 1.00 / 2.19 |
> | 8958661-8969155 | `TCTCTC`  | `TCTCTCTCTCTCTCTCTCTC` | `CTCTCTCTCTCTCCCTCTTT` | 1.00 / 1.74 |
>
> Contrast `8939955-8981795`: both outer-flank entropies 2.5-3.0, candidate motif `TCTCTTT` (not a trivial repeat) — this one should still pass.
>
> Solo-LTR false-positive pairing: `8957380-8970410` and `8957562-8976405` share a near-identical left LTR (~1 kb overlap) but have right LTRs ~5 kb apart. The shared left-LTR's own flanks (just 5' of `lLTR.start` and just 3' of `lLTR.end`) should be scanned for TSDs; a match there indicates a solo-LTR at that position and undermines both "full-length" calls.
>
> ---
>
> ### (A) Microsatellite-aware filter — what to implement
>
> Refactor `_find_tsd_seq` (line 2476) so that, before returning a motif, it applies ALL of the following rejection criteria. A candidate fails if ANY criterion triggers:
>
> 1. **TSD-motif is a simple repeat.** Reject if the motif can be expressed as `unit * n` where `len(unit) in {1, 2, 3}` and `n >= 2`. This alone rejects `CTCTCTC`, `CTCTC`, `TCTCTC`, `CAGCAGCAG`, `AAAAA`, etc. Implement as a helper `_is_simple_tandem(motif) -> bool`.
>
> 2. **Flank-context is low-complexity.** For each flank, compute Shannon dinucleotide entropy over a >=20-bp context centered at the boundary. Reject if either context's entropy is below a threshold (default `min_flank_entropy=2.0`; expose as kwarg). Implement as a helper `_dinuc_entropy(seq) -> float`. This criterion needs access to the longer context; see "caller plumbing" below.
>
> 3. **Motif is not unique in its flank.** In each 20-bp flank context, count exact occurrences of the candidate motif. Reject if the motif occurs `> max_motif_occurrences` times in either context (default 2). Implement as a helper `_count_occurrences(motif, seq) -> int`.
>
> 4. **Extensibility check.** From the matched k-mer's position in each flank, greedily extend the exact match inward and outward. If the extended exact match reaches `>= max_extensibility` bp in either flank, the sequence is repetitive rather than a TSD — reject. Default `max_extensibility=12`.
>
> **Caller plumbing**: Criteria 2-4 need a longer flank *context* than the short (~6-7 bp) windows some callers currently pass. Extend the signature:
>
> ```python
> def _find_tsd_seq(
>     left: str,
>     right: str,
>     min_len: int = 5,
>     left_context: Optional[str] = None,
>     right_context: Optional[str] = None,
>     min_flank_entropy: float = 2.0,
>     max_motif_occurrences: int = 2,
>     max_extensibility: int = 12,
> ) -> Optional[str]:
>     ...
> ```
>
> Keep the old behavior as the fallback when `left_context`/`right_context` is `None` (criterion #1 still applies; 2-4 skipped).
>
> Update each caller (`wfa_guided_tsd_names`, `rescue_nonautonomous_by_tsd_from_scn`, `tsd_positive_full_length_from_scn`) to pass 20-bp windows for context (carving them from the full `chrom` sequence they already load). Keep `_has_exact_tsd` delegating to the new signature.
>
> **Preserve**:
> - Longest-first search.
> - Exact-match requirement (no mismatches).
> - `>= 2 distinct bases` rule (kept as a cheap guard even though #1 subsumes it).
> - `min_len=5` default.
>
> ---
>
> ### (B) Solo-LTR TSD scan — what to implement
>
> Add a new function:
>
> ```python
> def scan_tsd_topology_from_scn(
>     stitched_scn: str,
>     genome_fa: str,
>     flank: int = 20,
>     min_len: int = 5,
> ) -> Dict[str, Dict[str, object]]:
>     """te_key -> {tsd_full, tsd_soloL, tsd_soloR, topology}"""
> ```
>
> Uses `_find_tsd_seq` (with the criteria from part A active) against three probe pairs derived from the SCN record's LTR boundaries (1-based inclusive: `lLTR_s, lLTR_e, rLTR_s, rLTR_e`):
>
> - `full_left`   = seq[`lLTR_s - 1 - flank` : `lLTR_s - 1`]
> - `full_right`  = seq[`rLTR_e`             : `rLTR_e + flank`]
> - `solo_L_in`   = seq[`lLTR_e`             : `lLTR_e + flank`]
> - `solo_R_in`   = seq[`rLTR_s - 1 - flank` : `rLTR_s - 1`]
>
> Probes:
>
> - `tsd_full`  = `_find_tsd_seq(full_left,  full_right,  ...)` -- intact signature
> - `tsd_soloL` = `_find_tsd_seq(full_left,  solo_L_in,   ...)` -- solo around left LTR
> - `tsd_soloR` = `_find_tsd_seq(solo_R_in,  full_right,  ...)` -- solo around right LTR
>
> Pass `left_context` and `right_context` to each probe so the microsatellite filters in (A) apply; these contexts can be the same 20-bp windows extended a further 10 bp for entropy calculation if desired (parameterize).
>
> Topology classification:
>
> - `intact`     : `tsd_full` present AND both solo probes absent.
> - `solo-left`  : `tsd_soloL` present AND others absent.
> - `solo-right` : `tsd_soloR` present AND others absent.
> - `solo-both`  : both solo probes present AND `tsd_full` absent.
> - `ambiguous`  : `tsd_full` present AND one or both solo probes also present (cannot resolve from TSDs alone — typically microsatellite residue; the filter in (A) should have made this rare).
> - `none`       : all three probes return `None`.
>
> ---
>
> ### TSV schema changes
>
> Existing header (`dedup_kmer2ltr_tsv`, line 2310):
>
> ```
> #name LTR_len aln_len subs ti tv raw_d raw_T JC69_d JC69_T K2P_d K2P_T left_trim right_trim tsd domains nest_status
> ```
>
> New header: insert two columns `tsd_topology` and `tsd_solo_motifs` between `tsd` and `domains`:
>
> ```
> #name ... tsd tsd_topology tsd_solo_motifs domains nest_status
> ```
>
> `tsd_topology` := one of the six tags above.
> `tsd_solo_motifs` := `L=<motif>;R=<motif>` with `.` for an absent side, or `.` when both are absent.
>
> Wire `scan_tsd_topology_from_scn` in `main()` (around line 2772; near where the merged SCN and full-length TSD scan are computed) and thread its output into `dedup_kmer2ltr_tsv` so the new columns are populated per element.
>
> Downstream:
>
> - `names_from_kmer2ltr_dedup` (line 2395) reads col1 only — unaffected.
> - `reconcile_nests.py` (if present in the repo) rewrites the last column (`nest_status`); ensure it still picks the LAST column by index, not a hardcoded column number. If it uses `cols[-1]`, that still works.
> - Check any other scripts that consume the TSV and fail loudly with a clear message if a parser hardcodes column positions that would shift.
>
> ---
>
> ### Do not change
>
> - Public signatures of `_has_exact_tsd`, `wfa_guided_tsd_names`, `rescue_nonautonomous_by_tsd_from_scn`, `tsd_positive_full_length_from_scn` beyond the optional kwargs added above. Defaults must keep prior behavior for callers that don't opt in.
> - Any function that doesn't touch TSDs or the TSV schema.
>
> ---
>
> ### Tests (new file `synLTR/module2/tests/test_tsd_enhancements.py`)
>
> Microsatellite filter (part A):
>
> - Passing: non-repeat flanks `CCCTCGCTGATCTCTTTCGC` and `CTCCATTCTCCCTCTCTTTG`, expected motif `TCTCTTT` (len 7).
> - Failing (#1 simple tandem): flanks `TCTCTCTCTCTCTCTCTCTC` and `CCTCTCTCTCTCTCTCTCTC`, naive match would be `TCTCTCTC`; expect `None`.
> - Failing (#1): flanks containing `CAGCAGCAGCAG...`; expect `None`.
> - Failing (#2 low entropy): flanks with dinucleotide entropy 1.5 bits but a non-tandem shared k-mer (e.g., inject `ATGCA` once into otherwise pure-CT context); expect `None` due to entropy threshold.
> - Failing (#3 multi-occurrence): flanks where the candidate motif appears 3+ times each; expect `None`.
> - Failing (#4 extensibility): match extends >= 12 bp contiguously; expect `None`.
>
> Solo-LTR topology (part B):
>
> - `intact` case: distinct non-tandem TSD motif (e.g. `ATGCA`) flanking only the outer boundaries; both solo probes see random, non-matching sequence. Expect `topology == "intact"`.
> - `solo-left` case: motif flanks the left LTR only (outer-L <-> inner-L); outer-R is random.
> - `solo-right` case: symmetric.
> - `solo-both` case: same motif hugs each LTR individually on the inner side; outer-L vs outer-R fails (different regions), so `tsd_full is None` but both solos hit.
> - `ambiguous` case: motif at all four boundaries AND all flanks high-entropy (not a microsatellite). Expect `topology == "ambiguous"`.
> - `none` case: all random. Expect `topology == "none"`.
>
> Integration:
>
> - Build a synthetic SCN file with 2 records and a small FASTA, run `scan_tsd_topology_from_scn`, assert returned tags.
>
> ---
>
> ### Verification on real data
>
> After implementation:
>
> 1. Extract `NC_135853.1` from `/home/chris/data/temp/Waustraliana.fa`.
> 2. Run `_find_tsd_seq` with 20-bp contexts at the three microsatellite failing cases:
>    - `8946922-8980902`, `8957380-8970410`, `8958661-8969155` — expect all three to return `None`.
> 3. Run `_find_tsd_seq` at the intact-looking case `8939955-8981795` (outer flanks `CCCTCGCTGATCTCTTTCGC` / `CTCCATTCTCCCTCTCTTTG`) — expect a non-`None` motif (should be `TCTCTTT` or similar 7-8 bp non-repeat string).
> 4. Run `scan_tsd_topology_from_scn` on the full stitched SCN for the failing region and print the six topology tags for the stack at 8.94-8.98 Mb. At least one of `8957380-8970410` / `8957562-8976405` should NOT tag `intact`. `8939955-8981795` is permitted to tag `intact`, `none`, or `ambiguous` but should not tag `solo-*`.
>
> Print the verification summary as part of the run so it shows up in logs.
>
> ---
>
> ### Deliverable checklist
>
> - [ ] `_find_tsd_seq` refactored with 4 criteria and optional long-context kwargs. Old behavior preserved when no context is supplied.
> - [ ] `_is_simple_tandem`, `_dinuc_entropy`, `_count_occurrences` helpers added and used.
> - [ ] All existing callers of `_find_tsd_seq` / `_has_exact_tsd` updated to pass long contexts where they already have the chromosome sequence loaded.
> - [ ] `scan_tsd_topology_from_scn` added and wired into `main()`.
> - [ ] `dedup_kmer2ltr_tsv` updated to accept topology info and emit `tsd_topology` + `tsd_solo_motifs` between `tsd` and `domains`.
> - [ ] TSV header updated.
> - [ ] `reconcile_nests.py` audited for positional column assumptions (should still work because it uses `cols[-1]`).
> - [ ] Tests added under `synLTR/module2/tests/`.
> - [ ] Verification run on the *Wolffia* data prints the expected results.
> - [ ] A short note added to the top-of-file docstring in `ltrharvest5.py` summarizing both enhancements.

## Suggested adjustments before sending

- Thresholds are starting points; tune on a dataset with known true
  TSDs + known false-positive regions.
- If the long-context plumbing feels heavy, the minimum-viable version
  is criterion #1 alone in (A) plus the full (B) scan. Criterion #1
  already rejects the three microsatellite failing cases, because
  `CTCTCTC` / `CTCTC` / `TCTCTC` are all pure dinucleotide tandems.
- The new columns are placed *between* `tsd` and `domains` to keep
  related fields adjacent. Trade-off: any downstream parser hardcoding
  `col[14]` for `tsd` still works, but `col[15]` (was `domains`) now
  means `tsd_topology`. Audit accordingly. If that's too risky, move
  both new columns to end-of-row instead.

## Notes

- The microsatellite-only version of this prompt lives in
  `memo_tsd_microsat_filter_prompt.md`.
- The solo-only version lives in `memo_tsd_solo_ltr_scan_prompt.md`.
- Plot of the stack that motivated both enhancements:
  `/data/chris/temp/nested_ltrs_NC_135853_8.94-8.98Mb.pdf`.
