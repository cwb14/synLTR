# Memo: prompt for AI to add microsatellite filter to TSD detection

**Date**: 2026-04-22
**Author/session**: Claude Code conversation (resume with the conversation name used in this dir if needed)
**Working dir**: `/data/chris/temp/`
**Purpose**: Hand-off prompt for a future AI task. The goal is to enhance the
stringent TSD detector in `synLTR/module2/ltrharvest5.py` so it rejects
"TSDs" that are really slices of a surrounding microsatellite tract. The
core longest-exact-substring approach is to be preserved; filters are
added on top.

## Why this is needed

While inspecting a 6-element putative nested-LTR stack on `NC_135853.1`
~8.94-8.98 Mb in *Wolffia australiana* (`/home/chris/data/temp/Waustraliana.fa`),
every inner call had a "TSD" reported (`CTCTC`, `CTCTCTC`, `TCTCTC`, ...)
but every outer-flank 20-bp window had dinucleotide entropy ≤ 1.95 bits
(random DNA ≈ 3.8-4.0 bits). The region is a CT-microsatellite hotspot
and the reported TSDs are just slices of the repeat. The current
`_find_tsd_seq` (at `ltrharvest5.py:2476`) only requires ≥2 distinct bases
in the k-mer, which `CTCTCTC` passes.

Contrast: the sole real-looking call at `8939955-8981795` has outer-flank
entropies 2.53 and 2.98 bits and a non-trivial motif `TCTCTTT`. That one
should still pass.

## Input context for the AI

| Item | Value |
|---|---|
| Repo | `/data/chris/temp/synLTR/` |
| File to edit | `synLTR/module2/ltrharvest5.py` |
| Core function | `_find_tsd_seq(left, right, min_len=5)` @ line 2476 |
| Wrapper | `_has_exact_tsd` @ line 2503 |
| Callers | `wfa_guided_tsd_names` @ 2507, `rescue_nonautonomous_by_tsd_from_scn` @ 2589, `tsd_positive_full_length_from_scn` @ 2677 |
| CLI arg | `--tsd-min-len` (used @ 3344, 3446) |
| Genome for verification | `/home/chris/data/temp/Waustraliana.fa` |
| Example failing coords | `NC_135853.1:8946922-8980902`, `:8957380-8970410`, `:8958661-8969155` |
| Example passing coords | `NC_135853.1:8939955-8981795` |

## The prompt (hand this to the AI verbatim)

> **Task**: Enhance TSD (target-site duplication) detection in `synLTR/module2/ltrharvest5.py` so that it rejects "TSDs" that are actually artifacts of surrounding microsatellite tracts. Do **not** change the core detection strategy (longest exact shared substring between outer-left and outer-right flanks). Add filters on top of it.
>
> **Context on TSDs**: A real LTR-RT insertion duplicates a short (typically 4-6 bp) target site, so identical short sequences flank the element on either side. The pipeline's existing detector is `_find_tsd_seq(left, right, min_len=5)` at `synLTR/module2/ltrharvest5.py:2476`, called from `wfa_guided_tsd_names` (`:2507`), `rescue_nonautonomous_by_tsd_from_scn` (`:2589`), and `tsd_positive_full_length_from_scn` (`:2677`). It returns the longest exact shared k-mer (>= min_len, requires >=2 distinct bases) between two short flank windows.
>
> **The failure mode**: When an LTR-RT is called inside or adjacent to a simple-repeat tract (most commonly CT/TC, CAG, or homopolymer), every base position on both sides of the element is part of the same microsatellite. Any k-mer like `CTCTCTC` will match on both sides by pure coincidence, so the detector reports a "TSD" that carries no insertion evidence. The existing >=2-distinct-bases rule is not enough -- `CTCTCTC` passes it.
>
> **Concrete failing cases** (checked by hand on `NC_135853.1` ~8.94-8.98 Mb in *Wolffia australiana*):
>
> | Element coords | Reported TSD | Outer-left 20-bp flank | Outer-right 20-bp flank | Dinuc entropy L / R (bits) |
> |---|---|---|---|---|
> | 8946922-8980902 | `CTCTC`   | `TCTCTCTCTCTCTCTCTCTC` | `CCTCTCTCTCTCTCTCTCTC` | 1.00 / 1.24 |
> | 8957380-8970410 | `CTCTCTC` | `CTCTCTCTCTCTCTCTCTCT` | `TCTCTCCCCCTCGCTCGCTC` | 1.00 / 2.19 |
> | 8958661-8969155 | `TCTCTC`  | `TCTCTCTCTCTCTCTCTCTC` | `CTCTCTCTCTCTCCCTCTTT` | 1.00 / 1.74 |
>
> Random DNA has dinucleotide entropy ≈ 3.8-4.0 bits; these flanks are near 1.0. The "TSDs" are just slices of the microsatellite and should be rejected. Contrast with a real-looking case on the same chromosome (`8939955-8981795`) where both outer flanks have entropy 2.5-3.0 and the candidate TSD motif `TCTCTTT` is not a trivial repeat.
>
> **What to implement**: In `_find_tsd_seq` (or via a helper it calls), before returning a motif, apply all of the following rejection criteria. A candidate fails if ANY criterion triggers:
>
> 1. **TSD-motif is a simple repeat.** Reject if the motif itself is a pure homopolymer, a pure dinucleotide tandem (e.g. `CTCTCT`, `TATATA`), or a pure trinucleotide tandem (e.g. `CAGCAGCAG`). Implement as: the motif is a simple repeat if it can be expressed as `unit * n` where `unit` has length 1, 2, or 3 and `n >= 2`.
>
> 2. **Flank-context is low-complexity.** For each flank window, compute Shannon dinucleotide entropy over a >=20-bp context centered at the boundary (extending beyond the short flank window currently passed to `_find_tsd_seq` -- the caller will need to provide these longer windows, or the helper should take the genome seq + boundary coords). Reject if either flank's entropy is below a threshold (start at **2.0 bits**; make it a configurable parameter).
>
> 3. **Motif is not unique in its flank.** In each 20-bp flank context, count exact occurrences of the candidate TSD motif. A real TSD appears essentially once per flank, at the position immediately abutting the LTR. Reject if the motif occurs **> 2 times in either flank's 20-bp context** (i.e. the "TSD" is just another instance of the local repeat unit).
>
> 4. **Extensibility check.** Check whether the matched motif extends into a much longer exact match (e.g. >= 12 bp) by growing the match on both sides. If the exact match extends >= 12 bp in either direction inside the flank context, the underlying sequence is repetitive rather than a TSD -- reject.
>
> **Preserve**:
> - The longest-first search.
> - The exact-match requirement (no mismatches).
> - The >=2-distinct-bases rule (now superseded by #1 but keep it as a guard).
> - `min_len` remains 5 by default.
> - Public signatures of `_find_tsd_seq`, `_has_exact_tsd`, `wfa_guided_tsd_names`, `rescue_nonautonomous_by_tsd_from_scn`, `tsd_positive_full_length_from_scn` -- extend with optional new kwargs (e.g. `min_flank_entropy=2.0`, `max_motif_occurrences=2`, `max_extensibility=12`, and a way to pass in the longer flank-context strings), default values should make the new filters active but conservative.
>
> **Output requirements**:
> - Add a short module-level docstring note explaining the microsatellite filter.
> - Update call sites so that when a longer flank context is available (the caller has the full genome sequence), the context is passed down; otherwise fall back to current behavior on just the short windows (the filters still run where they can: #1 always applies; #2/#3/#4 need the long context).
> - Add pytest-style tests in a new file `synLTR/module2/tests/test_tsd_microsat.py`. Include at minimum:
>   - a passing real-looking TSD case (non-repeat flanks, motif `ATGCA`),
>   - a failing pure-CT-microsatellite case (flanks `TCTCTCTCTCTCTCTCTCTC`, reported motif `CTCTCTC`),
>   - a failing trinucleotide-tandem case (flanks `CAGCAGCAG...`),
>   - a boundary case where the motif is non-trivial but flank entropy is borderline (document the expected behavior).
> - Do not change behavior of any function that doesn't touch TSDs.
>
> **Verification**: After the change, run the three rejection cases above (B/D/E at the listed coords) through the updated `_find_tsd_seq` with the true genomic flanks (extract from `/home/chris/data/temp/Waustraliana.fa`) and confirm all three now return `None`. Confirm the high-complexity case at `8939955-8981795` still returns a non-`None` motif.

## Suggested adjustments before sending

- Tighten/loosen thresholds (2.0 bits / 2 occurrences / 12 bp) to taste.
- Minimum-viable version: keep only criterion #1 (simple-repeat motif
  filter). That alone rejects all three failing cases above because
  `CTCTCTC`, `CTCTC`, `TCTCTC` are pure dinucleotide tandems. The other
  three criteria catch subtler cases but require plumbing a longer
  flank context through the call chain.

## Notes

- The 6-element nested-LTR stack that exposed this bug is on
  `NC_135853.1:8939955-8981795` in the depth4 bucket of
  `mafft_update_reconcile_depth*_ltr.{tsv,fa}` (also see the raw
  per-round files `mafft_update_r{1,2,3}_ltr.tsv`).
- None of the five inner calls had TEsorter domain hits, and the
  reported K2P divergences are biologically inconsistent (inner > outer),
  reinforcing that most of them are microsatellite-driven false positives.
- The plot `/data/chris/temp/nested_ltrs_NC_135853_8.94-8.98Mb.pdf`
  visualizes the stack for reference.
