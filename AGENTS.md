# Developing synLTR

A guide for developers of this repository, human or AI. Users start at [README.md](README.md).

## What is here

| Path | What it is | Rule |
|---|---|---|
| `src/synltr/` | The package, command `synltr` | All new work goes here |
| `tests/` | pytest suite; its folders mirror `src/synltr/` | Every change comes with a test on a small synthetic fixture |
| `tests/golden/` | Pinned outputs of four end-to-end toy runs | Change only on purpose (see Golden pins) |
| `docs/` | `outputs.md` (every output file), `benchmarks.md` (accuracy, time, memory) | Keep in step with the code |
| `examples/quickstart.sh` | Every command on toy data; also run by `tests/test_quickstart.py` | Keep in step with the README's Run section |
| `legacy/` | The earlier hand-built synLTR pipeline and the tarballs of LTRpav, LTRmix and ltrkin | Read-only reference (see Legacy) |

## Set up and check

```bash
mamba env create -f environment.yml && mamba activate synltr_pkg
pip install --no-deps -e .
ruff check src tests            # lint (line length 120)
pytest -m "not slow"            # fast tests
pytest tests/test_golden.py     # the four pinned end-to-end runs
pytest                          # everything, about 50 min
```

CI (`.github/workflows/tests.yml`) runs the first three checks. pytest turns every warning into an error. A test
skips when a program it needs (lastz, minimap2, samtools, bcftools, poppler) or msprime is missing.

## How the package is built

One command line (`cli.py`, table `COMMANDS`) dispatches to one function per command. Data flows upward through
layers, and `tests/test_layers.py` fails on any import that goes down or sideways:

| Rank | Units | Role |
|---|---|---|
| 0 | `base`, `io`, `seq`, `stats`, `plot.style`, `calls.states`, `calls.records` | Shared tools: checkpoints, logging, constants, readers of FASTA, samples sheets and annotations, sequence routines, statistics |
| 1 | `routes.assembly`, `routes.library`, `routes.reads` | Evidence routes: each turns its inputs into pair call records. A route never imports another route |
| 2 | `calls` | The call table: allele states, the locus catalog, its writers (`loci.tsv`, VCF, GFF3, BED), `manifest.json` |
| 3 | `analyses.geneflow`, `analyses.subgenomes`, `analyses.source` | Analyses on top of the calls |
| 4 | `plot`, `sim` | Figures; simulators with known truth and the scorer |
| 5, 6 | `cli`, `__main__` | The command line |

A module imports lower ranks, its own unit, or another rank-0 unit (no cycles). A new module must be listed in
`UNITS` in `tests/test_layers.py`. Each unit's `__init__.py` docstring says what its modules do: read it first.

Other rules the tests enforce:

- No copied code: `tests/test_dupes.py` fails when two functions of 5 or more statements are the same once names and
  constants are ignored. Move shared code to a rank-0 unit.
- Importing `synltr.cli` and building its parser load no numpy, numba or msprime (`tests/test_cli.py`): commands
  import their code on demand.
- Errors: raise `SynltrError` (exit 1), `UsageError` (exit 2: a bad option or input file) or `AnalysisError`
  (exit 3), all in `base/log.py`. The command line prints one `error:` line.
- Stages are checkpointed (`base/checkpoint.py`): a stage's key hashes its inputs, parameters, the version and the
  code it runs, so a rerun reuses a stage only when none of these changed.
- `-t` is the total CPU budget across all processes.
- Names users do not know stay out of everything a user reads (`tests/test_names.py`): the README, the changelog,
  `docs/`, `examples/`, help text, messages and output files never name the four programs the package was built
  from (LTRpav, LTRmix, ltrkin, ltrharvest_phasing10), its earlier name (LTRtrace) or LTRquest (now LTR_peel).
  Docstrings, comments and tests may: there a name says which program a rule came from (LTRpav: assemblies; LTRmix:
  the screen; ltrkin: reads; ltrharvest_phasing10 or "phasing": subgenomes), and `legacy/` holds each of them.

## Coordinates

Name the convention at every format boundary and test each conversion on a tiny fixture.

- 1-based closed: GFF3, VCF, SAM, the output tables, LTRquest tables, library record names `chrom:start-end`.
- 0-based half-open: BED, bedGraph, PAF, and everything inside `routes.reads`.
- `docs/outputs.md` gives the convention of every output file.

## Golden pins

`tests/test_golden.py` compares four end-to-end runs byte for byte with `tests/golden/` (without the lines that name
the version or the run's folder). A failing pin means a reported number changed. If the change is intended, replace
the pinned files with the new run's and say what changed, and why, in `CHANGELOG.md`. Never edit a pin to make a
test pass.

## Releasing a change

- The version is written in five places: `pyproject.toml`, `src/synltr/__init__.py`, `CITATION.cff`, the README's
  citation line and `tests/test_cli.py`.
- User-visible changes go in `CHANGELOG.md`. Before a push, check the README, `docs/outputs.md` and
  `docs/benchmarks.md` against the code and cut what no longer holds.
- A benchmark number is quoted only with its source (tool, version, date, where the record is).

## Legacy

`legacy/` holds the hand-built synLTR pipeline as it was before the merge, laid out as it was so that its wrappers
still run in place, plus the original LTRpav, LTRmix and ltrkin as tarballs. [legacy/README.md](legacy/README.md)
indexes every file, says which approaches the package has and which exist only there, and where each benchmark's
record is; `legacy/changes_from_original_tools.md` lists what changed when the four programs became commands.

- Do not edit, reformat, lint or delete anything under `legacy/`. It is the record of how each method was first
  worked out. The one file to keep current is `legacy/README.md`.
- The package never imports from `legacy/`, and no test depends on it.
- To reuse a legacy approach: read its row in `legacy/README.md` and the script itself (its header and closing notes
  often hold the reasoning and the dead ends), then write it anew in the right layer of `src/synltr/` with tests on
  synthetic fixtures. Where the legacy script still runs, use its output on a small input as the reference to
  reproduce, as the 0.1.0 merge did for the four tools. Record the port in `CHANGELOG.md` and set the file's status
  in `legacy/README.md` to `ported`.
- The numbered scripts (`ltrharvest_phasing.py` to `ltrharvest_phasing10.py`, `ltrharvest.py` to `ltrharvest5.py`)
  are successive versions: the highest number is the current one, the lower ones show what was tried before.
