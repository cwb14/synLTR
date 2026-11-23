"""examples/quickstart.sh end to end: the README's command lines, each run on its toy from tests/builders.py. The
assembly toy's tables and the subgenomes toy's files equal their pins (the same toys that tests/test_golden.py runs
in-process at 2 threads; here the console script at 4)."""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import builders
import pytest

ROOT = Path(__file__).resolve().parents[1]
GOLDEN = ROOT / "tests" / "golden"
# a few files of each step, by the folder the quick start runs it in
OUTPUTS = {
    "assemblies/out": ("loci.tsv", "pav.vcf.gz", "pav.vcf.gz.tbi", "A1.pav.bed", "pairs/C1__A1.calls.tsv.gz",
                       "stats.tsv", "findings.txt", "tracts.A1.bed", "overview.png", "overview.pdf", "manifest.json"),
    "assemblies/score": ("C1__A1.score.tsv", "C1__A1.confusion.tsv"),
    "libraries/screen": ("pairs/R__D.calls.tsv.gz", "blocks.tsv", "posterior.R.tsv", "tracts.R.bed", "stats.tsv",
                         "findings.txt", "overview.png", "manifest.json"),
    "reads/out": ("reads.genotypes.tsv", "reads.vcf.gz", "reads.sources.tsv", "reads.copies.tsv",
                  "reads.locality.tsv", "manifest.json"),
    "polyploid/out": ("chromosomes.tsv", "dating.tsv", "report.pdf", "report.png", "manifest.json"),
}


@pytest.mark.slow
def test_quickstart(tmp_path):
    bindir = Path(sys.executable).parent  # this environment's programs first, as an activated environment has them
    path = f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}"
    missing = [p for p in ("synltr", "lastz", "minimap2", "samtools") if shutil.which(p, path=path) is None]
    if missing:
        pytest.skip(f"not on PATH: {', '.join(missing)}")
    pytest.importorskip("msprime")
    qs = tmp_path / "qs"
    run = subprocess.run(["bash", str(ROOT / "examples" / "quickstart.sh"), str(qs)], env={**os.environ, "PATH": path},
                         capture_output=True, text=True, check=False)
    assert run.returncode == 0, run.stdout[-3000:] + run.stderr[-3000:]
    for folder, names in OUTPUTS.items():
        assert not [n for n in names if not (qs / folder / n).is_file()], folder
    verdicts = [line for line in (qs / "libraries/screen/findings.txt").read_text().splitlines() if "Verdict, " in line]
    assert len(verdicts) == 2 and all("a signal of gene flow between R and D" in v for v in verdicts)
    for pin in sorted((GOLDEN / "geneflow_toy").iterdir()):
        assert (qs / "assemblies/out" / pin.name).read_bytes() == pin.read_bytes(), pin.name
    for pin in sorted((GOLDEN / "subgenomes_toy").iterdir()):
        expect = builders.unversioned_gff3(pin.read_bytes()) if pin.suffix == ".gff3" else pin.read_bytes()
        assert builders.subgenomes_text(qs / "polyploid/out" / pin.name) == expect, pin.name
