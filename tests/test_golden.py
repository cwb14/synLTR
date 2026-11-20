"""Pinned end-to-end outputs, compared byte for byte, without the lines that name the run's folder or the package's
version (so that neither fails a pin). The quickstart toy: four genomes simulated from a planted backbone, `geneflow`
with the trio A1,B1,C1, against tests/golden/geneflow_toy/. The reads route: `genotype` on the diploid read set of
builders.read_set, against tests/golden/genotype_reads/ (the VCF is pinned as text without its ##reference and
##source lines). The subgenomes toy: `subgenomes` on a small simulated allotetraploid, against
tests/golden/subgenomes_toy/ (its tables, findings.txt without its header line, which names the version, its GFF3
with `ltrtrace` alone in the source column, and report_data.json). The screen toy: `geneflow --screen` on LTRmix's
admixed trio (builders.trio, seed 3), against tests/golden/screen_toy/ (its tables, findings.txt pinned without its
header line, which names the version and the run's folder, and its pair files compared decompressed)."""
import gzip
import shutil
from pathlib import Path

import builders
import pytest

from ltrtrace.cli import main
from tests.markers import needs_lastz

GOLDEN = Path(__file__).parent / "golden" / "geneflow_toy"
READS_GOLDEN = Path(__file__).parent / "golden" / "genotype_reads"
SUBGENOMES_GOLDEN = Path(__file__).parent / "golden" / "subgenomes_toy"
SCREEN_GOLDEN = Path(__file__).parent / "golden" / "screen_toy"
TABLES = ("loci.tsv", "pav_matrix.tsv", "pav.unplaced.tsv", "qc.tsv", "posterior.A1.tsv", "tracts.A1.bed",
          "decoy_tracts.B1.bed", "inversions.tsv")


@needs_lastz
@pytest.mark.slow
def test_geneflow_toy(tmp_path):
    pytest.importorskip("msprime")
    assert builders.write_backbone(tmp_path / "toy.fa", tmp_path / "toy.gff3") == 40
    assert main(["sim", "pop", "--fasta", str(tmp_path / "toy.fa"), "--annotation", str(tmp_path / "toy.gff3"),
                 "--region", "s1:1-300000", "--region", "s2:1-300000", "--region", "s3:1-300000",
                 "--per-taxon", "1", "--f", "0.2", "--ta", "2000", "--seed", "2", "--threads", "2",
                 "--outdir", str(tmp_path / "sim")]) == 0
    run = tmp_path / "run"
    assert main(["geneflow", str(tmp_path / "sim" / "samples.tsv"), "-o", str(run), "-t", "2",
                 "--trio", "A1,B1,C1"]) == 0
    for name in TABLES:
        assert (run / name).read_bytes() == (GOLDEN / name).read_bytes(), name


@pytest.mark.slow
def test_genotype_reads(tmp_path):
    if shutil.which("minimap2") is None or shutil.which("samtools") is None:
        pytest.skip("minimap2 and samtools (to align the simulated reads) not on PATH")
    rs = builders.read_set(tmp_path / "in")
    sheet = builders.reads_sheet(tmp_path / "samples.tsv", rs.ref, rs.gff, {"full": rs.bam})
    out = tmp_path / "out"
    assert main(["genotype", str(sheet), "-o", str(out), "-t", "2"]) == 0
    for name in ("reads.genotypes.tsv", "reads.unplaced.tsv"):
        assert (out / name).read_bytes() == (READS_GOLDEN / name).read_bytes(), name
    with gzip.open(out / "reads.vcf.gz", "rt") as fh:
        vcf = "".join(ln for ln in fh if not ln.startswith(("##reference=", "##source=")))
    assert vcf.encode() == (READS_GOLDEN / "reads.vcf").read_bytes()


def subgenomes_tables(run: Path) -> list[str]:
    """The files of a subgenomes run that the pin holds: the tables, findings.txt and report_data.json (not the
    figure, summary.json or the work folder, which name the run's folder, the version or timings)."""
    return sorted(p.name for p in run.iterdir() if p.is_file() and (
        p.suffix in (".tsv", ".bed", ".gff3", ".bedGraph") or p.name in ("findings.txt", "report_data.json")))


@pytest.mark.slow
def test_subgenomes_toy(tmp_path):
    lib, _ = builders.simulate(builders.SUBGENOMES_TOY, tmp_path / "in", tmp_path / "truth", seed=1)
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text("label\tlibrary\tannotation\tfasta\n"
                     f"toy\t{lib / 'ltr.fa'}\t{lib / 'annotation.gff3'}\t{lib / 'genome.fai'}\n")
    run = tmp_path / "run"
    assert main(["subgenomes", str(sheet), "-o", str(run), "-t", "2"]) == 0
    pinned = sorted(p.name for p in SUBGENOMES_GOLDEN.glob("*"))
    assert pinned, f"nothing pinned in {SUBGENOMES_GOLDEN}"
    assert subgenomes_tables(run) == pinned
    for name in pinned:
        pin = (SUBGENOMES_GOLDEN / name).read_bytes()
        assert builders.subgenomes_text(run / name) == (builders.unversioned_gff3(pin) if name.endswith(".gff3")
                                                          else pin), name


def screen_files(run: Path) -> list[str]:
    """The files of a library screen that the pin holds: the run folder's tables and findings.txt, and its pair files
    (pairs/*.tsv.gz); not the figure (overview.png, overview.pdf), the done markers or the work folder."""
    top = [p.name for p in run.iterdir() if p.is_file() and (p.suffix in (".tsv", ".bed") or p.name == "findings.txt")]
    return sorted(top + [f"pairs/{p.name}" for p in (run / "pairs").glob("*.tsv.gz")])


def screen_text(path: Path) -> bytes:
    """A screen file as pinned: a gzip file decompressed, findings.txt without its header line (which names the
    version and the run's folder); every other byte as written, so floats are compared exactly as written."""
    data = gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes()
    return data.split(b"\n", 1)[1] if path.name == "findings.txt" else data


@pytest.mark.slow
def test_screen_toy(tmp_path):
    """`ltrtrace geneflow --screen` on LTRmix's admixed trio (builders.trio, seed 3, as the screen tests build it),
    from an empty folder: the run writes exactly the pinned files and each equals its pin."""
    sc = builders.trio(tmp_path / "in", seed=3, admix=True)
    run = tmp_path / "run"
    assert main(["geneflow", "--screen", str(sc.sheet), "-o", str(run), "-t", "2"]) == 0
    pinned = sorted(p.relative_to(SCREEN_GOLDEN).as_posix() for p in SCREEN_GOLDEN.rglob("*") if p.is_file())
    assert pinned, f"nothing pinned in {SCREEN_GOLDEN}"
    assert screen_files(run) == pinned
    for name in pinned:
        pin = SCREEN_GOLDEN / name
        assert screen_text(run / name) == (gzip.decompress(pin.read_bytes()) if pin.suffix == ".gz"
                                           else pin.read_bytes()), name
