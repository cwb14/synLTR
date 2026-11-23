"""The reads route end to end (`synltr genotype` on a sheet of one reference and its read sets): ltrkin's
reads.genotypes.tsv and reads.vcf.gz for the run, every read set a sample, the elements it reads, its input checks,
the evidence checkpoint per read set, and the same outputs at any thread count."""
import csv
import gzip
import hashlib
import json
import os
import random
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

import builders
import pysam
import pytest

from synltr import cli
from synltr.io.seqio import Fasta
from synltr.routes.reads import genotype
from synltr.routes.reads.genotype import boundary_tolerance, relabel_other
from synltr.routes.reads.reads import Site, Target
from synltr.routes.reads.run import load_targets

OUTPUTS = ("reads.genotypes.tsv", "reads.vcf.gz", "reads.vcf.gz.tbi", "reads.unplaced.tsv")


def _run(sheet: Path, out: Path, threads: int = 2, *extra: str) -> Path:
    assert cli.main(["genotype", str(sheet), "-o", str(out), "-t", str(threads), *extra]) == 0
    return out


def _rows(path: Path, sample: str) -> dict[str, dict]:
    """The rows of read set `sample` in the table at `path`, by element."""
    with open(path) as fh:
        return {r["element_id"]: r for r in csv.DictReader(fh, delimiter="\t") if r["sample"] == sample}


def _calls(path: Path, sample: str = "full") -> dict[str, tuple]:
    return {k: (r["genotype"], r["gq"], r["n_F"], r["n_E"], r["n_S"]) for k, r in _rows(path, sample).items()}


def _genotype(data, tmp: Path, name: str, *, gff=None, bam=None, ref=None, label="full") -> dict[str, tuple]:
    """ltrkin's `genotype -r REF -a GFF -o <tmp>/<name> BAM`: the read set `label` (data.bam by default) genotyped
    into <tmp>/<name>; its calls by element."""
    sheet = builders.reads_sheet(tmp / f"{name}.samples.tsv", ref or data.ref, gff or data.gff,
                                 {label: bam or data.bam})
    out = _run(sheet, tmp / name)
    return _calls(out / "reads.genotypes.tsv", label)


def _files(out: Path) -> dict[str, bytes]:
    return {name: (out / name).read_bytes() for name in OUTPUTS}


@pytest.fixture(scope="module")
def two(data, tmp_path_factory):
    """Both read sets of the diploid (every record, primary records only) genotyped in one run, 4 threads."""
    d = tmp_path_factory.mktemp("two")
    sheet = builders.reads_sheet(d / "samples.tsv", data.ref, data.gff, {"full": data.bam, "primary": data.primary})
    return sheet, _run(sheet, d / "out", 4)


def test_run_manifest_names_the_reference(data, two):
    """manifest.json records the sheet's rows (the reference with its FASTA and annotation first), the parameters,
    the tools and the reference FASTA's checksum, which every read set's evidence is keyed on (a BAM is keyed by its
    size, modification time and header, never hashed); `synltr source` reads it to find the reference."""
    sheet, out = two
    doc = json.loads((out / "manifest.json").read_text())
    ref = [s for s in doc["samples"] if s["role"] == "reference"]
    assert len(ref) == 1 and ref[0]["fasta"] == str(data.ref) and ref[0]["annotation"] == str(data.gff)
    assert [s["label"] for s in doc["samples"] if s["role"] != "reference"] == ["full", "primary"]
    assert doc["params"]["ploidy"] == 2 and doc["command"][:3] == ["synltr", "genotype", str(sheet)]
    assert doc["params"]["min_mapq"] == 20 and doc["route"] == "reads" and list(doc["code"]) == ["calls.genotype"]
    assert [s["reads"] for s in doc["samples"]] == [None, str(data.bam), str(data.primary)]
    assert ref[0]["checksums"] == {"fasta": hashlib.sha256(Path(data.ref).read_bytes()).hexdigest()}
    assert [s["checksums"] for s in doc["samples"] if s["role"] != "reference"] == [{}, {}]
    assert ref[0]["skipped"] == {}


def test_failed_rerun_leaves_no_manifest(two, tmp_path, monkeypatch):
    """A read genotyping rerun in a finished run's folder that fails once its checks have passed leaves no
    manifest.json describing the earlier run; a usage error (--pair on read sets) leaves the record as it was."""
    from synltr.calls import genotype as command

    sheet, done = two
    out = tmp_path / "out"
    shutil.copytree(done, out)
    before = (out / "manifest.json").read_bytes()
    assert cli.main(["genotype", str(sheet), "-o", str(out), "--pair", "full,primary"]) == 2
    assert (out / "manifest.json").read_bytes() == before

    def boom(*args, **kwargs):
        raise RuntimeError("a planted failure")

    monkeypatch.setattr(command, "genotype_reads", boom)
    assert cli.main(["genotype", str(sheet), "-o", str(out), "-t", "2"]) == 1
    assert not (out / "manifest.json").exists()


# ---- the elements the route reads

def test_gff_ltrquest(data):
    els, _ = load_targets(str(data.gff), Fasta(data.ref), "ref")
    ids = {e.id for e in els}
    assert {"E1_het_empty", "E5_nested_outer", "E5_nested_outer_inner"} <= ids
    outer = next(e for e in els if e.id == "E5_nested_outer")
    inner = next(e for e in els if e.id == "E5_nested_outer_inner")
    assert outer.n_inner == 1 and inner.n_outer == 1
    assert outer.lltr[1] - outer.lltr[0] == 600
    assert outer.inner == ((inner.start, inner.end),)


def _fasta(path: Path, n: int) -> Fasta:
    return Fasta(builders.write_fasta(path, {"c1": builders.rand_seq(n, random.Random(1))}))


def test_gff_edta_style(tmp_path):
    g = tmp_path / "edta.gff3"
    g.write_text(
        "##gff-version 3\n"
        "c1\tEDTA\trepeat_region\t100\t5210\t.\t+\t.\tID=repeat_region_1\n"
        "c1\tEDTA\ttarget_site_duplication\t100\t104\t.\t+\t.\tID=tsd_1;Parent=repeat_region_1\n"
        "c1\tEDTA\tCopia_LTR_retrotransposon\t105\t5205\t.\t+\t.\tID=LTRRT_1;Parent=repeat_region_1;"
        "Name=TE_00001;Classification=LTR/Copia;tsd=ACGTA\n"
        "c1\tEDTA\tlong_terminal_repeat\t105\t504\t.\t+\t.\tID=lLTR_1;Parent=repeat_region_1\n"
        "c1\tEDTA\tlong_terminal_repeat\t4806\t5205\t.\t+\t.\tID=rLTR_1;Parent=repeat_region_1\n"
        "c1\tEDTA\ttarget_site_duplication\t5206\t5210\t.\t+\t.\tID=tsd_2;Parent=repeat_region_1\n")
    (e,), _ = load_targets(str(g), _fasta(tmp_path / "c1.fa", 6000), "ref")
    assert (e.start, e.end, e.lltr, e.rltr) == (104, 5205, (104, 504), (4805, 5205))
    assert e.family == "TE_00001" and e.clade == "Copia" and e.tsd == "ACGTA" and e.id == "LTRRT_1"


def test_gff_dante_ltr_style(tmp_path):
    g = tmp_path / "dante.gff3"
    cls = "Class_I|LTR|Ty3/gypsy|non-chromovirus|OTA|Athila"
    g.write_text(
        "##gff-version 3\n"
        f"c1\tdante_ltr\ttransposable_element\t1001\t9000\t.\t+\t.\tID=TE_1;Final_Classification={cls};Name={cls}\n"
        "c1\tdante_ltr\tlong_terminal_repeat\t1001\t1500\t.\t+\t.\tLTR=5LTR;Parent=TE_1;\n"
        "c1\tdante_ltr\tprotein_domain\t3001\t3500\t.\t+\t.\tName=RT;Parent=TE_1\n"
        "c1\tdante_ltr\tlong_terminal_repeat\t8501\t9000\t.\t+\t.\tLTR=3LTR;Parent=TE_1;\n"
        f"c1\tdante_ltr\ttransposable_element\t20001\t23000\t.\t+\t.\tID=TE_partial_2;Final_Classification={cls}\n")
    els, skipped = load_targets(str(g), _fasta(tmp_path / "c1.fa", 25000), "ref")
    assert len(els) == 1 and skipped["no_ltr_pair"] == 1  # the domain-only partial element is skipped
    (e,) = els
    assert (e.start, e.end, e.lltr, e.rltr) == (1000, 9000, (1000, 1500), (8500, 9000))
    assert e.clade == "Athila"


# ---- genotypes

@pytest.mark.parametrize("mode", ["full", "primary"])
def test_end_to_end(data, two, mode):
    _, out = two
    res = _rows(out / "reads.genotypes.tsv", mode)
    for eid, truth in data.truth.items():
        r = res[eid]
        if eid == "E5_nested_outer_inner":  # flanks are near-identical family-B sequence: a call may be withheld
            assert r["genotype"] in (truth, "N/N", "F/N") or int(r["gq"]) < 20, (eid, r["genotype"])
            continue
        if eid in ("E5_nested_outer", "E9_internal_del"):
            # hap2 carries the element with another internal length (no inner insertion / an internal deletion):
            # its spanning reads carry both element ends, so they count as full and the site is flagged
            assert "full_other_length" in r["flags"] and int(r["n_F_other"]) > 0, (mode, eid, r["flags"])
            # E9's other-length reads were the only hap2 evidence (v0.1 called F/N); E5's nested outer element has
            # junction reads on both haplotypes and its F/F vs F/N rests on depth (identical haplotypes, 15x each)
            assert r["genotype"] == "F/F" if eid == "E9_internal_del" else r["genotype"] in ("F/F", "F/N")
            continue
        assert r["genotype"] == truth, (mode, eid, r["genotype"], r["gq"], r["n_F"], r["n_E"], r["n_S"], r["flags"])
        # F/F vs F/N rests on read depth alone here (identical haplotypes, no phasing): 15x per haplotype supports
        # only a weak GQ; empty and solo alleles are seen directly and must be confident
        assert int(r["gq"]) >= (3 if truth in ("F/F", "F/N") else 20), (mode, eid, r["gq"])
    assert (out / "reads.vcf.gz").exists()


def test_two_read_sets_share_one_table_and_one_vcf(data, two):
    """ltrkin's run of several BAMs: one genotypes table, each read set's rows in turn (sheet order, elements in
    order), and one VCF with a sample column per read set, every element either holds; one evidence checkpoint per
    read set."""
    _, out = two
    with open(out / "reads.genotypes.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert [r["sample"] for r in rows] == ["full"] * 9 + ["primary"] * 9
    els = [r["element_id"] for r in rows[:9]]
    assert els == [r["element_id"] for r in rows[9:]] and set(els) == set(data.truth)
    assert [(r["chrom"], int(r["start"])) for r in rows[:9]] == sorted((r["chrom"], int(r["start"])) for r in rows[:9])
    with gzip.open(out / "reads.vcf.gz", "rt") as fh:
        lines = [ln.rstrip("\n").split("\t") for ln in fh if not ln.startswith("##")]
    assert lines[0][9:] == ["full", "primary"] and len(lines) == 1 + 9
    by = {(r["sample"], r["element_id"]): r for r in rows}
    for rec in lines[1:]:
        for s, cell in zip(("full", "primary"), rec[9:]):
            r = by[(s, rec[2])]
            assert cell.split(":")[:3] == [r["vcf_gt"], r["gq"], f"{r['n_F']},{r['n_E']},{r['n_S']}"]
    assert sorted(p.name for p in out.glob("reads.*.evidence.json.gz")) == [
        "reads.full.evidence.json.gz", "reads.primary.evidence.json.gz"]
    assert not list(out.glob("full.*")) and not list(out.glob("primary.*"))
    assert (out / "reads.unplaced.tsv").read_text() == "lid\tref_state\treason\n"  # no element at base 1 here


@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools (to write the CRAM) not on PATH")
def test_cram_genotypes_as_its_bam(data, two, tmp_path):
    """A read set given as CRAM (indexed, decoded with the reference row's FASTA, the one it was written against) is
    genotyped exactly as the same reads in BAM: the same rows of reads.genotypes.tsv."""
    cram = tmp_path / "full.cram"
    subprocess.run(["samtools", "view", "-C", "-T", str(data.ref), "-o", str(cram), str(data.bam)], check=True)
    subprocess.run(["samtools", "index", str(cram)], check=True)
    sheet = builders.reads_sheet(tmp_path / "samples.tsv", data.ref, data.gff, {"full": cram})
    out = _run(sheet, tmp_path / "out", 4)
    rows = _rows(out / "reads.genotypes.tsv", "full")
    assert cram.read_bytes()[:4] == b"CRAM" and len(rows) == len(data.truth)
    assert rows == _rows(two[1] / "reads.genotypes.tsv", "full")


def test_split_rescue_matters(two):
    """In the primary-only BAM some empty-site reads are single pieces; they must still be counted."""
    _, out = two
    full, prim = _rows(out / "reads.genotypes.tsv", "full"), _rows(out / "reads.genotypes.tsv", "primary")
    e_full, e_prim = int(full["E1_het_empty"]["n_E"]), int(prim["E1_het_empty"]["n_E"])
    assert e_prim >= 0.8 * e_full


def test_reproducible_across_hash_seeds(data, tmp_path):
    """Same inputs, different Python hash seeds -> identical genotype tables (set ordering must not leak)."""
    sheet = builders.reads_sheet(tmp_path / "samples.tsv", data.ref, data.gff, {"full": data.bam})
    outs = []
    for seed in ("1", "2"):
        out = tmp_path / f"h{seed}"
        env = dict(os.environ, PYTHONHASHSEED=seed)
        subprocess.run([sys.executable, "-m", "synltr", "genotype", str(sheet), "-o", str(out), "-t", "2"],
                       check=True, env=env, capture_output=True)
        outs.append((out / "reads.genotypes.tsv").read_text())
    assert outs[0] == outs[1]


def _out_of_order(monkeypatch) -> None:
    """One element per evidence job, and the job of element 0 a second late: forked workers then finish the jobs in
    another order than the one they were given."""
    monkeypatch.setattr(genotype, "EVIDENCE_CHUNK", 1)
    real = genotype._scan_chunk

    def late_first(chunk):
        if chunk[0][0] == 0:
            time.sleep(1.0)
        return real(chunk)

    monkeypatch.setattr(genotype, "_scan_chunk", late_first)


def test_tie_stable_across_threads(two, tmp_path, monkeypatch):
    """With forked workers finishing the evidence jobs out of order (9 one-element jobs, element 0's last), the
    evidence still comes back in element order: every output and both evidence checkpoints are the same at -t 1 and
    -t 4, and with the default jobs."""
    sheet, out_default = two
    _out_of_order(monkeypatch)
    out1 = _run(sheet, tmp_path / "t1", 1)
    out4 = _run(sheet, tmp_path / "t4", 4)
    assert _files(out1) == _files(out4) == _files(out_default)
    for label in ("full", "primary"):
        cache = f"reads.{label}.evidence.json.gz"
        data = [gzip.decompress((d / cache).read_bytes()) for d in (out1, out4, out_default)]
        assert data[0] == data[1] == data[2]
        keys = list(json.loads(data[1])["evidence"])
        assert keys == sorted(keys, key=int) and len(keys) == 9


def _renamed(data, d: Path, names: dict[str, str]) -> tuple[Path, Path, Path]:
    """The diploid's reference, annotation and every-record BAM with its sequences renamed by `names`, kept in
    their order (the BAM's records point at the same header entries)."""
    ref, gff, bam = d / "renamed.fa", d / "renamed.gff3", d / "renamed.bam"
    with pysam.FastaFile(str(data.ref)) as fa:
        builders.write_fasta(ref, {names[c]: fa.fetch(c) for c in fa.references})
    gff.write_text("##gff-version 3\n" + "".join(names[ln.split("\t", 1)[0]] + "\t" + ln.split("\t", 1)[1]
                                                  for ln in _gff_lines(data.gff)))
    with pysam.AlignmentFile(str(data.bam)) as inp:
        header = inp.header.to_dict()
        for sq in header["SQ"]:
            sq["SN"] = names[sq["SN"]]
        with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
            for x in inp:
                out.write(x)
    pysam.index(str(bam))
    return ref, gff, bam


def test_table_in_text_order_and_vcf_in_fai_order(data, tmp_path):
    """Sequence names whose text order (chr10, chr2, chr3) differs from the FASTA's (chr2, chr10, chr3): the table
    keeps ltrkin's element order (names as text), the VCF's contig lines and records follow the .fai."""
    ref, gff, bam = _renamed(data, tmp_path, {"chr1": "chr2", "chr2": "chr10", "chr3": "chr3"})
    out = _run(builders.reads_sheet(tmp_path / "s.tsv", ref, gff, {"full": bam}), tmp_path / "out")
    with open(out / "reads.genotypes.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert [r["chrom"] for r in rows] == ["chr10"] * 2 + ["chr2"] * 6 + ["chr3"]
    assert [r["element_id"] for r in rows[:2]] == ["E7_hom", "E8_decoy_hom"]
    with gzip.open(out / "reads.vcf.gz", "rt") as fh:
        lines = fh.read().splitlines()
    contigs = [ln.removeprefix("##contig=<ID=").split(",")[0] for ln in lines if ln.startswith("##contig")]
    assert contigs == ["chr2", "chr10", "chr3"]
    assert [ln.split("\t")[0] for ln in lines if not ln.startswith("#")] == ["chr2"] * 6 + ["chr10"] * 2 + ["chr3"]


def test_duplicate_sample_names_get_file_names(data, tmp_path):
    """Two BAMs with the same read-group SM (read sets of one plant) stay two samples, named by their labels."""
    b1, b2 = tmp_path / "first.bam", tmp_path / "second.bam"
    for src, dst in ((data.bam, b1), (data.primary, b2)):
        pysam.samtools.addreplacerg("-r", "ID:rg1", "-r", "SM:plant", "-o", str(dst), str(src), catch_stdout=False)
        pysam.index(str(dst))
    sheet = builders.reads_sheet(tmp_path / "samples.tsv", data.ref, data.gff, {"first": b1, "second": b2})
    out = _run(sheet, tmp_path / "multi")
    with open(out / "reads.genotypes.tsv") as fh:
        samples = {r["sample"] for r in csv.DictReader(fh, delimiter="\t")}
    assert samples == {"first", "second"}
    with gzip.open(out / "reads.vcf.gz", "rt") as fh:
        hdr = next(ln for ln in fh if ln.startswith("#CHROM")).rstrip("\n").split("\t")
    assert len(set(hdr[9:])) == 2


def test_boundary_tolerance_learned():
    """Exact ends -> no widening; offsets of ~40 bp at TSD-less elements -> about 3 x 40."""
    sites, ev = {}, {}
    for i in range(30):
        el = Target(f"e{i}", "c", 10000, 20000, "+", (10000, 10500), (19500, 20000))
        sites[i] = Site(el, 5)
        exp = sites[i].delta["E"]
        ev[i] = {"E": 0, "S": 0, "O": 3, "other_delta": [exp + 40 * (-1) ** i] * 3, "es_delta": []}
    tol, n = boundary_tolerance(ev, sites)
    assert n == 30 and 100 <= tol <= 140
    assert relabel_other(ev, sites, tol) == 90 and ev[0]["E"] == 3 and ev[0]["O"] == 0
    exact = {i: {"E": 3, "S": 0, "O": 0, "other_delta": [], "es_delta": [("E", sites[i].delta["E"])] * 3}
             for i in range(30)}
    assert boundary_tolerance(exact, sites)[0] == 0


# ---- inputs and the evidence checkpoint

def _gff_lines(path):
    with open(path) as fh:
        return [ln for ln in fh if ln.strip() and not ln.startswith("#")]


def test_evidence_cache_follows_the_annotation(data, tmp_path):
    """Same number of elements, another annotation: a rerun in the same outdir must not reuse the evidence."""
    _genotype(data, tmp_path, "run")
    lines = _gff_lines(data.gff)
    drop = [ln for ln in lines if "E8_decoy_hom" in ln]
    assert drop
    keep = [ln for ln in lines if "E8_decoy_hom" not in ln]
    add = ["chr1\tt\tLTR_retrotransposon\t140001\t145000\t.\t+\t.\tID=X_new;family=famA\n",
           "chr1\tt\tlong_terminal_repeat\t140001\t140400\t.\t+\t.\tParent=X_new\n",
           "chr1\tt\tlong_terminal_repeat\t144601\t145000\t.\t+\t.\tParent=X_new\n"]
    gff2 = tmp_path / "other.gff3"
    gff2.write_text("##gff-version 3\n" + "".join(keep + add))
    rerun = _genotype(data, tmp_path, "run", gff=gff2)
    fresh = _genotype(data, tmp_path, "fresh", gff=gff2)
    assert rerun == fresh


def test_truncated_cache_is_rescanned(data, tmp_path):
    first = _genotype(data, tmp_path, "t")
    cache = tmp_path / "t" / "reads.full.evidence.json.gz"
    raw = cache.read_bytes()
    cache.write_bytes(raw[: len(raw) // 2])
    assert _genotype(data, tmp_path, "t") == first
    assert cache.read_bytes() == raw


def test_unreadable_cache_with_its_size_is_rescanned(data, tmp_path, capsys):
    """A checkpoint whose bytes were damaged in place (same size): one warning, then a fresh scan."""
    first = _genotype(data, tmp_path, "u")
    cache = tmp_path / "u" / "reads.full.evidence.json.gz"
    raw = cache.read_bytes()
    cache.write_bytes(b"\0" * len(raw))
    capsys.readouterr()
    assert _genotype(data, tmp_path, "u") == first
    assert cache.read_bytes() == raw and "unreadable" in capsys.readouterr().err


def test_duplicate_element_ids_are_rejected(data, tmp_path, capsys):
    lines = _gff_lines(data.gff)
    e7 = [ln.replace("E7_hom", "E1_het_empty") for ln in lines if "E7_hom" in ln]  # another chromosome, same ID
    g = tmp_path / "dup.gff3"
    g.write_text("##gff-version 3\n" + "".join([ln for ln in lines if "E7_hom" not in ln] + e7))
    sheet = builders.reads_sheet(tmp_path / "s.tsv", data.ref, g, {"full": data.bam})
    assert cli.main(["genotype", str(sheet), "-o", str(tmp_path / "d")]) == 2
    assert "more than once" in capsys.readouterr().err


def test_missing_reference_sequence_is_skipped_with_ltrkins_warning(data, tmp_path, capsys):
    """An annotation sequence the reference lacks (an unplaced scaffold, say): its elements are skipped with
    ltrkin's warning line, word for word, the others genotyped as without it (the shared loader's "skip" rule; the
    default fails)."""
    g = tmp_path / "extra.gff3"
    g.write_text("##gff-version 3\n" + "".join(_gff_lines(data.gff)) +
                 "chrZ\tt\tLTR_retrotransposon\t1001\t6000\t.\t+\t.\tID=Z_1\n"
                 "chrZ\tt\tlong_terminal_repeat\t1001\t1400\t.\t+\t.\tParent=Z_1\n"
                 "chrZ\tt\tlong_terminal_repeat\t5601\t6000\t.\t+\t.\tParent=Z_1\n")
    capsys.readouterr()
    calls = _genotype(data, tmp_path, "z", gff=g)
    err = capsys.readouterr().err
    assert "warning: 1 annotation sequence(s) not in the reference, skipped: ['chrZ']" in err.splitlines()
    assert "Z_1" not in calls and calls == _genotype(data, tmp_path, "plain")


def test_bam_aligned_to_another_reference_is_rejected(data, tmp_path, capsys):
    fa = pysam.FastaFile(str(data.ref))
    ref2 = tmp_path / "ref_v2.fa"
    with open(ref2, "w") as fh:
        for c in fa.references:
            s = fa.fetch(c) + ("ACGT" * 250 if c == "chr1" else "")  # chr1 1 kb longer
            fh.write(f">{c}\n" + "\n".join(s[i:i + 80] for i in range(0, len(s), 80)) + "\n")
    sheet = builders.reads_sheet(tmp_path / "s.tsv", ref2, data.gff, {"full": data.bam})
    assert cli.main(["genotype", str(sheet), "-o", str(tmp_path / "r")]) == 2
    err = capsys.readouterr().err
    assert "another reference" in err and len(err.splitlines()) == 1


def test_element_past_the_sequence_end_is_rejected(data, tmp_path, capsys):
    g = tmp_path / "past.gff3"
    g.write_text("##gff-version 3\n" + "".join(_gff_lines(data.gff)) +
                 "chr3\tt\tLTR_retrotransposon\t39001\t47800\t.\t+\t.\tID=X_past\n"
                 "chr3\tt\tlong_terminal_repeat\t39001\t39400\t.\t+\t.\tParent=X_past\n"
                 "chr3\tt\tlong_terminal_repeat\t47401\t47800\t.\t+\t.\tParent=X_past\n")
    sheet = builders.reads_sheet(tmp_path / "s.tsv", data.ref, g, {"full": data.bam})
    assert cli.main(["genotype", str(sheet), "-o", str(tmp_path / "p")]) == 2
    assert "X_past" in capsys.readouterr().err


def test_duplicate_flagged_reads_are_ignored(data, tmp_path):
    """Reads flagged as duplicates must count nowhere: the calls equal those on the BAM without them."""
    flagged, removed = tmp_path / "flagged.bam", tmp_path / "removed.bam"
    with pysam.AlignmentFile(str(data.bam)) as inp, pysam.AlignmentFile(str(flagged), "wb", template=inp) as f, \
            pysam.AlignmentFile(str(removed), "wb", template=inp) as r:
        for x in inp:
            if sum(map(ord, x.query_name)) % 5 < 2:  # ~40% of reads, all their records
                x.flag |= 0x400
            else:
                r.write(x)
            f.write(x)
    for b in (flagged, removed):
        pysam.index(str(b))
    assert _genotype(data, tmp_path, "f", bam=flagged) == _genotype(data, tmp_path, "r", bam=removed)


def test_vcf_values_are_encoded(data, tmp_path):
    """An ID with a character that has a meaning in VCF is percent-encoded there, so the VCF stays valid."""
    g = tmp_path / "enc.gff3"
    g.write_text("##gff-version 3\n" + "".join(ln.replace("E3_hom", "E3%3Bhom") for ln in _gff_lines(data.gff)))
    calls = _genotype(data, tmp_path, "v", gff=g)
    assert "E3;hom" in calls  # the table keeps the decoded ID
    with pysam.VariantFile(str(tmp_path / "v" / "reads.vcf.gz")) as vf:
        ids = {rec.id for rec in vf}
    assert "E3%3Bhom" in ids
    with gzip.open(tmp_path / "v" / "reads.vcf.gz", "rt") as fh:
        assert all(";hom" not in ln for ln in fh)


# a genotype run whose second read set (primary) stalls in its evidence scan: argv = the command line
_SLOW_SECOND_SET = r"""
import sys, time
from synltr import cli
from synltr.routes.reads import genotype

real = genotype.scan


def scan(label, *args, **kw):
    if label == "primary":
        time.sleep(120)
    return real(label, *args, **kw)


genotype.scan = scan
sys.exit(cli.main(sys.argv[1:]))
"""


def test_evidence_resume(two, tmp_path, capsys):
    """A two-read-set run killed (SIGKILL) after the first set's evidence checkpoint, while the second is scanned:
    the rerun reuses the first set's evidence and writes the same files, byte for byte, as an uninterrupted run."""
    sheet, ref_out = two
    out = tmp_path / "out"
    cache = out / "reads.full.evidence.json.gz"
    done = cache.with_name(cache.name + ".done.json")
    script = tmp_path / "slow.py"
    script.write_text(_SLOW_SECOND_SET)
    proc = subprocess.Popen([sys.executable, str(script), "genotype", str(sheet), "-o", str(out), "-t", "2"],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        t0 = time.monotonic()
        while not done.exists():
            assert proc.poll() is None and time.monotonic() - t0 < 120
            time.sleep(0.05)
        time.sleep(1)  # the first set's table, then the second set's scan (stalled)
        os.kill(proc.pid, signal.SIGKILL)
        proc.wait()
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait()
    assert not (out / "reads.primary.evidence.json.gz").exists() and not (out / "reads.vcf.gz").exists()
    before = (cache.stat().st_mtime_ns, done.read_bytes())
    capsys.readouterr()
    _run(sheet, out, 2)
    err = capsys.readouterr().err
    assert "full: evidence done, skip" in err and "primary: evidence done" not in err
    assert (cache.stat().st_mtime_ns, done.read_bytes()) == before  # reused, not rewritten
    assert _files(out) == _files(ref_out)
