"""Outputs of `synltr sim nucleus`: an LTRquest-style library (LTR-RT FASTA, K2P table, GFF3, .fai) in one folder
and the truth tables in another.

Library: ltr.fa (records `chrom:start-end#class`, start-end 1-based closed, sequence in element orientation),
k2p.tsv (LTRquest's table; LTR offsets are 1-based positions along the copy's genome span, counted from its start,
the genome-left LTR first on either strand, as LTRquest writes them), annotation.gff3 (seqids and Name escaped as
GFF3 requires), genome.fai. Truth: copies.tsv, exchanges.tsv, chromosomes.tsv, orthologs.tsv, paralogs.tsv,
history.json. Times are years before present. In copies.tsv `root` identifies the original insertion a copy descends
from (copies inherited across a split or moved by an exchange share it), and `source` is the `root` of the copy this
insertion was transposed from, or "." when its source was a master.
"""
from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

from synltr.base.constants import GFF3_SOURCE, SIM_NUCLEUS_K2P_RATE, SIM_NUCLEUS_SEQID
from synltr.io.annotation.gff3 import escape_seqid, escape_value
from synltr.sim.nucleus.model import DUP, Seqs, Sim, k2p, merge_defaults

# ---------------------------------------------------------------------------------------------- outputs
BASES = np.frombuffer(b"ACGT", np.uint8)
TSV_COLS = ["seq_id", "seq_len", "status", "ltr5_start", "ltr5_end", "ltr3_start", "ltr3_end", "ltr5_len", "ltr3_len",
            "flank5_len", "flank3_len", "aln_len", "n_sites", "n_ts", "n_tv", "n_gapcols", "identity", "p_dist", "k2p",
            "k2p_se", "bitscore", "flank_margin_bits", "cigar", "motif", "k2p_time", "orientation", "tsd",
            "tsd_offset", "tsd_input", "strand", "family", "domains", "nest_status"]


def _fmt(x):
    return "NA" if isinstance(x, float) and math.isnan(x) else (f"{x:.6g}" if isinstance(x, float) else str(x))


def _mrca_split(sim, a, b):
    def path(i):
        out = []
        while i >= 0:
            out.append(i)
            i = sim.lin[i].parent
        return out
    pa = path(a)
    for i in path(b):
        if i in pa:
            return sim.lin[i].t_to
    return float("nan")


def _unique_coordinates(sim, surv, tok):
    """Survivors with distinct output coordinates (an annotator never reports two elements at identical
    coordinates; a chance collision keeps the first) and the number dropped."""
    seen, keep = set(), []
    for c in surv:
        k = (tok[(sim.lineage[c], sim.chrom[c])], sim.pos[c], sim._elem_len(sim.fam[c]))
        if k not in seen:
            seen.add(k)
            keep.append(c)
    return keep, len(surv) - len(keep)


def input_family_labels(sim, surv):
    """Family labels as an annotator would give them: the true family, except a share of copies (more often the
    older ones: probability 2 * share * age rank / n) written as one-copy families."""
    n = len(surv)
    rank = np.argsort(np.argsort([sim.tins[c] for c in surv], kind="stable"), kind="stable")
    p = 2 * sim.sc["singleton_share"] * rank / max(n - 1, 1)
    single = np.random.default_rng([sim.seed, 7]).random(n) < p
    out, k = {}, 0
    for c, sg in zip(surv, single):
        if sg:
            k += 1
            out[c] = f"single{k:07d}"
        else:
            out[c] = f"fam{sim.fam[c]:04d}"
    return out


def write_outputs(sim, outdir, truthdir):
    outdir, truthdir = Path(outdir), Path(truthdir)
    outdir.mkdir(parents=True, exist_ok=True)
    truthdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng([sim.seed, 99])
    pairs = [(l, c) for l in sim.leaves for c in range(sim.n_chr)]
    tok = {p: SIM_NUCLEUS_SEQID.format(j + 1) for j, p in zip(rng.permutation(len(pairs)), pairs)}
    with open(outdir / "genome.fai", "w") as f:
        f.writelines(f"{tok[p]}\t{sim.chr_len}\t0\t0\t0\n" for p in sorted(pairs, key=tok.get))
    with open(truthdir / "chromosomes.tsv", "w") as f:
        f.write("token\tsubgenome\thomoeolog\tlength\n")
        f.writelines(f"{tok[p]}\t{sim.lin[p[0]].name}\t{p[1]}\t{sim.chr_len}\n" for p in sorted(pairs, key=tok.get))
    q = Seqs(sim)
    surv = sorted(sim.survivors(), key=lambda c: (tok[(sim.lineage[c], sim.chrom[c])], sim.pos[c]))
    surv, n_dup = _unique_coordinates(sim, surv, tok)
    label = input_family_labels(sim, surv)
    fam_size = defaultdict(int)
    for c in surv:
        fam_size[label[c]] += 1
    cid = {}
    leafname = {i: sim.lin[i].name for i in sim.leaves}
    with open(outdir / "ltr.fa", "w") as ffa, open(outdir / "k2p.tsv", "w") as ftsv, \
            open(outdir / "annotation.gff3", "w") as fgff, open(truthdir / "copies.tsv", "w") as ftr:
        ftsv.write("#" + "\t".join(TSV_COLS) + "\n")
        fgff.write("##gff-version 3\n")
        fgff.writelines(f"##sequence-region\t{escape_seqid(tok[p])}\t1\t{sim.chr_len}\n"
                        for p in sorted(pairs, key=tok.get))
        ftr.write("copy_id\tchrom\tsubgenome\thomoeolog\tfamily\tsuperfamily\tclade\tt_insert\tt_birth\torigin\tmode\t"
                  "post_merger\texchange_id\troot\tn_sites\tn_ts\tn_tv\tk2p\tsource\n")
        for n, c in enumerate(surv, 1):
            f = sim.fam[c]
            st = q.state(c, 0.0)
            codes = np.concatenate([st["ltr5"], st["int"], st["ltr3"]])   # element orientation, as LTRquest writes
            L, lt = codes.size, int(sim.ltr_len[f])
            left = right = lt    # LTRquest's offsets run along the genome: the left LTR first, on either strand
            chrom = tok[(sim.lineage[c], sim.chrom[c])]
            s0, e0 = sim.pos[c] + 1, sim.pos[c] + L
            name = f"{chrom}:{s0}-{e0}"
            cid[c] = name
            cls = f"LTR/{sim.superfamily[f]}/{sim.clade[f]}"
            d, se, ts, tv = k2p(st["ltr5"], st["ltr3"])
            ffa.write(f">{name}#{cls}\n{BASES[codes].tobytes().decode()}\n")
            ident = 1 - (ts + tv) / lt
            row = {"seq_id": f"{name}#{cls}", "seq_len": L, "status": "pass", "ltr5_start": 1, "ltr5_end": left,
                       "ltr3_start": L - right + 1, "ltr3_end": L, "ltr5_len": left, "ltr3_len": right, "flank5_len": 0, "flank3_len": 0,
                       "aln_len": lt, "n_sites": lt, "n_ts": ts, "n_tv": tv, "n_gapcols": 0, "identity": ident, "p_dist": (ts + tv) / lt,
                       "k2p": d, "k2p_se": se, "bitscore": "NA", "flank_margin_bits": "NA", "cigar": ".", "motif": "tg...ca",
                       "k2p_time": "NA" if math.isnan(d) else round(d / (2 * SIM_NUCLEUS_K2P_RATE)), "orientation": "+" if sim.strand[c] else "-",
                       "tsd": ".", "tsd_offset": "NA", "tsd_input": ".", "strand": "+" if sim.strand[c] else "-",
                       "family": label[c], "domains": ".", "nest_status": "."}
            ftsv.write("\t".join(_fmt(row[k]) for k in TSV_COLS) + "\n")
            sd, gchrom = "+" if sim.strand[c] else "-", escape_seqid(chrom)
            fgff.write(f"{gchrom}\t{GFF3_SOURCE}\tLTR_retrotransposon\t{s0}\t{e0}\t.\t{sd}\t.\tID=sim{n:07d};"
                       f"Name={escape_value(name)};"
                       f"classification={cls};order=LTR;superfamily={sim.superfamily[f]};clade={sim.clade[f]};"
                       f"family={label[c]};family_size={fam_size[label[c]]};depth=0;ltr5_len={lt};ltr3_len={lt};aln_len={lt};"
                       f"identity={ident:.6g};subs={ts + tv};ti={ts};tv={tv};K2P_d={_fmt(d)}\n")
            fgff.write(f"{gchrom}\t{GFF3_SOURCE}\tlong_terminal_repeat\t{s0}\t{s0 + lt - 1}\t.\t{sd}\t.\t"
                       f"ID=sim{n:07d}.lLTR;Parent=sim{n:07d}\n")
            fgff.write(f"{gchrom}\t{GFF3_SOURCE}\tlong_terminal_repeat\t{e0 - lt + 1}\t{e0}\t.\t{sd}\t.\t"
                       f"ID=sim{n:07d}.rLTR;Parent=sim{n:07d}\n")
            leaf = leafname[sim.lineage[c]]
            mode = "own" if sim.origin[c] == leaf else "foreign" if sim.origin[c] in leafname.values() else "ancestral"
            post = sim.tins[c] < sim.merge_time[sim.lineage[c]]
            src = sim.source_root(c)
            ftr.write(f"{name}\t{chrom}\t{leaf}\t{sim.chrom[c]}\tfam{f:04d}\t{sim.superfamily[f]}\t{sim.clade[f]}\t"
                      f"{sim.tins[c]:.6g}\t{sim.tb[c]:.6g}\t{sim.origin[c]}\t{mode}\t{'yes' if post else 'no'}\t"
                      f"{sim.exid[c]}\t{sim.root[c]}\t{lt}\t{ts}\t{tv}\t{_fmt(d)}\t"
                      f"{'.' if src < 0 else src}\n")
    with open(truthdir / "orthologs.tsv", "w") as f:
        f.write("copy_a\tcopy_b\tt_insert\tsplit_time\n")
        by_root = defaultdict(list)
        for c in surv:
            if sim.exid[c] < 0:
                by_root[sim.root[c]].append(c)
        for v in by_root.values():
            for i in range(len(v)):
                for j in range(i + 1, len(v)):
                    a, b = v[i], v[j]
                    if sim.lineage[a] != sim.lineage[b]:
                        f.write(f"{cid[a]}\t{cid[b]}\t{sim.tins[a]:.6g}\t{_mrca_split(sim, sim.lineage[a], sim.lineage[b]):.6g}\n")
    with open(truthdir / "paralogs.tsv", "w") as f:
        f.write("copy_dup\tcopy_src\texchange_id\texchange_time\n")
        for c in surv:
            p = sim.par[c]
            if sim.exid[c] >= 0 and sim.mode[c] == DUP and p in cid:
                f.write(f"{cid[c]}\t{cid[p]}\t{sim.exid[c]}\t{sim.tb[c]:.6g}\n")
    with open(truthdir / "exchanges.tsv", "w") as f:
        f.write("index\ttype\ttime\trecipient\tdonor\tchrom\tstart\tend\tn_moved\n")
        for e in sim.exchanges:
            R = sim.leaf_of[e["recipient"]]
            rows = [(R, e["recipient"], e["donor"])]
            if e["type"] == "reciprocal":
                rows.append((sim.leaf_of[e["donor"]], e["donor"], e["recipient"]))
            f.writelines(f"{e['index']}\t{e['type']}\t{e['time']:.6g}\t{rec}\t{don}\t{tok[(l, e['chrom'])]}\t"
                        f"{e['start']}\t{e['end']}\t{e['n_moved']}\n" for l, rec, don in rows)
    true_size = defaultdict(int)
    for c in surv:
        true_size[sim.fam[c]] += 1
    hist = {"scenario": sim.sc, "seed": sim.seed,
                "splits": [{"node": l.name, "time": l.t_to, "children": [sim.lin[c].name for c in l.children]}
                        for l in sim.lin if l.children],
                "mergers": sim.sc["mergers"], "merge_time": {sim.lin[l].name: t for l, t in sim.merge_time.items()},
                "families": [{"name": f"fam{f:04d}", "ltr_len": int(sim.ltr_len[f]), "int_len": int(sim.int_len[f]),
                               "mu": float(sim.mu_f[f]), "superfamily": sim.superfamily[f], "clade": sim.clade[f],
                               "n_copies": true_size.get(f, 0)} for f in range(len(sim.fam_w))],
                "n_coordinate_collisions_dropped": n_dup,
                "n_copies": {sim.lin[l].name: sum(sim.lineage[c] == l for c in surv) for l in sim.leaves}}
    (truthdir / "history.json").write_text(json.dumps(hist, indent=1) + "\n")
    return len(surv)


def simulate(scenario: dict, outdir, truthdir, seed: int) -> int:
    """Simulate `scenario` (its keys override the defaults) with `seed`, write the library to `outdir` and the
    truth to `truthdir`; the number of LTR-RTs written."""
    return write_outputs(Sim(merge_defaults(scenario), seed).run(), outdir, truthdir)
