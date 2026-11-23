import logging
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parent.parent / "src" / "synltr"


@pytest.fixture
def src_root() -> Path:
    """The package source folder."""
    return SRC


@pytest.fixture(autouse=True)
def _pristine_log():
    """After each test, the synltr logger as before any command ran (cli.main's log.setup gives it a handler
    and stops propagation, which would hide later tests' records from caplog)."""
    yield
    log = logging.getLogger("synltr")
    log.handlers.clear()
    log.setLevel(logging.NOTSET)
    log.propagate = True


@pytest.fixture(scope="session")
def toy_fastas(tmp_path_factory):
    """One directory of N-only FASTAs shared by the session: catalog_toy.write_run writes
    <dir>/<layout>/<label>.fa once per seqid layout and reuses it after."""
    return tmp_path_factory.mktemp("toy_fastas")


@pytest.fixture(scope="session")
def admixed(tmp_path_factory):
    """LTRmix's admixed synthetic trio (builders.trio, seed 3, as its own tests ran it) screened with 2 threads:
    (scenario, run folder)."""
    import builders

    from synltr.io.samples import read_samples
    from synltr.routes.library.run import screen_pairs

    d = tmp_path_factory.mktemp("admixed")
    sc = builders.trio(d / "in", seed=3, admix=True)
    screen_pairs(read_samples(sc.sheet), d / "run", 2)
    return sc, d / "run"


def _screened(d: Path, sc) -> tuple:
    """Scenario `sc`'s sheet screened into d/run with 2 threads: (scenario, run folder)."""
    from synltr.io.samples import read_samples
    from synltr.routes.library.run import screen_pairs

    screen_pairs(read_samples(sc.sheet), d / "run", 2)
    return sc, d / "run"


@pytest.fixture(scope="session")
def null_screen(tmp_path_factory):
    """LTRmix's trio without gene flow (builders.trio, seed 3, admix False) screened with 2 threads: (scenario, run
    folder)."""
    import builders

    d = tmp_path_factory.mktemp("null")
    return _screened(d, builders.trio(d / "in", seed=3, admix=False))


@pytest.fixture(scope="session")
def four_screen(tmp_path_factory):
    """builders.trio_second_donor (seed 5: R, S, D and D's sister E) screened with 2 threads: (scenario, run
    folder)."""
    import builders

    d = tmp_path_factory.mktemp("four")
    return _screened(d, builders.trio_second_donor(d / "in", seed=5))


@pytest.fixture(scope="session")
def unplaced_screen(tmp_path_factory):
    """The admixed trio (seed 3) with 4 records whose names hold no coordinates (mutated copies of placed records, no
    table row: tests/golden/ltrmix_stats.npz e2e_unplaced_*) appended to R's library and 4 to D's, screened with 2
    threads: (the scenario with that sheet, run folder)."""
    import builders
    import numpy as np

    d = tmp_path_factory.mktemp("unplaced")
    sc = builders.trio(d / "in", seed=3, admix=True)
    with np.load(Path(__file__).resolve().parent / "golden" / "ltrmix_stats.npz") as z:
        stored = {k: z[k] for k in z.files if k.startswith("e2e_unplaced_")}
    lines = ["label\tlibrary"]
    for g in "RSD":
        fa, tb = sc.libraries[g].split(",")
        seq, offs = stored[f"e2e_unplaced_{g}_seq"], stored[f"e2e_unplaced_{g}_offs"]
        extra = "".join(f">contig_{k}\n{seq[offs[k]:offs[k + 1]].tobytes().decode()}\n" for k in range(offs.size - 1))
        (d / "in" / f"{g}u.fa").write_text(Path(fa).read_text() + extra)
        lines.append(f"{g}\t{d / 'in' / f'{g}u.fa'},{tb}")
    sheet = d / "in" / "unplaced.tsv"
    sheet.write_text("\n".join(lines) + "\n")
    return _screened(d, builders.TrioScenario(sheet, sc.libraries, sc.truth, sc.tracts, sc.inversion))
