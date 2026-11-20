import builders
import pytest


@pytest.fixture(scope="session")
def data(tmp_path_factory):
    """Synthetic diploid with planted states (identical haplotypes apart from the events): builders.read_set."""
    return builders.read_set(tmp_path_factory.mktemp("synth"), 2)
