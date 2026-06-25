import pytest
from pathlib import Path


@pytest.fixture
def test_data_dir():
    """Path to the test_data directory containing input files."""
    return Path(__file__).parent / "test_data"


@pytest.fixture
def paf_file(test_data_dir):
    """PAF alignment file from test_data."""
    return test_data_dir / "test.paf"


@pytest.fixture
def fasta_file(test_data_dir):
    """FASTA file from test_data."""
    return test_data_dir / "test.fa"


@pytest.fixture
def chromosome_list_file(test_data_dir):
    """Chromosome list CSV file from test_data."""
    return test_data_dir / "test.chromosome.list.csv"