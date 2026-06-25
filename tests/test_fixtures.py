"""Test fixtures validation."""


def test_fixtures_exist(test_data_dir, paf_file, fasta_file):
    """Test that all fixtures point to existing files."""
    assert test_data_dir.exists()
    assert test_data_dir.is_dir()
    assert paf_file.exists()
    assert fasta_file.exists()