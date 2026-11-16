from pathlib import Path

import tomllib

from rename_and_orient import cli


def test_version_matches_pyproject():
    pyproject = Path(__file__).parent.parent / "pyproject.toml"
    assert cli.__version__ == tomllib.loads(pyproject.read_text())["project"]["version"]
