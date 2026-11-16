"""FASTA reading/writing and sequence helpers."""
import gzip
import logging
import os
import shutil
import tempfile
import weakref
from collections.abc import Iterator, Mapping
from pathlib import Path
from typing import BinaryIO, NamedTuple, TextIO

from .models import FinalChromosomeAssignment, UnlocMapping
from .names import _unloc_output_name, group_unlocs_by_parent

logger = logging.getLogger(__name__)


_WRITE_BLOCK_LINES = 8192  # wrapped lines joined per write() call


class FastaRecord(NamedTuple):
    """Location of one FASTA record inside the (plain) file.

    Attributes:
        offset: Byte offset of the first byte after the header line.
        nbytes: Number of bytes spanning all sequence lines of the record.
        length: Sequence length (bytes) without newlines and edge whitespace.
        edge_ws: True if some sequence line has leading/trailing whitespace
            other than the line terminator (forces the slow per-line read path).
    """

    offset: int
    nbytes: int
    length: int
    edge_ws: bool


def _unlink_quietly(path: str) -> None:
    """Remove a temporary file, ignoring a missing file."""
    try:
        os.unlink(path)
    except FileNotFoundError:
        pass


class FastaIndex(Mapping[str, str]):
    """Random-access index of a FASTA file; sequences are read lazily, one at a time.

    Behaves as a read-only mapping name -> sequence (in file order), but nothing is
    kept in memory except the index: ``index[name]`` / ``index.get(name)`` read just
    that record from disk. Iterating over the index yields names. Use it as a context
    manager (or call :meth:`close`) so that the temporary decompressed copy of a
    gzipped input is removed.
    """

    def __init__(self, path: Path, records: dict[str, FastaRecord], tmp_path: str | None) -> None:
        """Create an index; use :func:`index_fasta` instead of calling this directly."""
        self.path = path
        self._records = records
        self._tmp_path = tmp_path
        self._read_path = tmp_path if tmp_path is not None else str(path)
        self._fh: BinaryIO | None = None
        self._finalizer = (
            weakref.finalize(self, _unlink_quietly, tmp_path) if tmp_path is not None else None
        )

    @property
    def names(self) -> list[str]:
        """Sequence names in file order."""
        return list(self._records)

    @property
    def total_length(self) -> int:
        """Total length of all sequences."""
        return sum(r.length for r in self._records.values())

    def length(self, name: str) -> int:
        """Return the sequence length of one record without reading it."""
        return self._records[name].length

    def __len__(self) -> int:
        """Return the number of records."""
        return len(self._records)

    def __iter__(self) -> Iterator[str]:
        """Iterate over names in file order."""
        return iter(self._records)

    def __contains__(self, name: object) -> bool:
        """Return True if a record with this name exists (does not read the file)."""
        return name in self._records

    def __getitem__(self, name: str) -> str:
        """Read the sequence of one record from disk."""
        return self.read_sequence(name)

    def read_sequence(self, name: str) -> str:
        """Read one record's sequence (lines joined, line edges stripped).

        Args:
            name: Sequence name.

        Returns:
            The sequence as a string; empty for a record without sequence lines.

        Raises:
            KeyError: If the name is not in the index.
            ValueError: If the index is closed and the data source is gone.
        """
        rec = self._records[name]
        if rec.nbytes == 0:
            return ""
        if self._fh is None:
            self._fh = open(self._read_path, "rb")
        self._fh.seek(rec.offset)
        data = self._fh.read(rec.nbytes)
        if rec.edge_ws:
            # Same semantics as the old text reader: strip every line, skip blank ones.
            return b"".join(line.strip() for line in data.splitlines()).decode()
        data = data.translate(None, b"\r\n")  # rebinding frees the raw block early
        return data.decode()

    def close(self) -> None:
        """Close the file handle and delete the temporary decompressed copy, if any."""
        if self._fh is not None:
            self._fh.close()
            self._fh = None
        if self._finalizer is not None:
            self._finalizer()

    def __enter__(self) -> "FastaIndex":
        """Enter the context manager."""
        return self

    def __exit__(self, *exc_info: object) -> None:
        """Close the index on leaving the context."""
        self.close()


def _decompress_to_temp(path: Path) -> str:
    """Decompress a gzip file to a temporary file in the default temp dir ($TMPDIR)."""
    fd, tmp_path = tempfile.mkstemp(prefix="rename_and_orient_", suffix=".fa")
    try:
        with os.fdopen(fd, "wb") as out, gzip.open(path, "rb") as src:
            shutil.copyfileobj(src, out, length=1 << 20)
    except BaseException:
        _unlink_quietly(tmp_path)
        raise
    return tmp_path


def index_fasta(fasta_path: Path) -> FastaIndex:
    """Index a FASTA file in a single pass without keeping sequences in memory.

    A gzipped file (name ends with ``.gz``) is decompressed once to a temporary file
    (in ``$TMPDIR``), which is removed by :meth:`FastaIndex.close`. Plain files are
    read in place.

    Args:
        fasta_path: Path to a FASTA file (plain or gzip).

    Returns:
        A :class:`FastaIndex` (use it as a context manager).

    Raises:
        ValueError: On an empty header name or duplicate sequence names.
    """
    fasta_path = Path(fasta_path)
    tmp_path = _decompress_to_temp(fasta_path) if str(fasta_path).endswith(".gz") else None
    read_path = tmp_path if tmp_path is not None else str(fasta_path)
    records: dict[str, FastaRecord] = {}
    try:
        with open(read_path, "rb") as f:
            name: str | None = None
            start = pos = 0
            length = 0
            edge_ws = False
            for line in f:
                end = pos + len(line)
                stripped = line.strip()
                if stripped.startswith(b">"):
                    if name is not None:
                        records[name] = FastaRecord(start, pos - start, length, edge_ws)
                    fields = stripped[1:].split()
                    if not fields:
                        raise ValueError(
                            f"{fasta_path}: empty sequence name in header at byte {pos}"
                        )
                    name = fields[0].decode()
                    if name in records:
                        raise ValueError(f"{fasta_path}: duplicate sequence name '{name}'")
                    start = end
                    length = 0
                    edge_ws = False
                elif stripped:
                    length += len(stripped)
                    if len(stripped) != len(line.rstrip(b"\r\n")):
                        edge_ws = True
                pos = end
            if name is not None:
                records[name] = FastaRecord(start, pos - start, length, edge_ws)
    except BaseException:
        if tmp_path is not None:
            _unlink_quietly(tmp_path)
        raise
    return FastaIndex(fasta_path, records, tmp_path)


def read_fasta(fasta_path: Path) -> dict[str, str]:
    """Read a whole FASTA file into memory (supports gzip compression).

    Prefer :func:`index_fasta` for large genomes: this loads every sequence.

    Args:
        fasta_path: Path to FASTA file

    Returns:
        Dictionary mapping sequence names to sequences

    Raises:
        ValueError: If sequence names are duplicated.
    """
    with index_fasta(fasta_path) as index:
        return {name: index.read_sequence(name) for name in index}


# Full IUPAC nucleotide complement table (ambiguity codes included), both cases.
# Same mapping as IUPAC_COMPLEMENT in agp-tpf-utils (tola/fasta/simple.py, MIT).
# Characters outside the table (e.g. U, '-', '*') pass through unchanged.
_RC_TABLE = str.maketrans(
    "ACGTRYMKSWHBVDNacgtrymkswhbvdn",
    "TGCAYRKMSWDVBHNtgcayrkmswdvbhn",
)


def reverse_complement(seq: str) -> str:
    """Return reverse complement of a DNA sequence.

    IUPAC ambiguity codes are complemented (e.g. ``R`` <-> ``Y``, ``K`` <-> ``M``,
    ``H`` <-> ``D``, ``B`` <-> ``V``); case is preserved; any other character is
    passed through unchanged.

    Args:
        seq: DNA sequence string

    Returns:
        Reverse complement sequence
    """
    return seq.translate(_RC_TABLE)[::-1]


def calculate_genome_length(fasta_path: Path) -> int:
    """Calculate total genome length from FASTA file (streaming, nothing is accumulated).

    Args:
        fasta_path: Path to FASTA file (plain or gzip)

    Returns:
        Total length of all sequences
    """
    is_gzipped = str(fasta_path).endswith(".gz")
    total = 0
    in_record = False
    with (gzip.open(fasta_path, "rb") if is_gzipped else open(fasta_path, "rb")) as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith(b">"):
                in_record = True
            elif in_record:
                total += len(stripped)
    return total


def _write_wrapped(f: TextIO, seq: str, line_width: int) -> None:
    """Write a sequence wrapped at line_width, joining many lines per write() call."""
    block = line_width * _WRITE_BLOCK_LINES
    for start in range(0, len(seq), block):
        chunk = seq[start:start + block]
        lines = [chunk[i:i + line_width] for i in range(0, len(chunk), line_width)]
        f.write("\n".join(lines) + "\n")


def write_fasta(
    sequences: Mapping[str, str],
    assignments: list[FinalChromosomeAssignment],
    unloc_mappings: list[UnlocMapping],
    output_path: Path,
    output_prefix: str = "SUPER_",
    line_width: int = 60
) -> None:
    """Write FASTA file with renamed chromosomes and reverse complement where needed.

    Orientation is taken from a.needs_reverse_complement (chromosomes) and
    unloc.needs_reverse_complement (unlocalized contigs).

    Records are fetched one at a time with ``sequences.get(name, '')``, so a
    :class:`FastaIndex` keeps peak memory near the size of the largest record.

    Args:
        sequences: Original sequences: a dict, or a FastaIndex (lazy mapping)
        assignments: List of FinalChromosomeAssignment (sorted for output)
        unloc_mappings: List of UnlocMapping for unloc contigs
        output_path: Path to output FASTA file
        output_prefix: Prefix for output chromosome names
        line_width: Line width for sequence output (default 60)
    """
    parent_suffix_lookup = {a.original_name: a.new_suffix for a in assignments}
    processed_sequences = set()

    unloc_by_parent = group_unlocs_by_parent(unloc_mappings)

    with open(output_path, 'w') as f:
        for a in assignments:
            seq = sequences.get(a.original_name, '')
            if not seq:
                logger.warning("No sequence found for %s", a.original_name)
                continue

            processed_sequences.add(a.original_name)

            if a.needs_reverse_complement:
                seq = reverse_complement(seq)

            f.write(f">{a.new_name}\n")
            _write_wrapped(f, seq, line_width)
            del seq  # release the record before the next one is read

            for unloc in unloc_by_parent.get(a.original_name, []):
                unloc_seq = sequences.get(unloc.contig_name, '')
                if not unloc_seq:
                    logger.warning("No sequence found for %s", unloc.contig_name)
                    continue

                processed_sequences.add(unloc.contig_name)

                if unloc.needs_reverse_complement:
                    unloc_seq = reverse_complement(unloc_seq)

                new_name, _ = _unloc_output_name(unloc, parent_suffix_lookup, output_prefix)

                f.write(f">{new_name}\n")
                _write_wrapped(f, unloc_seq, line_width)
                del unloc_seq

        for seq_name in sequences:
            if seq_name in processed_sequences:
                continue

            f.write(f">{seq_name}\n")
            _write_wrapped(f, sequences[seq_name], line_width)

    logger.info("FASTA written to: %s", output_path)
