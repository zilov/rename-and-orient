"""Chromosome name parsing and classification utilities."""
import re
from typing import Tuple

SEX_CHROMOSOME_SUFFIXES = {'W', 'X', 'Y', 'Z', 'B'}
_HAP_SUFFIX_RE = re.compile(r'_HAP\d+$', re.IGNORECASE)


def strip_hap_suffix(name: str) -> str:
    """Strip _HAPJ suffix from a scaffold name (e.g. SUPER_1_HAP1 -> SUPER_1)."""
    return _HAP_SUFFIX_RE.sub('', name)


def is_hap_contig(name: str) -> bool:
    """Return True if name carries a _HAPJ suffix produced by pretext-to-asm."""
    return bool(_HAP_SUFFIX_RE.search(name))


def is_unloc_contig(name: str) -> bool:
    """Check if sequence name is an unlocalized contig (SUPER_N_unloc_M format)."""
    return "_unloc_" in name


def is_sex_chromosome_suffix(suffix: str) -> bool:
    """Check if a suffix indicates a sex chromosome (W, Z, X, Y, Z1, Z2, etc.)."""
    suffix = strip_hap_suffix(suffix)
    if not suffix:
        return False
    first_char = suffix[0].upper()
    if first_char in SEX_CHROMOSOME_SUFFIXES:
        remaining = suffix[1:]
        return remaining == '' or remaining.isdigit()
    return False


def is_autosome_suffix(suffix: str) -> bool:
    """Check if suffix represents an autosome (purely numeric)."""
    return strip_hap_suffix(suffix).isdigit()


def extract_chromosome_suffix(name: str, prefix: str = "SUPER_") -> str:
    """Extract chromosome suffix from name by stripping the given prefix."""
    if name.startswith(prefix):
        return name[len(prefix):]
    return name


def parse_unloc_name(name: str) -> Tuple[str, int]:
    """Parse unlocalized contig name into (parent_chromosome, unloc_number).

    Handles haplotype-tagged unloc contigs (e.g. SUPER_1_unloc_2_HAP1), where
    the _HAPn suffix trails the unloc number rather than the parent name.
    """
    hap_match = _HAP_SUFFIX_RE.search(name)
    hap_suffix = hap_match.group(0) if hap_match else ""
    base_name = strip_hap_suffix(name)
    parts = base_name.split("_unloc_")
    parent = parts[0] + hap_suffix
    unloc_num = int(parts[1]) if len(parts) > 1 else 0
    return parent, unloc_num
