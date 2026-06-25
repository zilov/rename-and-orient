"""Unit tests for individual functions in rename_and_orient.py."""
import pytest
from pathlib import Path
from rename_and_orient import (
    read_fasta,
    parse_paf,
    filter_paf_records,
    validate_paf_fasta_consistency,
    build_chromosome_mappings,
    reverse_complement,
    extract_chromosome_suffix,
    strip_hap_suffix,
    is_hap_contig,
    write_fasta,
    write_chromosome_list,
    build_unloc_mappings,
    sort_assignments_for_output,
    is_sex_chromosome_suffix,
    is_autosome_suffix,
    resolve_chromosome_assignments,
    detect_reference_prefix,
    merge_intervals,
    calculate_target_alignments,
    load_mapping_table_assignments,
    PAFRecord,
    ChromosomeMapping,
    FinalChromosomeAssignment,
)


@pytest.fixture
def sample_fasta_content():
    """Sample FASTA content for testing."""
    return """>SUPER_1
ATCGATCG
>SUPER_2
GCTAGCTA
>SUPER_W
NNNNAAAA
"""


@pytest.fixture
def sample_fasta_dict():
    """Expected FASTA dictionary from sample_fasta_content."""
    return {
        "SUPER_1": "ATCGATCG",
        "SUPER_2": "GCTAGCTA",
        "SUPER_W": "NNNNAAAA"
    }


@pytest.fixture
def sample_paf_content():
    """Sample PAF content for testing."""
    return """SUPER_1\t8\t0\t8\t+\tchr_1\t100\t0\t8\t8\t8\t60
SUPER_2\t8\t0\t8\t-\tchr_2\t100\t0\t8\t8\t8\t60
SUPER_W\t8\t0\t8\t+\tchr_W\t100\t0\t8\t8\t8\t60
"""


@pytest.fixture
def sample_paf_records():
    """Expected PAF records from sample_paf_content."""
    return [
        PAFRecord("SUPER_1", 8, 0, 8, "+", "chr_1", 100, 0, 8, 8, 8, 60),
        PAFRecord("SUPER_2", 8, 0, 8, "-", "chr_2", 100, 0, 8, 8, 8, 60),
        PAFRecord("SUPER_W", 8, 0, 8, "+", "chr_W", 100, 0, 8, 8, 8, 60),
    ]


class TestReadFasta:
    """Test read_fasta function."""

    def test_read_fasta_basic(self, tmp_path, sample_fasta_content, sample_fasta_dict):
        """Test basic FASTA reading."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(sample_fasta_content)

        result = read_fasta(fasta_file)
        assert result == sample_fasta_dict

    def test_read_fasta_gzipped(self, tmp_path, sample_fasta_content, sample_fasta_dict):
        """Test reading gzipped FASTA."""
        import gzip
        fasta_file = tmp_path / "test.fa.gz"
        with gzip.open(fasta_file, 'wt') as f:
            f.write(sample_fasta_content)

        result = read_fasta(fasta_file)
        assert result == sample_fasta_dict


class TestParsePaf:
    """Test parse_paf function."""

    def test_parse_paf_basic(self, tmp_path, sample_paf_content, sample_paf_records):
        """Test basic PAF parsing."""
        paf_file = tmp_path / "test.paf"
        paf_file.write_text(sample_paf_content)

        result = parse_paf(paf_file)
        assert len(result) == 3
        assert result[0].query_name == "SUPER_1"
        assert result[0].query_length == 8
        assert result[0].strand == "+"
        assert result[0].target_name == "chr_1"


class TestFilterPafRecords:
    """Test filter_paf_records function."""

    def test_filter_with_auto_detect(self, sample_paf_records):
        """Test filtering with auto-detected target prefix."""
        filtered, prefix = filter_paf_records(sample_paf_records, "SUPER_")
        assert len(filtered) == 3
        assert prefix == "chr_"

    def test_filter_with_custom_prefix(self, sample_paf_records):
        """Test filtering with custom target prefix."""
        # Create records with different target prefix
        records = [
            PAFRecord("SUPER_1", 8, 0, 8, "+", "scaffold_1", 100, 0, 8, 8, 8, 60),
        ]
        filtered, prefix = filter_paf_records(records, "SUPER_", "scaffold_")
        assert len(filtered) == 1
        assert prefix == "scaffold_"


class TestValidatePafFastaConsistency:
    """Test validate_paf_fasta_consistency function."""

    def test_consistent_data(self, sample_paf_records, sample_fasta_dict):
        """Test with consistent PAF and FASTA data."""
        is_valid, in_paf_not_fasta, in_fasta_not_paf = validate_paf_fasta_consistency(
            sample_paf_records, sample_fasta_dict, "SUPER_"
        )
        assert is_valid
        assert len(in_paf_not_fasta) == 0
        assert len(in_fasta_not_paf) == 0

    def test_inconsistent_data(self, sample_paf_records, sample_fasta_dict):
        """Test with inconsistent PAF and FASTA data."""
        # Add a chromosome in PAF but not in FASTA
        extra_record = PAFRecord("SUPER_3", 8, 0, 8, "+", "chr_3", 100, 0, 8, 8, 8, 60)
        paf_records = sample_paf_records + [extra_record]

        is_valid, in_paf_not_fasta, in_fasta_not_paf = validate_paf_fasta_consistency(
            paf_records, sample_fasta_dict, "SUPER_"
        )
        assert not is_valid
        assert "SUPER_3" in in_paf_not_fasta


class TestReverseComplement:
    """Test reverse_complement function."""

    def test_reverse_complement_basic(self):
        """Test basic reverse complement."""
        seq = "ATCG"
        result = reverse_complement(seq)
        assert result == "CGAT"

    def test_reverse_complement_lower_case(self):
        """Test reverse complement with lowercase."""
        seq = "atcg"
        result = reverse_complement(seq)
        assert result == "cgat"

    def test_reverse_complement_with_n(self):
        """Test reverse complement with N bases."""
        seq = "ATCGN"
        result = reverse_complement(seq)
        assert result == "NCGAT"

    def test_reverse_complement_empty(self):
        """Test reverse complement of empty string."""
        result = reverse_complement("")
        assert result == ""


class TestChromosomeSuffixFunctions:
    """Test chromosome suffix related functions."""

    @pytest.mark.parametrize("name,prefix,expected", [
        ("SUPER_1", "SUPER_", "1"),
        ("chr_W", "chr_", "W"),
        ("scaffold_10", "scaffold_", "10"),
        ("chrZ1", "chr", "Z1"),
        # HAP suffix is preserved by extract_chromosome_suffix (stripping is for classification only)
        ("SUPER_1_HAP1", "SUPER_", "1_HAP1"),
        ("SUPER_W_HAP2", "SUPER_", "W_HAP2"),
    ])
    def test_extract_chromosome_suffix(self, name, prefix, expected):
        """Test extract_chromosome_suffix with various inputs."""
        result = extract_chromosome_suffix(name, prefix)
        assert result == expected

    @pytest.mark.parametrize("name,expected", [
        ("SUPER_1_HAP1", "SUPER_1"),
        ("SUPER_W_HAP2", "SUPER_W"),
        ("SUPER_1", "SUPER_1"),
        ("1_HAP1", "1"),
        ("W_HAP2", "W"),
    ])
    def test_strip_hap_suffix(self, name, expected):
        assert strip_hap_suffix(name) == expected

    @pytest.mark.parametrize("name,expected", [
        ("SUPER_1_HAP1", True),
        ("SUPER_W_HAP2", True),
        ("SUPER_1", False),
        ("SUPER_1_unloc_1", False),
    ])
    def test_is_hap_contig(self, name, expected):
        assert is_hap_contig(name) == expected

    @pytest.mark.parametrize("suffix,expected", [
        ("W", True),
        ("Z", True),
        ("X", True),
        ("Y", True),
        ("Z1", True),
        ("Z2", True),
        ("1", False),
        ("10", False),
        ("ABC", False),
        # HAP-carrying suffixes (as returned by extract_chromosome_suffix)
        ("W_HAP1", True),
        ("Z1_HAP2", True),
        ("1_HAP1", False),
    ])
    def test_is_sex_chromosome_suffix(self, suffix, expected):
        """Test is_sex_chromosome_suffix."""
        result = is_sex_chromosome_suffix(suffix)
        assert result == expected

    @pytest.mark.parametrize("suffix,expected", [
        ("1", True),
        ("10", True),
        ("123", True),
        ("W", False),
        ("Z1", False),
        ("ABC", False),
        # HAP-carrying suffixes (as returned by extract_chromosome_suffix)
        ("1_HAP1", True),
        ("10_HAP2", True),
        ("W_HAP1", False),
    ])
    def test_is_autosome_suffix(self, suffix, expected):
        """Test is_autosome_suffix."""
        result = is_autosome_suffix(suffix)
        assert result == expected


class TestBuildChromosomeMappings:
    """Test build_chromosome_mappings function."""

    def test_build_mappings_basic(self, sample_paf_records):
        """Test basic chromosome mapping building."""
        mappings = build_chromosome_mappings(sample_paf_records, 0.5, "chr_")

        assert len(mappings) == 3
        # Check first mapping
        mapping = mappings[0]
        assert isinstance(mapping, ChromosomeMapping)
        assert mapping.query_name == "SUPER_1"
        assert mapping.target_name == "chr_1"
        assert mapping.coverage == 1.0  # 8/8
        assert not mapping.needs_reverse_complement  # strand +

    def test_build_mappings_with_rc(self, sample_paf_records):
        """Test mapping with reverse complement needed."""
        # SUPER_2 has strand -, so should need RC
        mappings = build_chromosome_mappings(sample_paf_records, 0.5, "chr_")

        super2_mapping = next(m for m in mappings if m.query_name == "SUPER_2")
        assert super2_mapping.needs_reverse_complement

    def test_build_mappings_coverage_filter(self):
        """Test coverage filtering."""
        # Create record with low coverage
        low_coverage_record = PAFRecord("SUPER_1", 100, 0, 10, "+", "chr_1", 100, 0, 10, 10, 10, 60)
        mappings = build_chromosome_mappings([low_coverage_record], 0.5, "chr_")

        # Coverage is 10/100 = 0.1, which is below 0.5 threshold
        assert len(mappings) == 0


class TestResolveChromosomeAssignments:
    """Test resolve_chromosome_assignments function."""

    def test_resolve_basic_autosomes(self, sample_fasta_dict):
        """Test basic autosome assignment resolution."""
        mappings = [
            ChromosomeMapping("SUPER_1", 8, "chr_1", 8, 1.0, 8, 0, False, "chr_"),
            ChromosomeMapping("SUPER_2", 8, "chr_2", 8, 1.0, 0, 8, True, "chr_"),
        ]

        assignments, rc_lookup = resolve_chromosome_assignments(
            mappings, sample_fasta_dict, "SUPER_", "SUPER_"
        )

        assert len(assignments) == 3  # SUPER_1, SUPER_2, SUPER_W (unmapped)
        assert assignments[0].original_name == "SUPER_1"
        assert assignments[0].new_name == "SUPER_1"
        assert not assignments[0].needs_reverse_complement

        assert assignments[1].original_name == "SUPER_2"
        assert assignments[1].new_name == "SUPER_2"
        assert assignments[1].needs_reverse_complement

        # SUPER_W is unmapped, keeps original
        super_w = next(a for a in assignments if a.original_name == "SUPER_W")
        assert super_w.new_name == "SUPER_W"
        assert super_w.is_sex_chromosome

    def test_resolve_sex_chromosome_mapping(self, sample_fasta_dict):
        """Test sex chromosome mapping."""
        mappings = [
            ChromosomeMapping("SUPER_W", 8, "chr_W", 8, 1.0, 8, 0, False, "chr_"),
        ]

        assignments, rc_lookup = resolve_chromosome_assignments(
            mappings, sample_fasta_dict, "SUPER_", "SUPER_"
        )

        super_w = next(a for a in assignments if a.original_name == "SUPER_W")
        assert super_w.new_name == "SUPER_W"  # Keeps original suffix
        assert super_w.is_sex_chromosome

    def test_resolve_hap_autosomes(self):
        """HAP-suffixed autosomes are renamed correctly; _HAPJ is preserved in new_name."""
        fasta = {
            "SUPER_1_HAP2": "ATCGATCG",
            "SUPER_2_HAP2": "GCTAGCTA",
        }
        mappings = [
            ChromosomeMapping("SUPER_1_HAP2", 8, "chr_1", 8, 1.0, 8, 0, False, "chr_"),
            ChromosomeMapping("SUPER_2_HAP2", 8, "chr_2", 8, 1.0, 0, 8, True, "chr_"),
        ]

        assignments, rc_lookup = resolve_chromosome_assignments(
            mappings, fasta, "SUPER_", "SUPER_"
        )

        a1 = next(a for a in assignments if a.original_name == "SUPER_1_HAP2")
        assert a1.new_name == "SUPER_1_HAP2"
        assert a1.new_suffix == "1"  # bare for sorting
        assert not a1.is_sex_chromosome
        assert not a1.needs_reverse_complement

        a2 = next(a for a in assignments if a.original_name == "SUPER_2_HAP2")
        assert a2.new_name == "SUPER_2_HAP2"
        assert a2.new_suffix == "2"
        assert a2.needs_reverse_complement

    def test_resolve_hap_sex_chromosome(self):
        """HAP-suffixed sex chromosomes are classified correctly; new_suffix is bare."""
        fasta = {
            "SUPER_W_HAP2": "NNNNAAAA",
            "SUPER_1_HAP2": "ATCGATCG",
        }
        mappings = [
            ChromosomeMapping("SUPER_W_HAP2", 8, "chr_W", 8, 1.0, 8, 0, False, "chr_"),
            ChromosomeMapping("SUPER_1_HAP2", 8, "chr_1", 8, 1.0, 8, 0, False, "chr_"),
        ]

        assignments, rc_lookup = resolve_chromosome_assignments(
            mappings, fasta, "SUPER_", "SUPER_"
        )

        w = next(a for a in assignments if a.original_name == "SUPER_W_HAP2")
        assert w.new_name == "SUPER_W_HAP2"
        assert w.new_suffix == "W"   # bare for sorting/unloc
        assert w.is_sex_chromosome


class TestDetectReferencePrefix:
    """Test detect_reference_prefix function."""

    @pytest.mark.parametrize("target_names,expected_prefix", [
        # Standard chr_ (Ensembl / insect assemblies)
        (["chr_1", "chr_2", "chr_3", "chr_W"], "chr_"),
        # Standard chr (UCSC)
        (["chr1", "chr2", "chr3", "chrW"], "chr"),
        # SUPER_ reference
        (["SUPER_1", "SUPER_2", "SUPER_3", "SUPER_W"], "SUPER_"),
        # scaffold_ reference
        (["scaffold_1", "scaffold_2", "scaffold_3"], "scaffold_"),
        # Mixed: most common wins
        (["chr_1", "chr_2", "chr_3", "unplaced_1"], "chr_"),
    ])
    def test_detect_prefix_from_target_names(self, target_names, expected_prefix):
        """Test that correct prefix is detected from target names."""
        records = [
            PAFRecord(f"query_{i}", 100, 0, 100, "+", name, 100, 0, 100, 100, 100, 60)
            for i, name in enumerate(target_names)
        ]
        result = detect_reference_prefix(records)
        assert result == expected_prefix

    def test_detect_prefix_empty_records(self):
        """Test that ValueError is raised when no records are provided."""
        with pytest.raises(ValueError, match="reference chromosome prefix"):
            detect_reference_prefix([])

    def test_detect_prefix_no_alpha_prefix(self):
        """Test that ValueError is raised when targets have no alphabetic prefix."""
        records = [
            PAFRecord("q1", 100, 0, 100, "+", "1", 100, 0, 100, 100, 100, 60),
        ]
        with pytest.raises(ValueError, match="reference chromosome prefix"):
            detect_reference_prefix(records)

    def test_detect_prefix_many_unlocs_dont_override_chromosomes(self):
        """Test that 1000 unloc contigs don't override 10 main chromosomes.
        
        Unloc names like SUPER_1_unloc_1 still yield prefix 'SUPER_' (regex
        stops at first digit), so counts are additive and SUPER_ still wins.
        """
        chr_records = [
            PAFRecord(f"q_{i}", 100, 0, 100, "+", f"SUPER_{i}", 100, 0, 100, 100, 100, 60)
            for i in range(1, 11)  # 10 chromosomes: SUPER_1 .. SUPER_10
        ]
        unloc_records = [
            PAFRecord(f"qu_{i}", 100, 0, 100, "+", f"SUPER_{i // 100 + 1}_unloc_{i % 100 + 1}", 100, 0, 100, 100, 100, 60)
            for i in range(1000)  # 1000 unloc contigs
        ]
        result = detect_reference_prefix(chr_records + unloc_records)
        assert result == "SUPER_"


class TestFilterPafRecordsAutoDetect:
    """Test filter_paf_records auto-detection of reference prefix."""

    def test_autodetect_super_reference(self):
        """Test auto-detection when reference uses SUPER_ prefix."""
        records = [
            PAFRecord("scaffold_1", 100, 0, 100, "+", "SUPER_1", 100, 0, 100, 100, 100, 60),
            PAFRecord("scaffold_2", 100, 0, 100, "-", "SUPER_2", 100, 0, 100, 100, 100, 60),
            PAFRecord("scaffold_3", 100, 0, 100, "+", "SUPER_W", 100, 0, 100, 100, 100, 60),
        ]
        filtered, prefix = filter_paf_records(records, "scaffold_")
        assert prefix == "SUPER_"
        assert len(filtered) == 3

    def test_autodetect_chr_reference(self):
        """Test auto-detection when reference uses chr prefix (no underscore)."""
        records = [
            PAFRecord("SUPER_1", 100, 0, 100, "+", "chr1", 100, 0, 100, 100, 100, 60),
            PAFRecord("SUPER_2", 100, 0, 100, "-", "chr2", 100, 0, 100, 100, 100, 60),
        ]
        filtered, prefix = filter_paf_records(records, "SUPER_")
        assert prefix == "chr"
        assert len(filtered) == 2

    def test_autodetect_chr_underscore_reference(self):
        """Test auto-detection when reference uses chr_ prefix."""
        records = [
            PAFRecord("SUPER_1", 100, 0, 100, "+", "chr_1", 100, 0, 100, 100, 100, 60),
            PAFRecord("SUPER_2", 100, 0, 100, "-", "chr_2", 100, 0, 100, 100, 100, 60),
        ]
        filtered, prefix = filter_paf_records(records, "SUPER_")
        assert prefix == "chr_"
        assert len(filtered) == 2

    def test_provided_prefix_overrides_autodetect(self):
        """Test that explicitly provided prefix overrides auto-detection."""
        records = [
            PAFRecord("SUPER_1", 100, 0, 100, "+", "chr_1", 100, 0, 100, 100, 100, 60),
        ]
        # Provide wrong prefix on purpose — should use it, result in 0 filtered
        filtered, prefix = filter_paf_records(records, "SUPER_", "scaffold_")
        assert prefix == "scaffold_"
        assert len(filtered) == 0

    def test_filter_excludes_non_matching_query_prefix(self):
        """Test that records with wrong query prefix are excluded."""
        records = [
            PAFRecord("SUPER_1", 100, 0, 100, "+", "chr_1", 100, 0, 100, 100, 100, 60),
            PAFRecord("contig_1", 100, 0, 100, "+", "chr_2", 100, 0, 100, 100, 100, 60),
        ]
        filtered, prefix = filter_paf_records(records, "SUPER_")
        assert len(filtered) == 1
        assert filtered[0].query_name == "SUPER_1"

class TestMergeIntervals:
    """Test merge_intervals function."""

    def test_no_overlap(self):
        """Non-overlapping intervals are returned unchanged."""
        intervals = [(0, 10), (20, 30), (40, 50)]
        result = merge_intervals(intervals)
        assert result == [(0, 10), (20, 30), (40, 50)]

    def test_fully_nested(self):
        """Interval fully contained within another is absorbed."""
        intervals = [(0, 100), (10, 50), (20, 80)]
        result = merge_intervals(intervals)
        assert result == [(0, 100)]
        assert sum(e - s for s, e in result) == 100

    def test_partial_overlap(self):
        """Partially overlapping intervals are merged into one."""
        intervals = [(0, 60), (40, 100)]
        result = merge_intervals(intervals)
        assert result == [(0, 100)]

    def test_adjacent_intervals(self):
        """Adjacent (touching) intervals are merged."""
        intervals = [(0, 50), (50, 100)]
        result = merge_intervals(intervals)
        assert result == [(0, 100)]

    def test_identical_intervals(self):
        """Many identical intervals collapse to one."""
        intervals = [(40302, 110580)] * 20
        result = merge_intervals(intervals)
        assert result == [(40302, 110580)]
        assert sum(e - s for s, e in result) == 110580 - 40302

    def test_unsorted_input(self):
        """Unsorted input is handled correctly."""
        intervals = [(50, 100), (0, 60)]
        result = merge_intervals(intervals)
        assert result == [(0, 100)]

    def test_empty_input(self):
        """Empty list returns empty list."""
        assert merge_intervals([]) == []

    def test_single_interval(self):
        """Single interval is returned as-is."""
        assert merge_intervals([(5, 15)]) == [(5, 15)]


class TestCalculateTargetAlignmentsNested:
    """Test that calculate_target_alignments correctly handles nested/overlapping
    alignments — the real-world bug where SUPER_19 was inflated 37x on chr_19."""

    def _make_record(self, query_name, query_len, qstart, qend, strand, target):
        return PAFRecord(
            query_name=query_name,
            query_length=query_len,
            query_start=qstart,
            query_end=qend,
            strand=strand,
            target_name=target,
            target_length=100_000_000,
            target_start=0,
            target_end=qend - qstart,
            num_matches=qend - qstart,
            alignment_length=qend - qstart,
            mapping_quality=60,
        )

    def test_identical_records_count_once(self):
        """20 identical alignments to the same region must count as one region."""
        records = [
            self._make_record("SUPER_19", 55_000_000, 40302, 110580, "+", "chr_19")
            for _ in range(20)
        ]
        stats = calculate_target_alignments(records)
        expected = 110580 - 40302  # 70278 bp
        assert stats["chr_19"]["total"] == expected

    def test_nested_alignments_not_double_counted(self):
        """Alignments where one interval contains another count unique bases only."""
        records = [
            self._make_record("SUPER_19", 55_000_000, 0, 100_000, "+", "chr_19"),   # outer
            self._make_record("SUPER_19", 55_000_000, 10_000, 50_000, "+", "chr_19"),  # inner
            self._make_record("SUPER_19", 55_000_000, 20_000, 80_000, "+", "chr_19"),  # inner
        ]
        stats = calculate_target_alignments(records)
        assert stats["chr_19"]["total"] == 100_000  # only the outer range

    def test_correct_winner_with_nested_inflated_loser(self):
        """chr_17 wins over chr_19 when chr_19 has nested duplicate alignments.

        Mirrors the production bug: SUPER_19 had 37x inflated chr_19 score,
        but chr_17 had genuine unique coverage > chr_19 actual unique coverage.
        """
        query_len = 55_000_000

        # chr_19: same 70 kbp region repeated 37 times → raw 2.6M, unique 70 kbp
        chr19_records = [
            self._make_record("SUPER_19", query_len, 40302, 110580, "-", "chr_19")
            for _ in range(37)
        ]

        # chr_17: real unique 53 Mbp coverage
        chr17_records = [
            self._make_record("SUPER_19", query_len, i * 1_000_000, (i + 1) * 1_000_000, "+", "chr_17")
            for i in range(53)
        ]

        stats = calculate_target_alignments(chr19_records + chr17_records)

        assert stats["chr_19"]["total"] == 110580 - 40302  # 70 278 bp — not 37x
        assert stats["chr_17"]["total"] == 53_000_000

        # chr_17 must win
        best = max(stats, key=lambda t: stats[t]["total"])
        assert best == "chr_17"

    def test_strand_totals_also_deduplicated(self):
        """Plus and minus strand totals are also based on unique intervals."""
        records = [
            self._make_record("SUPER_1", 1_000_000, 0, 50_000, "+", "chr_1"),
            self._make_record("SUPER_1", 1_000_000, 0, 50_000, "+", "chr_1"),  # duplicate +
            self._make_record("SUPER_1", 1_000_000, 60_000, 100_000, "-", "chr_1"),
            self._make_record("SUPER_1", 1_000_000, 60_000, 100_000, "-", "chr_1"),  # duplicate -
        ]
        stats = calculate_target_alignments(records)
        assert stats["chr_1"]["plus"] == 50_000
        assert stats["chr_1"]["minus"] == 40_000
        assert stats["chr_1"]["total"] == 90_000  # non-overlapping: 0-50k + 60k-100k


class TestUnlocOrderInOutput:
    """Test that unloc contigs appear immediately after their parent chromosome
    in both the FASTA output and the chromosome list CSV.

    Before the fix, write_fasta() wrote all chromosomes first and then all unloc
    contigs as a single block at the end (just before shrapnel).  The expected
    behaviour is:
        SUPER_1
        SUPER_1_unloc_1
        SUPER_1_unloc_2
        SUPER_2
        SUPER_2_unloc_1
        SUPER_3          <- no unloc
        <shrapnel>
    """

    # ------------------------------------------------------------------ helpers

    def _make_assignment(self, original, new_name, suffix, rc=False, sex=False):
        return FinalChromosomeAssignment(
            original_name=original,
            new_name=new_name,
            new_suffix=suffix,
            needs_reverse_complement=rc,
            is_sex_chromosome=sex,
        )

    def _make_unloc(self, contig, parent, num, rc=False):
        from rename_and_orient import UnlocMapping
        return UnlocMapping(
            contig_name=contig,
            parent_chromosome=parent,
            unloc_number=num,
            needs_reverse_complement=rc,
        )

    def _sequences(self):
        return {
            "SUPER_1":         "AAAA",
            "SUPER_1_unloc_1": "CCCC",
            "SUPER_1_unloc_2": "GGGG",
            "SUPER_2":         "TTTT",
            "SUPER_2_unloc_1": "ACAC",
            "SUPER_3":         "TGTG",
            "shrapnel_1":      "NNNN",
        }

    # ------------------------------------------------------------------ fasta

    def test_unloc_immediately_follows_parent_in_fasta(self, tmp_path):
        """Headers in output FASTA must follow parent → unloc(s) → next parent order."""
        assignments = sort_assignments_for_output([
            self._make_assignment("SUPER_1", "SUPER_1", "1"),
            self._make_assignment("SUPER_2", "SUPER_2", "2"),
            self._make_assignment("SUPER_3", "SUPER_3", "3"),
        ])
        unloc_mappings = [
            self._make_unloc("SUPER_1_unloc_1", "SUPER_1", 1),
            self._make_unloc("SUPER_1_unloc_2", "SUPER_1", 2),
            self._make_unloc("SUPER_2_unloc_1", "SUPER_2", 1),
        ]
        out_fasta = tmp_path / "out.fa"
        write_fasta(self._sequences(), assignments, unloc_mappings,
                    out_fasta, output_prefix="SUPER_")

        headers = [line.strip()[1:] for line in out_fasta.read_text().splitlines()
                   if line.startswith(">")]

        expected_order = [
            "SUPER_1",
            "SUPER_1_unloc_1",
            "SUPER_1_unloc_2",
            "SUPER_2",
            "SUPER_2_unloc_1",
            "SUPER_3",
            "shrapnel_1",
        ]
        assert headers == expected_order, (
            f"Unexpected FASTA order.\nGot:      {headers}\nExpected: {expected_order}"
        )

    def test_unloc_not_grouped_at_end_of_chromosomes(self, tmp_path):
        """Regression: unloc contigs must NOT appear after all chromosomes."""
        assignments = sort_assignments_for_output([
            self._make_assignment("SUPER_1", "SUPER_1", "1"),
            self._make_assignment("SUPER_2", "SUPER_2", "2"),
        ])
        unloc_mappings = [
            self._make_unloc("SUPER_1_unloc_1", "SUPER_1", 1),
        ]

        out_fasta = tmp_path / "out.fa"
        write_fasta(self._sequences(), assignments, unloc_mappings,
                    out_fasta, output_prefix="SUPER_")

        headers = [line.strip()[1:] for line in out_fasta.read_text().splitlines()
                   if line.startswith(">")]

        super2_idx = headers.index("SUPER_2")
        unloc_idx = headers.index("SUPER_1_unloc_1")

        assert unloc_idx < super2_idx, (
            "SUPER_1_unloc_1 must appear before SUPER_2, not after all chromosomes. "
            f"Got order: {headers}"
        )

    # ------------------------------------------------------------------ csv

    def test_unloc_immediately_follows_parent_in_csv(self, tmp_path):
        """Rows in chromosome list CSV must follow parent → unloc(s) → next parent."""
        assignments = sort_assignments_for_output([
            self._make_assignment("SUPER_1", "SUPER_1", "1"),
            self._make_assignment("SUPER_2", "SUPER_2", "2"),
            self._make_assignment("SUPER_3", "SUPER_3", "3"),
        ])
        unloc_mappings = [
            self._make_unloc("SUPER_1_unloc_1", "SUPER_1", 1),
            self._make_unloc("SUPER_1_unloc_2", "SUPER_1", 2),
            self._make_unloc("SUPER_2_unloc_1", "SUPER_2", 1),
        ]

        out_csv = tmp_path / "out.csv"
        write_chromosome_list(assignments, unloc_mappings, out_csv, output_prefix="SUPER_")

        names = [line.split(",")[0] for line in out_csv.read_text().splitlines() if line]

        expected_order = [
            "SUPER_1",
            "SUPER_1_unloc_1",
            "SUPER_1_unloc_2",
            "SUPER_2",
            "SUPER_2_unloc_1",
            "SUPER_3",
        ]
        assert names == expected_order, (
            f"Unexpected CSV order.\nGot:      {names}\nExpected: {expected_order}"
        )


class TestLoadMappingTable:
    """Tests for load_mapping_table_assignments (--mapping-table mode)."""

    HEADER = "query\ttarget\trenamed_to\tquery_length\talignment_length\tcoverage\tplus_strand\tminus_strand\tneeds_reverse_complement\n"

    def _write_table(self, tmp_path, rows: list[str]) -> Path:
        p = tmp_path / "mapping.tsv"
        p.write_text(self.HEADER + "".join(rows))
        return p

    # ------------------------------------------------------------------ basic

    def test_basic_rename_and_rc(self, tmp_path):
        """Chromosomes are renamed and RC flag applied from the table."""
        table = self._write_table(tmp_path, [
            "SUPER_9\tchr_1\tSUPER_1\t1000\t900\t0.9\t100\t800\tyes\n",
            "SUPER_3\tchr_2\tSUPER_2\t1000\t900\t0.9\t800\t100\tno\n",
        ])
        sequences = {"SUPER_9": "A" * 100, "SUPER_3": "T" * 100}
        assignments, unloc_mappings = load_mapping_table_assignments(
            table, sequences, query_chromosome_prefix="SUPER_", output_prefix="SUPER_"
        )
        by_orig = {a.original_name: a for a in assignments}
        assert by_orig["SUPER_9"].new_name == "SUPER_1"
        assert by_orig["SUPER_9"].needs_reverse_complement is True
        assert by_orig["SUPER_3"].new_name == "SUPER_2"
        assert by_orig["SUPER_3"].needs_reverse_complement is False
        assert unloc_mappings == []

    def test_sex_chromosome_detected(self, tmp_path):
        """Sex chromosome suffix is correctly detected from renamed_to."""
        table = self._write_table(tmp_path, [
            "SUPER_X\tchr_X\tSUPER_X\t500\t450\t0.9\t100\t350\tyes\n",
        ])
        sequences = {"SUPER_X": "N" * 50}
        assignments, _ = load_mapping_table_assignments(table, sequences)
        assert assignments[0].is_sex_chromosome is True
        assert assignments[0].new_suffix == "X"

    def test_unloc_in_fasta_mapped_via_parent(self, tmp_path):
        """Unloc contigs unique to haplotype 2 are mapped via their parent."""
        table = self._write_table(tmp_path, [
            "SUPER_5\tchr_1\tSUPER_1\t1000\t900\t0.9\t800\t100\tno\n",
        ])
        # haplotype 2 has an extra unloc not present in haplotype 1's table
        sequences = {"SUPER_5": "A" * 100, "SUPER_5_unloc_1": "T" * 50}
        assignments, unloc_mappings = load_mapping_table_assignments(
            table, sequences, query_chromosome_prefix="SUPER_", output_prefix="SUPER_"
        )
        assert len(unloc_mappings) == 1
        u = unloc_mappings[0]
        assert u.contig_name == "SUPER_5_unloc_1"
        assert u.parent_chromosome == "SUPER_5"
        assert u.unloc_number == 1
        assert u.needs_reverse_complement is False  # unlocs never RC'd

    def test_unloc_never_reverse_complemented(self, tmp_path):
        """Even when parent is RC'd, unloc itself must not be."""
        table = self._write_table(tmp_path, [
            "SUPER_2\tchr_1\tSUPER_1\t1000\t900\t0.9\t100\t800\tyes\n",
        ])
        sequences = {"SUPER_2": "A" * 100, "SUPER_2_unloc_1": "T" * 50}
        _, unloc_mappings = load_mapping_table_assignments(
            table, sequences, query_chromosome_prefix="SUPER_", output_prefix="SUPER_"
        )
        assert unloc_mappings[0].needs_reverse_complement is False

    def test_chr_not_in_table_kept_with_warning(self, tmp_path, capsys):
        """Chromosome missing from table keeps its original suffix, RC=False."""
        table = self._write_table(tmp_path, [
            "SUPER_1\tchr_1\tSUPER_1\t1000\t900\t0.9\t800\t100\tno\n",
        ])
        sequences = {"SUPER_1": "A" * 100, "SUPER_99": "C" * 100}
        assignments, _ = load_mapping_table_assignments(table, sequences)
        by_orig = {a.original_name: a for a in assignments}
        assert by_orig["SUPER_99"].new_name == "SUPER_99"
        assert by_orig["SUPER_99"].needs_reverse_complement is False
        assert "Warning" in capsys.readouterr().out

    def test_missing_column_raises(self, tmp_path):
        """ValueError raised when a required column is missing."""
        bad = tmp_path / "bad.tsv"
        bad.write_text("query\ttarget\tquery_length\n"
                       "SUPER_1\tchr_1\t1000\n")
        with pytest.raises(ValueError, match="missing columns"):
            load_mapping_table_assignments(bad, {"SUPER_1": "A"})

    def test_empty_table_raises(self, tmp_path):
        """ValueError raised when table has header but no data rows."""
        table = tmp_path / "empty.tsv"
        table.write_text(self.HEADER)
        with pytest.raises(ValueError, match="empty"):
            load_mapping_table_assignments(table, {"SUPER_1": "A"})

    def test_output_prefix_applied(self, tmp_path):
        """output_prefix is applied to new chromosome names."""
        table = self._write_table(tmp_path, [
            "SUPER_3\tchr_1\tSUPER_1\t1000\t900\t0.9\t800\t100\tno\n",
        ])
        sequences = {"SUPER_3": "A" * 100}
        assignments, _ = load_mapping_table_assignments(
            table, sequences, query_chromosome_prefix="SUPER_", output_prefix="chr_"
        )
        assert assignments[0].new_name == "chr_1"
