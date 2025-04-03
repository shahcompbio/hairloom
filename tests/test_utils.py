import pandas as pd
from pysam import AlignedSegment

import pytest

from hairloom.datatypes import BreakpointChain, Breakpoint
from hairloom.utils import (enumerate_breakpoints, get_secondaries, make_split_read_table, make_seg_table, make_brk_table, make_tra_table,
    is_breakpoints_not_sorted, reverse_complement, melt_tra_table, melt_brk_table)


class MockAlignment:
    """A mock class to simulate alignment objects for testing"""
    def __init__(self, read_name, refname, start, end, strand, clip1, match, clip2, pclip1):
        self.read_name = read_name
        self.refname = refname
        self.start = start
        self.end = end
        self.strand = strand
        self.clip1 = clip1
        self.match = match
        self.clip2 = clip2
        self.pclip1 = pclip1


class MockBreakpoint:
    """A mock class to simulate Breakpoint objects."""
    def __init__(self, chrom, pos, ori):
        self.chrom = chrom
        self.pos = pos
        self.ori = ori


class MockSegment:
    """A mock class to simulate Segment objects."""
    def __init__(self, brk1, brk2):
        self.brk1 = brk1
        self.brk2 = brk2


class MockBreakpointChain:
    """A mock class to simulate BreakpointChain objects."""
    def __init__(self, segs=None, tras=None):
        self.segs = segs if segs else []
        self.tras = tras if tras else []


def test_enumerate_breakpoints():
    df = pd.DataFrame({
        'chrom': {0: 'chr6', 1: 'chr1', 2: 'chr6'},
        'start': {0: 26424060, 1: 4596, 2: 152942012},
        'end': {0: 26424268, 1: 4995, 2: 152944889},
        'strand': {0: '+', 1: '-', 2: '+'},
        'clip1': {0: 30, 1: 2886, 2: 631},
        'match': {0: 209, 1: 398, 2: 2878},
        'clip2': {0: 3280, 1: 235, 2: 10},
        'pclip1': {0: 30, 1: 235, 2: 631}
    })
    brks = enumerate_breakpoints(df)
    assert len(brks) % 2 == 0, brks
    assert brks[0].pos == 26424268, brks[0]
    assert brks[0].ori == '+', brks[0]
    assert brks[1].pos == 4995, brks[1]
    assert brks[1].ori == '+', brks[1]
    assert brks[2].pos == 4596, brks[2]
    assert brks[2].ori == '-', brks[2]
    assert brks[3] == brks[-1], brks[3]
    assert brks[3].pos == 152942012, brks[3]
    assert brks[3].ori == '-', brks[3]


def test_get_secondaries():
    # Case 1: Read with a valid 'SA' tag
    read_with_sa = AlignedSegment()
    read_with_sa.tags = [('SA', 'chr1,100,+,60M,60;chr2,200,-,60M,60;')]
    expected_result = ['chr1,100,+,60M,60', 'chr2,200,-,60M,60']
    assert get_secondaries(read_with_sa) == expected_result

    # Case 2: Read with no 'SA' tag
    read_without_sa = AlignedSegment()
    read_without_sa.tags = [('NM', 1)]  # No 'SA' tag
    assert get_secondaries(read_without_sa) == []

    # Case 3: Read with an empty 'SA' tag
    read_with_empty_sa = AlignedSegment()
    read_with_empty_sa.tags = [('SA', '')]
    assert get_secondaries(read_with_empty_sa) == []

    # Case 4: Read with improperly formatted 'SA' tag (edge case)
    read_with_partial_sa = AlignedSegment()
    read_with_partial_sa.tags = [('SA', 'chr1,100,+,60M,60;;')]
    expected_result_partial = ['chr1,100,+,60M,60']
    assert get_secondaries(read_with_partial_sa) == expected_result_partial

    # Case 5: Read with multiple tags including 'SA'
    read_with_multiple_tags = AlignedSegment()
    read_with_multiple_tags.tags = [
        ('NM', 1),
        ('SA', 'chr3,300,+,50M,30;chr4,400,-,40M,40;')
    ]
    expected_result_multiple = ['chr3,300,+,50M,30', 'chr4,400,-,40M,40']
    assert get_secondaries(read_with_multiple_tags) == expected_result_multiple


def test_make_split_read_table():
    # Test case 1: Valid alignments
    alignments = [
        MockAlignment("read1", "chr1", 100, 150, "+", 10, 40, 5, 10),
        MockAlignment("read2", "chr2", 200, 250, "-", 15, 35, 10, 15),
        MockAlignment("read1", "chr1", 100, 150, "+", 10, 40, 5, 10),  # Duplicate entry
    ]
    expected_data = {
        "qname": ["read1", "read2"],
        "chrom": ["chr1", "chr2"],
        "start": [100, 200],
        "end": [150, 250],
        "strand": ["+", "-"],
        "clip1": [10, 15],
        "match": [40, 35],
        "clip2": [5, 10],
        "pclip1": [10, 15],
    }
    expected_df = pd.DataFrame(expected_data)
    result_df = make_split_read_table(alignments)

    pd.testing.assert_frame_equal(result_df, expected_df)

    # Test case 2: Empty alignments list
    alignments = []
    expected_df = pd.DataFrame(columns=["qname", "chrom", "start", "end", "strand", "clip1", "match", "clip2", "pclip1"])
    result_df = make_split_read_table(alignments)

    pd.testing.assert_frame_equal(result_df, expected_df)

    # Test case 3: Single alignment
    alignments = [
        MockAlignment("read1", "chr1", 300, 350, "+", 5, 45, 10, 5),
    ]
    expected_data = {
        "qname": ["read1"],
        "chrom": ["chr1"],
        "start": [300],
        "end": [350],
        "strand": ["+"],
        "clip1": [5],
        "match": [45],
        "clip2": [10],
        "pclip1": [5],
    }
    expected_df = pd.DataFrame(expected_data)
    result_df = make_split_read_table(alignments)

    pd.testing.assert_frame_equal(result_df, expected_df)


def test_make_seg_table():
    # Prepare test data
    brk1 = Breakpoint("chr1", 100, "+")
    brk2 = Breakpoint("chr1", 200, "+")
    brk3 = Breakpoint("chr1", 300, "-")
    brk4 = Breakpoint("chr1", 400, "-")
    bundle = [BreakpointChain([brk1, brk2, brk3, brk4])]

    # Test function
    result = make_seg_table(bundle)
    expected = pd.DataFrame({"chrom": ["chr1"], "pos1": [200], "pos2": [300], "support": [1]})
    pd.testing.assert_frame_equal(result, expected)


def test_make_brk_table():
    # Prepare test data
    brk1 = MockBreakpoint("chr1", 100, "+")
    brk2 = MockBreakpoint("chr1", 200, "-")
    bundle = [[brk1, brk2]]

    brk_supports = {("chr1", 100, "+"): 1, ("chr1", 200, "-"): 1}

    # Test function
    result = make_brk_table(bundle)
    expected = pd.DataFrame({
        "chrom": ["chr1", "chr1"],
        "pos": [100, 200],
        "ori": ["+", "-"],
        "support": [1, 1]
    })
    pd.testing.assert_frame_equal(result, expected)


def test_make_tra_table():
    # Prepare test data
    brk1 = Breakpoint("chr1", 100, "+")
    brk2 = Breakpoint("chr2", 200, "-")
    bundle = [BreakpointChain([brk1, brk2])]

    # Test function
    result = make_tra_table(bundle)
    expected = pd.DataFrame({
        "chrom1": ["chr1"],
        "pos1": [100],
        "ori1": ["+"],
        "chrom2": ["chr2"],
        "pos2": [200],
        "ori2": ["-"],
        "support": [1]
    })
    pd.testing.assert_frame_equal(result, expected)


def test_is_breakpoints_not_sorted():
    chrom_order = ["chr1", "chr2", "chr3", "chrX", "chrY"]

    # Test for unsorted chromosomes
    assert is_breakpoints_not_sorted("chr2", 100, "chr1", 200, chrom_order) is True

    # Test for sorted chromosomes
    assert is_breakpoints_not_sorted("chr1", 200, "chr2", 100, chrom_order) is False

    # Test for same chromosome, unsorted positions
    assert is_breakpoints_not_sorted("chr1", 300, "chr1", 200, chrom_order) is True

    # Test for same chromosome, sorted positions
    assert is_breakpoints_not_sorted("chr1", 100, "chr1", 200, chrom_order) is False

    # Test for edge case of equal chromosomes and positions
    assert is_breakpoints_not_sorted("chr1", 200, "chr1", 200, chrom_order) is False


def test_reverse_complement():
    # Test for a valid DNA sequence
    assert reverse_complement("ATCG") == "CGAT"

    # Test for a sequence with 'N' (undefined base)
    assert reverse_complement("ATGCN") == ".GCAT"

    # Test for an empty sequence
    assert reverse_complement("") == ""

    # Test for a palindrome sequence
    assert reverse_complement("ATGCAT") == "ATGCAT"

    # Test for case sensitivity
    with pytest.raises(KeyError):
        reverse_complement("atcg")  # Lowercase input should raise a KeyError


def test_melt_tra_table():
    df = pd.DataFrame([
        ['chr1', 10, '+', 'chr1', 100, '+', 10],
        ['chr2', 20, '+', 'chr2', 200, '-', 9],
    ], columns=['chrom1', 'pos1', 'ori1', 'chrom2', 'pos2', 'ori2', 'support'])
    expected = {
        Breakpoint('chr1', 10, '+'),
        Breakpoint('chr1', 100, '+'),
        Breakpoint('chr2', 20, '+'),
        Breakpoint('chr2', 200, '-'),
    }
    breakpoints = melt_tra_table(df)
    assert expected == breakpoints, breakpoints


def test_melt_brk_table():
    df = pd.DataFrame([
        ['chr1', 10, '+', 10],
        ['chr1', 100, '-', 5],
        ['chr2', 20, '+', 2],
    ], columns=['chrom', 'pos', 'ori', 'support'])
    expected = {
        Breakpoint('chr1', 10, '+'),
        Breakpoint('chr1', 100, '-'),
        Breakpoint('chr2', 20, '+'),
    }
    breakpoints = melt_brk_table(df)
    assert expected == breakpoints, breakpoints