import pytest
import importlib.resources as pkg_resources

import pyfaidx
import pandas as pd

from hairloom.datatypes import get_breakpoint_seqs
from hairloom.datatypes import (Breakpoint, BreakpointPair, BreakpointChain, Transitions, Segments)
from hairloom.utils import enumerate_breakpoints


def get_fasta_path():
    """Get the full path of the test_datatypes.fasta file."""
    return str(pkg_resources.files("hairloom").joinpath('data/test_datatypes.fasta'))


@pytest.fixture
def mock_genome():
    """Mock genome as a dictionary of sequences.
    ```
    $ cat tests/data/test_datatypes.fasta
    >chr1
    AACCGGTT
    >chr2
    ACACGTGT
    >6
    TCGATCGA
    >Circular
    GAAACCCT
    ```
    """
    fasta_path = get_fasta_path()
    return pyfaidx.Fasta(fasta_path)

@pytest.fixture
def read_table():
    return pd.DataFrame(
        [
            ["testread", "chr6", 2, 3, "+"],
            ["testread", "chr1", 3, 4, "-"],
            ["testread", "chr6", 4, 5, "+"],
            ["testread", "chr2", 5, 6, "+"],
        ],
        columns=["qname", "chrom", "start", "end", "strand"]
    )


def test_get_breakpoint_seqs(mock_genome):
    up, dn = get_breakpoint_seqs('chr1', 3, margin=2, genome=mock_genome)
    assert up == "AA", f"Upstream sequence mismatch: {up}"
    assert dn == "CC", f"Downstream sequence mismatch: {dn}"

    up, dn = get_breakpoint_seqs('chr2', 5, margin=3, genome=mock_genome)
    assert up == "CAC", f"Upstream sequence mismatch: {up}"
    assert dn == "GTG", f"Downstream sequence mismatch: {dn}"

    up, dn = get_breakpoint_seqs('Circular', 8, margin=3, genome=mock_genome)
    assert up == "CCC", f"Upstream sequence mismatch: {up}"
    assert dn == "TGA", f"Downstream sequence mismatch: {dn}"

    up, dn = get_breakpoint_seqs('Circular', 1, margin=3, genome=mock_genome)
    assert up == "CCT", f"Upstream sequence mismatch: {up}"
    assert dn == "GAA", f"Downstream sequence mismatch: {dn}"


class TestBreakpoint:
    """Test suite for the Breakpoint class."""
    def test_breakpoint_comparison(self):
        brk1 = Breakpoint('chr1', 1000, '+')
        brk2 = Breakpoint('chr1', 2000, '-')
        assert brk1 < brk2, (brk1, brk2)

        brk1 = Breakpoint('chr1',  1000, '+')
        brk2 = Breakpoint('chr10', 1000, '-')
        assert brk1 < brk2, (brk1, brk2)

        brk1 = Breakpoint('chr1', 1000, '+')
        brk2 = Breakpoint('chrX', 1000, '-')
        assert brk1 < brk2, (brk1, brk2)

        brk1 = Breakpoint('chr9',  1000, '+')
        brk2 = Breakpoint('chr10', 1000, '-')
        assert brk1 < brk2, (brk1, brk2)


    def test_get_breakpoint_seqs(self, mock_genome):
        brk = Breakpoint('chr1', 3, '-')
        
        brk.get_breakpoint_seqs(margin=2, genome=mock_genome)
        assert brk.upstream == "AA", f"Upstream sequence mismatch: {brk.upstream}"
        assert brk.downstream == "CC", f"Downstream sequence mismatch: {brk.downstream}"
        assert brk.seq_rearranged == "CC", f"Rearranged sequence mismatch: {brk.seq_rearranged}"
        assert brk.seq_removed == "AA", f"Removed sequence mismatch: {brk.seq_removed}"
        
        brk = Breakpoint('Circular', 8, '+')
        brk.get_breakpoint_seqs(margin=3, genome=mock_genome)
        assert brk.upstream == "CCC", f"Upstream sequence mismatch: {brk.upstream}"
        assert brk.downstream == "TGA", f"Downstream sequence mismatch: {brk.downstream}"
        assert brk.seq_rearranged == "CCC", f"Rearranged sequence mismatch: {brk.seq_rearranged}"
        assert brk.seq_removed == "TGA", f"Removed sequence mismatch: {brk.seq_removed}"


    def test_equality(self):
        bp1 = Breakpoint('chr1', 100, '+')
        bp2 = Breakpoint('chr1', 100, '+')
        bp3 = Breakpoint('chr1', 200, '+')
        bp4 = Breakpoint('chr2', 100, '+')
        bp5 = Breakpoint('chr1', 100, '-')

        # Equality checks
        assert bp1 == bp2  # Same attributes
        assert bp1 != bp3  # Different position
        assert bp1 != bp4  # Different chromosome
        assert bp1 != bp5  # Different strand
        assert bp1 != "not_a_breakpoint"  # Different type should not be equal

    def test_hashing(self):
        bp1 = Breakpoint('chr1', 100, '+')
        bp2 = Breakpoint('chr1', 100, '+')
        bp3 = Breakpoint('chr1', 200, '+')

        # Hash checks
        assert hash(bp1) == hash(bp2)  # Same attributes, same hash
        assert hash(bp1) != hash(bp3)  # Different attributes, different hash

    def test_set_inclusion(self):
        bp1 = Breakpoint('chr1', 100, '+')
        bp2 = Breakpoint('chr1', 100, '+')
        bp3 = Breakpoint('chr1', 200, '+')

        # Set behavior
        bp_set = {bp1, bp2, bp3}
        assert len(bp_set) == 2  # bp1 and bp2 are considered the same
        assert bp1 in bp_set
        assert bp3 in bp_set


class TestBreakpointChain:
    """Test suite for the BreakpointChain class."""
    def test_sort(self):
        brks = BreakpointChain([
            Breakpoint('chr2', 3, '+'),
            Breakpoint('chr1', 2, '+'),
            Breakpoint('chr1', 4, '+'),
            Breakpoint('chr1', 2, '+'),
        ])
        expected_tra0 = 'chr1:2:+-chr2:3:+'
        expected_tra1 = 'chr1:2:+-chr1:4:+'
        assert str(brks.tras[0]) == expected_tra0, f"Transition 0 mismatch: {brks.tras[0]}"
        assert str(brks.tras[1]) == expected_tra1, f"Transition 1 mismatch: {brks.tras[1]}"


    def test_length_validation(self):
        # Valid chain with a multiple of 2 length
        brk1 = Breakpoint("chr1", 100, "+")
        brk2 = Breakpoint("chr1", 200, "-")
        valid_chain = BreakpointChain([brk1, brk2])  # No error expected
        assert len(valid_chain) % 2 == 0, "BreakpointChain length must be a multiple of two."

        # Invalid chain with an odd length
        brk3 = Breakpoint("chr2", 300, "+")
        with pytest.raises(ValueError, match="BreakpointChain length must be a multiple of two."):
            BreakpointChain([brk1, brk2, brk3])  # Should raise an error


    def test_set_operations(self):
        # Create breakpoints
        brk1 = Breakpoint("chr1", 100, "+")
        brk2 = Breakpoint("chr1", 200, "-")
        brk3 = Breakpoint("chr2", 300, "+")
        brk4 = Breakpoint("chr2", 400, "-")
        brk5 = Breakpoint("chr3", 500, "+")

        # Create chains
        chain1 = BreakpointChain([brk1, brk2, brk3, brk4])
        chain2 = BreakpointChain([brk1, brk5, brk3, brk4])

        # Perform set operations
        intersection = set(chain1) & set(chain2)
        union = set(chain1) | set(chain2)
        difference_chain1 = set(chain1) - set(chain2)
        difference_chain2 = set(chain2) - set(chain1)

        # Assertions
        assert intersection == {brk1, brk3, brk4}, f"Expected {intersection} to be {{brk1, brk3, brk4}}"
        assert union == {brk1, brk2, brk3, brk4, brk5}, f"Expected {union} to be {{brk1, brk2, brk3, brk4, brk5}}"
        assert difference_chain1 == {brk2}, f"Expected {difference_chain1} to be {{brk2}}"
        assert difference_chain2 == {brk5}, f"Expected {difference_chain2} to be {{brk5}}"


def test_get_transitions(read_table):
    # Run enumerate_breakpoints on the dataframe
    brks = enumerate_breakpoints(read_table)
    
    # Check the number of transitions
    assert len(brks.tras) == 3, f"Transitions count mismatch: {len(brks.tras)} != 3"
    
    # Compare the actual transitions list with the expected list
    expected_transitions = [
        (('chr1', 4, '+'), ('chr6', 3, '+')), # sorted form
        (('chr1', 3, '-'), ('chr6', 4, '-')),
        (('chr2', 5, '-'), ('chr6', 5, '+')), # sorted form
    ]
    actual_transitions = [(tra.brk1.chrom, tra.brk1.pos, tra.brk1.ori,
                           tra.brk2.chrom, tra.brk2.pos, tra.brk2.ori) for tra in brks.tras]

    assert actual_transitions == [
        (et[0][0], et[0][1], et[0][2], et[1][0], et[1][1], et[1][2]) 
        for et in expected_transitions
    ], f"Transitions mismatch: {actual_transitions}"


def test_get_transition_list(read_table):
    transitions = Transitions(read_table)
    assert len(transitions.list) == 3, f"Transitions list mismatch: {transitions.list}"
    assert transitions.list == [
        (('chr6', 3, '+'), ('chr1', 4, '+')), # unsorted form
        (('chr1', 3, '-'), ('chr6', 4, '-')),
        (('chr6', 5, '+'), ('chr2', 5, '-')), # unsorted form
    ], f"Transitions list content mismatch: {transitions.list}"


def test_get_segments(read_table):
    brks = enumerate_breakpoints(read_table)

    # Validate the number of segments
    assert len(brks.segs) == 2, f"Segments mismatch: {brks.segs}"

    # Validate segment details
    assert brks.segs[0].brk1.pos == 3, f"Segment 0 start position mismatch: {brks.segs[0].brk1.pos}"
    assert brks.segs[0].brk2.pos == 4, f"Segment 0 end position mismatch: {brks.segs[0].brk2.pos}"
    assert brks.segs[1].brk1.pos == 4, f"Segment 1 start position mismatch: {brks.segs[1].brk1.pos}"
    assert brks.segs[1].brk2.pos == 5, f"Segment 1 end position mismatch: {brks.segs[1].brk2.pos}"


def test_get_segment_list(read_table):
    segments = Segments(read_table)

    # Validate the segments DataFrame shape
    assert segments.df.shape[0] == 4, f"Segments DataFrame mismatch: {segments.df}"

    # Ensure 'qname' is present in columns
    assert "qname" in segments.df.columns, "Missing 'qname' column in segments DataFrame"

    # Validate the segment list
    expected_segments = [
        ("chr1", 3, 4),
        ("chr6", 4, 5),
    ]
    assert segments.list == expected_segments, f"Segments content mismatch: {segments.list}"


class TestBreakpointPair:
    """Test suite for the BreakpointPair class."""

    def test_initialization(self):
        """Test BreakpointPair initialization."""
        brk1 = Breakpoint("chr1", 100, "+")
        brk2 = Breakpoint("chr2", 200, "-")

        pair = BreakpointPair(brk1, brk2)

        # Check attributes
        assert pair.brk1 == brk1, "Mismatch in brk1"
        assert pair.brk2 == brk2, "Mismatch in brk2"
        assert pair.aln_segment is False, "Mismatch in aln_segment default value"

    def test_repr(self):
        """Test the __repr__ method."""
        brk1 = Breakpoint("chr1", 100, "+")
        brk2 = Breakpoint("chr2", 200, "-")

        pair = BreakpointPair(brk1, brk2)
        expected_repr = "chr1:100:+-chr2:200:-"

        # Check string representation
        assert repr(pair) == expected_repr, f"String representation mismatch: {repr(pair)}"

    def test_attributes_are_independent(self):
        """Ensure that modifying the breakpoints does not affect the BreakpointPair instance."""
        brk1 = Breakpoint("chr1", 100, "+")
        brk2 = Breakpoint("chr2", 200, "-")
        pair = BreakpointPair(brk1, brk2)

        # Modify brk1 and brk2 attributes
        brk1.pos = 300
        brk2.chrom = "chr3"

        # Ensure changes to brk1 and brk2 do not affect the pair
        assert pair.brk1.pos == 300, "Mismatch in brk1 position"
        assert pair.brk2.chrom == "chr3", "Mismatch in brk2 chromosome"

    def test_equality(self):
        brk1 = Breakpoint("chr1", 100, "+")
        brk2 = Breakpoint("chr2", 200, "-")
        brk3 = Breakpoint("chr3", 300, "+")

        pair1 = BreakpointPair(brk1, brk2)
        pair2 = BreakpointPair(brk1, brk2)
        pair3 = BreakpointPair(brk2, brk1)  # Same pair but reverse order
        pair4 = BreakpointPair(brk1, brk3)

        # Equality checks
        assert pair1 == pair2  # Exact match
        assert pair1 == pair3  # Unordered equality
        assert pair1 != pair4  # Different breakpoints

    def test_hashing(self):
        brk1 = Breakpoint("chr1", 100, "+")
        brk2 = Breakpoint("chr2", 200, "-")
        brk3 = Breakpoint("chr3", 300, "+")

        pair1 = BreakpointPair(brk1, brk2)
        pair2 = BreakpointPair(brk1, brk2)
        pair3 = BreakpointPair(brk2, brk1)  # Same pair but reverse order
        pair4 = BreakpointPair(brk1, brk3)

        # Hash checks
        assert hash(pair1) == hash(pair2)  # Exact match
        assert hash(pair1) == hash(pair3)  # Unordered hash equality
        assert hash(pair1) != hash(pair4)  # Different pairs, different hash

    def test_set_inclusion(self):
        brk1 = Breakpoint("chr1", 100, "+")
        brk2 = Breakpoint("chr2", 200, "-")
        brk3 = Breakpoint("chr3", 300, "+")

        pair1 = BreakpointPair(brk1, brk2)
        pair2 = BreakpointPair(brk2, brk1)  # Same pair, reverse order
        pair3 = BreakpointPair(brk1, brk3)

        # Add to a set
        bp_set = {pair1, pair2, pair3}

        # Ensure only unique pairs are kept
        assert len(bp_set) == 2
        assert pair1 in bp_set
        assert pair3 in bp_set
