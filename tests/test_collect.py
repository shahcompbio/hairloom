import os
import pytest
from collections import Counter
from unittest.mock import MagicMock
import importlib.resources as pkg_resources

import pysam
import numpy as np
import pandas as pd

from hairloom.datatypes import Breakpoint, BreakpointPair, BreakpointChain

from hairloom.collect import (map_similar_coordinate_to_higher_rank, fix_lower_support_coordinates, 
    get_breakpoint_support_from_bundle, normalize_sv_table, pull_sv_supporting_reads_from_bundle, find_presence_of_matching_sv,
    extract_read_data, make_bundle, get_svtype)

@pytest.fixture
def bundle():
    brk_chain1 = [
        Breakpoint('chr1', 1000, '+'),
        Breakpoint('chr2', 2000, '-'),
        Breakpoint('chr1', 1000, '-'),
        Breakpoint('chr3', 1000, '+'),
    ]
    brk_chain2 = [
        Breakpoint('chr1', 1001, '+'),
        Breakpoint('chr1', 1011, '+'),
    ]
    brk_chain3 = [
        Breakpoint('chr1', 1001, '+'),
        Breakpoint('chr1', 1012, '+'),
    ]
    bundle = [BreakpointChain(brk_chain1), BreakpointChain(brk_chain2), BreakpointChain(brk_chain3)]
    return bundle

@pytest.fixture
def mock_bam():
    # Create a mock object for `pysam.AlignmentFile`
    bam = MagicMock()

    # Mock the `.fetch()` method to yield mocked reads
    def fetch(contig=None, start=None, end=None):
        # Simulate reads based on parameters
        reads = [
            MagicMock(
                reference_name="chr1",
                pos=100,
                cigarstring="100M",
                is_reverse=False,
                qname="read1",
                get_tag=MagicMock(return_value="NM:i:1")
            ),
            MagicMock(
                reference_name="chr1",
                pos=200,
                cigarstring="50M50S",
                is_reverse=True,
                qname="read2",
                get_tag=MagicMock(return_value="NM:i:2")
            ),
        ]
        for read in reads:
            yield read

    bam.fetch.side_effect = fetch
    return bam

@pytest.fixture
def mock_enumerate_breakpoints(monkeypatch):
    # Mock `enumerate_breakpoints` to return predefined `BreakpointChain`
    def mock_function(df):
        breakpoints = [
            Breakpoint("chr1", 100, "+"),
            Breakpoint("chr1", 200, "-"),
        ]
        return BreakpointChain(breakpoints)

    # Patch the original `enumerate_breakpoints` function
    monkeypatch.setattr("hairloom.collect.enumerate_breakpoints", mock_function)
    return mock_function

# def test_extract_read_data_with_start_and_end(mock_bam):
#     # Mock reads with secondaries
#     mock_reads = [
#         MagicMock(reference_name="chr1", pos=100, cigarstring="100M", is_reverse=False, qname="read1"),
#         MagicMock(reference_name="chr1", pos=200, cigarstring="50M50S", is_reverse=True, qname="read2")
#     ]
#     mock_bam.fetch.return_value = iter(mock_reads)

#     # Call the function
#     df = extract_read_data(mock_bam, contig="chr1", start=100, end=200, max_reads=10)

#     # Validate the DataFrame
#     assert isinstance(df, pd.DataFrame)
#     assert len(df) == len(mock_reads)

# def test_extract_read_data_no_start_or_end(mock_bam):
#     mock_bam.fetch.return_value = iter([])
#     df = extract_read_data(mock_bam, contig="chr1", start=None, end=None)
#     assert df.empty

def test_make_bundle_empty():
    reads = pd.DataFrame(columns=["qname", "alignment"])
    bundle = make_bundle(reads)
    assert bundle == []

def test_make_bundle_with_data(mock_enumerate_breakpoints):
    reads = pd.DataFrame({
        "qname": ["read1", "read1", "read2"],
        "alignment": ["data1", "data2", "data3"]
    })
    mock_enumerate_breakpoints.return_value = MagicMock(tras=[], qname=None, info={})
    bundle = make_bundle(reads)
    assert len(bundle) == 2

def test_get_breakpoint_support_from_bundle(bundle):
    expected_support = Counter({
        "chr1:1000:+": 1,
        "chr2:2000:-": 1,
        "chr1:1000:-": 1,
        "chr3:1000:+": 1,
        "chr1:1001:+": 2,
        "chr1:1011:+": 1,
        "chr1:1012:+": 1,
    })
    result = get_breakpoint_support_from_bundle(bundle)
    assert result == expected_support

def test_map_similar_coordinate_to_higher_rank(bundle):
    margin = 10
    print(f"Called with margin={margin}, bundle size={len(bundle)}")
    breakpoint_support = get_breakpoint_support_from_bundle(bundle)
    coord_map, coord_map_log = map_similar_coordinate_to_higher_rank(bundle, breakpoint_support, margin=margin)
    assert coord_map['chr1:1000:+'] == 'chr1:1001:+', coord_map
    assert coord_map['chr1:1000:-'] == 'chr1:1000:-', coord_map
    assert coord_map['chr1:1011:+'] == 'chr1:1001:+', coord_map
    assert coord_map['chr1:1012:+'] == 'chr1:1012:+', coord_map

# def test_map_similar_coordinate_to_higher_rank():
#     bundle = [
#         [Breakpoint("chr1", 100, "+"), Breakpoint("chr1", 110, "+")],
#         [Breakpoint("chr1", 200, "-")]
#     ]
#     breakpoint_support = Counter({
#         "chr1:100:+": 5,
#         "chr1:110:+": 3,
#         "chr1:200:-": 1
#     })
#     coord_map, coord_map_log = map_similar_coordinate_to_higher_rank(bundle, breakpoint_support, margin=10)
#     assert coord_map["chr1:110:+"] == "chr1:100:+"
#     assert coord_map["chr1:200:-"] == "chr1:200:-"

def test_fix_lower_support_coordinates(bundle):
    breakpoint_support = get_breakpoint_support_from_bundle(bundle)
    coord_map, coord_map_log = map_similar_coordinate_to_higher_rank(bundle, breakpoint_support, margin=10)
    fixed_bundle = fix_lower_support_coordinates(bundle, coord_map)
    expected_output1 = BreakpointPair(Breakpoint('chr1', 1001, '+'), Breakpoint('chr2', 2000, '-')) # 1000 fixed to 1001
    expected_output2 = BreakpointPair(Breakpoint('chr1', 1001, '+'), Breakpoint('chr1', 1001, '+')) # 1011 fixed to 1001
    assert str(fixed_bundle[0].tras[0]) == str(expected_output1), fixed_bundle[0].tras[0]
    assert str(fixed_bundle[1].tras[0]) == str(expected_output2), fixed_bundle[1].tras[0]

def test_normalize_sv_table():
    input_sv = pd.DataFrame({
    0: {
        'chromosome_1': '11',
        'position_1': 2000,
        'strand_1': '+',
        'chromosome_2': '2',
        'position_2': 1000,
        'strand_2': '-'
    },
    1: {
        'chromosome_1': '1',
        'position_1': 2000,
        'strand_1': '+',
        'chromosome_2': '1',
        'position_2': 1000,
        'strand_2': '-'
    }}).T
    expected_sv = pd.DataFrame({
    0: {
        'chromosome_1': '2',
        'position_1': 1000,
        'strand_1': '-',
        'chromosome_2': '11',
        'position_2': 2000,
        'strand_2': '+'
    },
    1: {
        'chromosome_1': '1',
        'position_1': 1000,
        'strand_1': '-',
        'chromosome_2': '1',
        'position_2': 2000,
        'strand_2': '+'
    }}).T
    output_sv = normalize_sv_table(input_sv)
    assert (output_sv != expected_sv).sum().sum() == 0, output_sv

def test_pull_sv_supporting_reads_from_bundle(bundle):
    sv = {'chromosome_1':'chr1', 'position_1':1000, 'strand_1':'+',
          'chromosome_2':'chr2', 'position_2':2000, 'strand_2':'-',}
    sv_reads = pull_sv_supporting_reads_from_bundle(sv, bundle)
    expected_output = [bundle[0]]
    assert sv_reads == expected_output, (expected_output, sv_reads)

def test_find_presence_of_matching_sv():
    ix_cols = ['chromosome_1', 'position_1', 'strand_1', 'chromosome_2', 'position_2', 'strand_2']
    sv1 = pd.DataFrame([
        ['chr1', 1000, '+', 'chr2', 2000, '+'], 
        ['chr3', 3000, '+', 'chr4', 4000, '+'], 
    ], columns=ix_cols)
    sv2 = pd.DataFrame([
        ['chr1', 1049, '+', 'chr2', 2049, '+'], 
        ['chr3', 3000, '+', 'chr4', 4000, '-'], 
    ], columns=ix_cols)
    expected = pd.DataFrame([
        ['chr1', 1000, '+', 'chr2', 2000, '+', True], 
        ['chr3', 3000, '+', 'chr4', 4000, '+', False], 
    ], columns=ix_cols+['match'])
    sv1['match'] = find_presence_of_matching_sv(sv1, sv2, margin=50)
    assert (sv1 != expected).sum().sum() == 0, sv1

def test_extract_read_data():
    bam_path = str(pkg_resources.files("hairloom").joinpath('data/test.bam'))
    assert os.path.exists(bam_path), f'{bam_path} does not exist.'
    bam = pysam.AlignmentFile(bam_path)
    df = extract_read_data(bam, contig='PBEF1NeoTransposon', start=1, end=2)
    assert df.shape[0] == 0, df

    bam_path = str(pkg_resources.files("hairloom").joinpath('data/test.bam'))
    assert os.path.exists(bam_path), f'{bam_path} does not exist.'
    bam = pysam.AlignmentFile(bam_path)
    df = extract_read_data(bam, contig='PBEF1NeoTransposon', start=1470, end=1477)
    assert df.shape[0] == 0, df

    expected = np.array([['02ce28f5-83e5-53a4-a7ed-96f331c6b305', 'chr10', 51339923,
        51340887, '+', 35, 960, 3763, 35],
       ['02ce28f5-83e5-53a4-a7ed-96f331c6b305', 'PBEF1NeoTransposon',
        1478, 4996, '-', 288, 3480, 990, 990],
       ['02ce28f5-83e5-53a4-a7ed-96f331c6b305', 'chr10', 51340883,
        51341166, '+', 4466, 281, 11, 4466]])
    bam_path = str(pkg_resources.files("hairloom").joinpath('data/test.bam'))
    assert os.path.exists(bam_path), f'{bam_path} does not exist.'
    bam = pysam.AlignmentFile(bam_path)
    df = extract_read_data(bam, contig='PBEF1NeoTransposon', start=1477, end=1478)
    assert np.all(df.to_numpy().astype(str) == expected.astype(str))

def test_get_svtype_translocation():
    # Different chromosomes
    brk1 = Breakpoint("chr1", 100, "+")
    brk2 = Breakpoint("chr2", 200, "-")
    tra = BreakpointPair(brk1, brk2)
    assert get_svtype(tra) == "TRA"

def test_get_svtype_inversion():
    # Same chromosome, same orientation
    brk1 = Breakpoint("chr1", 100, "+")
    brk2 = Breakpoint("chr1", 200, "+")
    tra = BreakpointPair(brk1, brk2)
    assert get_svtype(tra) == "INV"

    brk1 = Breakpoint("chr1", 200, "-")
    brk2 = Breakpoint("chr1", 100, "-")
    tra = BreakpointPair(brk1, brk2)
    assert get_svtype(tra) == "INV"

def test_get_svtype_deletion():
    # Same chromosome, '+-' orientation
    brk1 = Breakpoint("chr1", 100, "+")
    brk2 = Breakpoint("chr1", 200, "-")
    tra = BreakpointPair(brk1, brk2)
    assert get_svtype(tra) == "DEL"

    tra = BreakpointPair(brk2, brk1)
    assert get_svtype(tra) == "DEL"

def test_get_svtype_duplication():
    # Same chromosome, '-+' orientation
    brk1 = Breakpoint("chr1", 100, "-")
    brk2 = Breakpoint("chr1", 200, "+")
    tra = BreakpointPair(brk1, brk2)
    assert get_svtype(tra) == "DUP"

    tra = BreakpointPair(brk2, brk1)
    assert get_svtype(tra) == "DUP"

def test_get_svtype_invalid_orientation():
    # Same chromosome but unsupported orientation
    brk1 = Breakpoint("chr1", 100, "+")
    brk2 = Breakpoint("chr1", 200, "*")  # Invalid orientation
    tra = BreakpointPair(brk1, brk2)
    with pytest.raises(ValueError):
        get_svtype(tra)