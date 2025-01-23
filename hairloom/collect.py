from collections import Counter

import tqdm
import pandas as pd
import pysam

from hairloom.utils import get_secondaries, make_split_read_table, enumerate_breakpoints
from hairloom.datatypes import SplitAlignment, Breakpoint, BreakpointPair

def extract_split_alignments(reads, max_reads=500) -> list:
    """Extract ``SplitAlignment`` objects from IteratorRow with a max_reads parameter

    Args:
        reads (pysam.IteratorRow): Reads fetched from a pysam.Alignmentfile
        max_reads (int, optional): Number of reads to extract at maximum. Defaults to 500.

    Returns:
        list: list of ``SplitAlignment`` objects
    """
    alignments = []
    for i, read in enumerate(reads):
        if i >= max_reads:
            break
        secondary_list = get_secondaries(read)
        if len(secondary_list) < 1: # only consider reads with secondaries (chimeric reads)
            continue # return None
        strand = '-' if read.is_reverse else '+'
        tpos = read.pos + 1
        cigarstring = read.cigarstring
        qname = read.qname
        refname = read.reference_name
        alignment = SplitAlignment(cigarstring, qname, refname, tpos, strand)
        alignments.append(alignment)

        for _, secondary in enumerate(secondary_list):
            refname_s, tpos_s, strand_s, cigarstring_s, _, _ = secondary.split(',')
            tpos_s = int(tpos_s)
            alignment = SplitAlignment(cigarstring_s, qname, refname_s, tpos_s, strand_s)
            alignments.append(alignment)
    return alignments

def extract_read_data(bam:pysam.AlignmentFile, contig:str, start=None, end=None, max_reads=500) -> pd.DataFrame:
    """
    Extract alignment tables for split reads and concatenate them into a single DataFrame.

    This function retrieves reads from a specified region in a BAM file, extracts split 
    alignments, and organizes them into a structured pandas DataFrame.

    Args:
        bam (pysam.AlignmentFile): The input BAM file opened with `pysam.AlignmentFile`.
        contig (str): The contig (e.g., chromosome) to extract reads from.
        start (int, optional): 1-based start position of the region. If None, reads 
            are fetched from the beginning of the contig. Defaults to None.
        end (int, optional): 1-based end position of the region. If None, reads 
            are fetched to the end of the contig. Defaults to None.
        max_reads (int, optional): The maximum number of reads to extract. Defaults to 500.

    Returns:
        pd.DataFrame: A DataFrame containing alignment data for all split reads in the region,
        concatenated and organized.

    Notes:
        - The `start` and `end` positions are converted to 0-based coordinates for 
          compatibility with `pysam.AlignmentFile.fetch`.
        - The `make_split_read_table` function is used to create the DataFrame from 
          the extracted alignments.

    Example:
        >>> import pysam
        >>> bam = pysam.AlignmentFile("example.bam", "rb")
        >>> df = extract_read_data(bam, contig="chr1", start=100, end=1000, max_reads=100)
        >>> print(df.head())

    """
    if start is not None:
        start -= 1
    reads = bam.fetch(contig=contig, start=start, end=end) # convert to 0-based pos
    alignments = extract_split_alignments(reads, max_reads=max_reads)
    df = make_split_read_table(alignments)
    return df
    

def pull_breakpoints_from_reads_in_sv_regions(bam, tra, get_read_table=False, min_n_breakpoint=2, margin=10) -> list:
    """Extract and append ``BreakpointChain`` objects from a bam file and a table of SVs

    Args:
        bam (pysam.AlignmentFile): BAM file
        tra (pandas.DataFrame): Table of SVs
        get_read_table (bool, optional): Return table of read alignment stats. Defaults to False.
        min_n_breakpoint (int, optional): Minimum number of breakpoints required to be saved. Useful in selecting complex rearrangements if the number is high. Defaults to 3.
        margin (int, optional): Margin (bp) from breakpoints to fetch reads. Defaults to 10.

    Returns:
        list: List of ``BreakpointChain``
    """
    bundle = []
    saved_qnames = set()
    canonical_chroms = [str(c) for c in range(1, 22+1)] + ['X', 'Y']
    canonical_chroms += ['chr'+c for c in canonical_chroms]
    read_info = pd.DataFrame()
    for _, row in tqdm.tqdm(tra.iterrows(), total=tra.shape[0]):
        # row = tra.iloc[6] # breakpoint pair: SV
        chrom1 = row['chromosome_1']
        pos1 = row['position_1']
        chrom2 = row['chromosome_2']
        pos2 = row['position_2']
        reads1 = bam.fetch(chrom1, pos1-margin-1, pos1+margin)
        reads2 = bam.fetch(chrom2, pos2-margin-1, pos2+margin)
        for reads in [reads1, reads2]:
            alignments = extract_split_alignments(reads)
            df = make_split_read_table(alignments)
            df = df[~df['qname'].isin(saved_qnames)] # don't double-count breakpoints already saved
            read_info = pd.concat([read_info, df])
            df = df[df['chrom'].isin(canonical_chroms)] # remove segments in decoys
            _qnames = set(df['qname'].unique().tolist())
            saved_qnames.update(_qnames)
            _bundle = make_bundle(df)
            for brks in _bundle:
                if len(brks) >= min_n_breakpoint:
                    bundle.append(brks)
    if get_read_table:
        return bundle, read_info
    return bundle

def make_bundle(reads_df):
    """Make a list of ``BreapointChain`` based on alignment table

    Args:
        reads_df (pandas.DataFrame): Table of read alignment statistics

    Returns:
        list: List of ``BreakpointChain``
    """
    bundle = []
    for qname, qdf in reads_df.groupby('qname'):
        brks = enumerate_breakpoints(qdf)
        brks.qname = qname
        brks.info = {'sv':brks.tras, 'qname':brks.qname}
        bundle.append(brks)
    return bundle

def get_breakpoint_support_from_bundle(bundle):
    """Get breakpoint support count

    Args:
        bundle (list): List of ``BreakpointChain``

    Returns:
        collections.Counter: Support for str(``Breakpoint``) coordinates
    """
    breakpoint_support = Counter()
    for brks in bundle:
        for brk in brks:
            breakpoint_support[str(brk)] += 1
    return breakpoint_support

def map_similar_coordinate_to_higher_rank(bundle, breakpoint_support, margin=10):
    """Make mapping of close-by coordinates, with breakpoints of higher support taking priority

    Args:
        bundle (list): List of ``BreakpointChain``
        breakpoint_support (dict | collections.Counter): Support for breakpoint coordinates
        margin (int, optional): Margin (bp) to merge close-by coordinates. Defaults to 10.

    Returns:
        tuple: tuple containing:
        
            coord_map (dict): src -> dst coordinate

            coord_map_log (tuple): (max_coord, src_count, max_count) [only for debugging]
    """
    coord_map = {}
    coord_map_log = {}
    for _, brks in enumerate(bundle):
        for _, brk in enumerate(brks):
            chrom, pos, ori = brk.chrom, brk.pos, brk.ori 
            src_coord = str(brk)
            src_count = breakpoint_support[src_coord]
            max_count = src_count
            max_coord = src_coord
            for _pos in range(pos-margin, pos+margin):
                _coord = f'{chrom}:{_pos}:{ori}'
                if _coord in breakpoint_support:
                    _count = breakpoint_support[_coord]
                    if _count > max_count:
                        max_count = _count
                        max_coord = _coord
            coord_map[src_coord] = max_coord
            coord_map_log[src_coord] = (max_coord, src_count, max_count)
    return coord_map, coord_map_log

def fix_lower_support_coordinates(bundle, coord_map):
    """Map breakpoint of lower support to close-by breakpoint with higher support

    Args:
        bundle (list): List of ``BreakpointChain``
        coord_map (dict): Map of str(``Breakpoint``) coordinates

    Returns:
        list: List of ``BreakpointChain``, mapped to fixed coordinates
    """
    for cix, brks in enumerate(bundle):
        for bix, brk in enumerate(brks):
            src_coord = str(brk)
            dst_coord = coord_map[src_coord]
            dst_chrom, dst_pos, dst_ori = dst_coord.split(':')
            dst_pos = int(dst_pos)
            bundle[cix][bix].pos = dst_pos
    return bundle

def normalize_sv_table(sv, chrom1_col='chromosome_1', chrom2_col='chromosome_2', 
                      pos1_col='position_1', pos2_col='position_2', 
                      ori1_col='strand_1', ori2_col='strand_2', chroms=None):
    """Sort breakpoint1 and breakpoint2 of a SV table

    Args:
        sv (pandas.DataFrame): Table of SVs
        chrom1_col (str, optional): Defaults to 'chromosome_1'.
        chrom2_col (str, optional): Defaults to 'chromosome_2'.
        pos1_col (str, optional): Defaults to 'position_1'.
        pos2_col (str, optional): Defaults to 'position_2'.
        ori1_col (str, optional): Defaults to 'strand_1'.
        ori2_col (str, optional): Defaults to 'strand_2'.
        chroms (list, optional): List of input contigs for coordinate sorting. Defaults to None.

    Returns:
        pandas.DataFrame: Sorted (normalized) SV table
    """
    sv = sv.copy()
    if chroms is None:
        chroms = [str(c) for c in range(1, 22+1)]+['X', 'Y']
        chroms += ['chr'+c for c in chroms] # chr prefices
    chrom_map = dict(zip(chroms, range(len(chroms))))
    assert (~sv[chrom1_col].isin(chroms)).sum() == 0, sv[chrom1_col].unique()
    assert (~sv[chrom2_col].isin(chroms)).sum() == 0, sv[chrom2_col].unique()
    flag_inverse = (sv[chrom1_col].map(chrom_map) > sv[chrom2_col].map(chrom_map))
    flag_inverse |= (sv[chrom1_col]==sv[chrom2_col]) & (sv[pos1_col] > sv[pos2_col])
    sv.loc[flag_inverse, [chrom1_col, chrom2_col]] = sv.loc[flag_inverse, [chrom2_col, chrom1_col]].values
    sv.loc[flag_inverse, [pos1_col, pos2_col]] = sv.loc[flag_inverse, [pos2_col, pos1_col]].values
    sv.loc[flag_inverse, [ori1_col, ori2_col]] = sv.loc[flag_inverse, [ori2_col, ori1_col]].values
    return sv

def get_svtype(tra:BreakpointPair):
    """Get SV type string for a given ``BreakpointPair``

    Args:
        tra (``BreakpointPair``): Paired breakpoint object

    Raises:
        ValueError: If no SV type has been assigned

    Returns:
        str: SV type string
    """
    translocation = 'TRA'
    inversion = 'INV'
    duplication = 'DUP'
    deletion = 'DEL'
    chrom1 = tra.brk1.chrom
    chrom2 = tra.brk2.chrom
    if chrom1 != chrom2:
        return translocation
    else: # same chrom
        pos1, pos2 = tra.brk1.pos, tra.brk2.pos
        ori1, ori2 = tra.brk1.ori, tra.brk2.ori
        if pos1 > pos2:
            pos1, pos2 = pos2, pos1
            ori1, ori2 = ori2, ori1
        if ori1 == ori2:
            return inversion
        else:
            if ori1 == '+':
                if ori2 == '-':
                    return deletion
            elif ori1 == '-':
                if ori2 == '+':
                    return duplication
    raise ValueError(f'tra:{tra}')

def pull_sv_supporting_reads_from_bundle(sv, bundle):
    """Filter bundle to include ``BreakpointChain`` objects that have breakpoints matching that of the input sv table

    Args:
        sv (pandas.DataFrame): SV table
        bundle (list): list of ``BreapointChain``

    Returns:
        list: Filtered list of ``BreakpointChain``
    """
    output_bundle = []
    sv_prop = {
        (sv['chromosome_1'], sv['position_1'], sv['strand_1']), 
        (sv['chromosome_2'], sv['position_2'], sv['strand_2']), 
    }
    for brks in bundle:
        brks.contains_sv_support = False
        for tra in brks.tras:
            tra.supports_sv = False
            brk1 = tra.brk1
            brk2 = tra.brk2
            tra_prop = {
                (brk1.chrom, brk1.pos, brk1.ori),
                (brk2.chrom, brk2.pos, brk2.ori),
            }
            if sv_prop == tra_prop:
                tra.supports_sv = True
                brks.contains_sv_support = True
        if brks.contains_sv_support:
            output_bundle.append(brks)
    return output_bundle

def find_presence_of_matching_sv(sv1, sv2, margin=50):
    """Check overlap of sv2 for sv1 table

    Args:
        sv1 (pandas.DataFrame): SV table to label matching SVs
        sv2 (pandas.DataFrame): SV table reference to check presence of overlap
        margin (int, optional): Margin (bp) of breakpoint coordinate difference. Defaults to 50.

    Returns:
        pd.Series: `{True, False}` list of matches. Length equal to sv1 row size.
    """
    match = []
    for i, (_, row) in enumerate(sv2.iterrows()):
        _match = (
            (sv1['chromosome_1'] == row['chromosome_1']) &
            (sv1['chromosome_2'] == row['chromosome_2']) &
            (abs(sv1['position_1']-row['position_1']) < margin) &
            (abs(sv1['position_2']-row['position_2']) < margin) &
            (sv1['strand_1'] == row['strand_1']) &
            (sv1['strand_2'] == row['strand_2'])
        )
        if i == 0:
            match = _match
        else:
            match |= _match
    return match

def make_tumor_sv_table(bundle:list, sv=None, margin=10, get_support=True) -> pd.DataFrame:
    """Make SV table from list of ``BreakpointChain``

    Args:
        bundle (list): List of ``BreakpointChain``
        sv (pandas.DataFrame, optional): Table of source SVs as reference for `in_source` flag. Defaults to None
        margin (int, optional): Margin (bp) for merging clustered breakpoints. Defaults to 10.
        get_support (bool, optional): Merge breakpoints with same coordinates and add count as `support`. Defaults to True.

    Returns:
        pandas.DataFrame: SV table from bundle [, with `in_source` labels] [, collapsed by coordinate with support counts]
    """
    breakpoint_support = get_breakpoint_support_from_bundle(bundle)
    coord_map, coord_map_log = map_similar_coordinate_to_higher_rank(bundle, breakpoint_support, margin=margin)
    bundle = fix_lower_support_coordinates(bundle, coord_map)
    # breakpoint_support_post = get_breakpoint_support_from_bundle(bundle)
    sv_cols = ['chromosome_1', 'position_1', 'strand_1', 'chromosome_2', 'position_2', 'strand_2', 'type']
    data = []
    for brks in bundle:
        for tra in brks.tras:
            svtype = get_svtype(tra)
            chrom1, chrom2 = tra.brk1.chrom, tra.brk2.chrom
            pos1, pos2 = tra.brk1.pos, tra.brk2.pos
            ori1, ori2 = tra.brk1.ori, tra.brk2.ori
            field = [chrom1, pos1, ori1, chrom2, pos2, ori2, svtype]
            data.append(field)
    tumor_sv = pd.DataFrame(data, columns=sv_cols)
    tumor_sv = normalize_sv_table(tumor_sv)
    ix_cols = ['chromosome_1', 'position_1', 'strand_1', 'chromosome_2', 'position_2', 'strand_2', 'type']
    if get_support:
        tumor_sv = tumor_sv.groupby(ix_cols).value_counts().reset_index()
        tumor_sv.rename(columns={'count':'support'}, inplace=True)
    if type(sv) == pd.DataFrame:
        tumor_sv['in_source'] = find_presence_of_matching_sv(tumor_sv, sv, margin=50)
    return tumor_sv

def get_normalized_sv(tra:BreakpointPair) -> tuple:
    """
    Sorts (normalizes) a ``BreakpointPair`` based on chromosome and position order.

    This function ensures a consistent ordering of breakpoints in a ``BreakpointPair``
    by sorting them based on chromosome precedence and genomic position.

    Args:
        tra (``BreakpointPair``): A pair of breakpoints to normalize. Each breakpoint 
            in the pair is expected to have the attributes `chrom`, `pos`, and `ori`.

    Returns:
        tuple: A flattened tuple of the normalized breakpoint coordinates in the format:
               [chrom1, pos1, ori1, chrom2, pos2, ori2].

    Notes:
        - Chromosomes are sorted based on their order in ``Breakpoint.chroms``.
        - If the chromosomes are the same, the breakpoints are sorted by position.
        - The orientation (`ori`) remains associated with its respective breakpoint.

    Example:
        >>> from your_module import Breakpoint, BreakpointPair, get_normalized_sv
        >>> brk1 = Breakpoint("chr2", 100, "+")
        >>> brk2 = Breakpoint("chr1", 200, "-")
        >>> pair = BreakpointPair(brk1, brk2)
        >>> get_normalized_sv(pair)
        ('chr1', 200, '-', 'chr2', 100, '+')
    """
    chrom1, pos1, ori1 = tra.brk1.chrom, tra.brk1.pos, tra.brk1.ori
    chrom2, pos2, ori2 = tra.brk2.chrom, tra.brk2.pos, tra.brk2.ori
    if Breakpoint.chroms.index(chrom1) < Breakpoint.chroms.index(chrom1):
        chrom1, pos1, ori1, chrom2, pos2, ori2 = chrom2, pos2, ori2, chrom1, pos1, ori1
    elif chrom1 == chrom2:
        if pos1 > pos2:
            chrom1, pos1, ori1, chrom2, pos2, ori2 = chrom2, pos2, ori2, chrom1, pos1, ori1
    return (chrom1, pos1, ori1, chrom2, pos2, ori2)

